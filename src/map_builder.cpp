/** Defines a class with a local occupancy grid
*/

#include <cmath>
#include <fstream>
#include <math.h>  // for round(), std::round() is since C++11.

#include <map_builder.h>
#include <chrono>
namespace local_map
{

const double g_default_p_occupied_when_laser = 0.9;
const double g_default_p_occupied_when_no_laser = 0.3;
const double g_default_large_log_odds = 100;
const double g_default_max_log_odds_for_belief = 20;



MapBuilder::MapBuilder(int width, int height, double resolution) :
  angle_resolution_(M_PI / 560),
  p_occupied_when_laser_(g_default_p_occupied_when_laser),
  p_occupied_when_no_laser_(g_default_p_occupied_when_no_laser),
  large_log_odds_(g_default_large_log_odds),
  max_log_odds_for_belief_(g_default_max_log_odds_for_belief),
  has_frame_id_(false)
{
  map_frame_id_ = "/local_map";
  map_.header.frame_id = map_frame_id_;
  map_.info.width = width;
  map_.info.height = height;
  map_.info.resolution = resolution;
  map_.info.origin.position.x = -static_cast<double>(width) / 2 * resolution;
  map_.info.origin.position.y = -static_cast<double>(height) / 2 * resolution;
  map_.info.origin.orientation.w = 1.0;
  map_.data.assign(width * height, -1);  // Fill with "unknown" occupancy.
  // log_odds = log(occupancy / (1 - occupancy); prefill with
  // occupancy = 0.5, equiprobability between occupied and free.
  log_odds_.assign(width * height, 0);

  map_width_ = width;
  map_height_ = height;
  resolution_ = resolution;

  ros::NodeHandle private_nh("~");
  private_nh.getParam("angle_resolution", angle_resolution_);
  private_nh.getParam("p_occupied_when_laser", p_occupied_when_laser_);
  if (p_occupied_when_laser_ <=0 || p_occupied_when_laser_ >= 1)
  {
    ROS_ERROR_STREAM("Parameter "<< private_nh.getNamespace() <<
        "/p_occupied_when_laser must be within ]0, 1[, setting to default (" <<
        g_default_p_occupied_when_laser << ")");
    p_occupied_when_laser_ = g_default_p_occupied_when_laser;
  }
  private_nh.getParam("p_occupied_when_no_laser", p_occupied_when_no_laser_);
  if (p_occupied_when_no_laser_ <=0 || p_occupied_when_no_laser_ >= 1)
  {
    ROS_ERROR_STREAM("Parameter "<< private_nh.getNamespace() <<
        "/p_occupied_when_no_laser must be within ]0, 1[, setting to default (" <<
        g_default_p_occupied_when_no_laser << ")");
    p_occupied_when_no_laser_ = g_default_p_occupied_when_no_laser;
  }
  private_nh.getParam("large_log_odds", large_log_odds_);
  if (large_log_odds_ <=0)
  {
    ROS_ERROR_STREAM("Parameter "<< private_nh.getNamespace() << "/large_log_odds must be positive, setting to default (" <<
        g_default_large_log_odds << ")");
    large_log_odds_ = g_default_large_log_odds;
  }
  private_nh.getParam("max_log_odds_for_belief", max_log_odds_for_belief_);
  try
  {
    std::exp(max_log_odds_for_belief_);
  }
  catch (std::exception)
  {
    ROS_ERROR_STREAM("Parameter "<< private_nh.getNamespace() << "/max_log_odds_for_belief too large, setting to default (" <<
        g_default_max_log_odds_for_belief << ")");
    max_log_odds_for_belief_ = g_default_max_log_odds_for_belief;
  }


  // Fill in the lookup cache.
  const double angle_start = -M_PI;
  const double angle_end = angle_start + 2 * M_PI - 1e-6;
  for (double a = angle_start; a <= angle_end; a += angle_resolution_)
  {
    ray_caster_.getRayCastToMapBorder(a, height, width, 0.9 * angle_resolution_);
  }

  first_scan_ = true;

  map_to_laser_.setIdentity();

  correlative_scan_matcher_ = new correlative_scan_math::correlativeScanMatcher(map_width_, map_height_, resolution);

  float linear_seacher_window = 0.3;
  float linear_seacher_step = 0.05;
  float angular_seacher_window = 0.3;
  float angular_seacher_step = 0.03;
  int max_depth = 0;
  private_nh.getParam("linear_seacher_window", linear_seacher_window);
  private_nh.getParam("linear_seacher_step", linear_seacher_step);
  private_nh.getParam("angular_seacher_window", angular_seacher_window);
  private_nh.getParam("angular_seacher_step", angular_seacher_step);
  private_nh.getParam("max_depth", max_depth);
  correlative_scan_matcher_->setSearchParameters(linear_seacher_window, 
                                                linear_seacher_step, 
                                                angular_seacher_window, 
                                                angular_seacher_step,
                                                max_depth);

}

/** Update occupancy and log odds for a point
 *
 * @param[in] occupied true if the point was measured as occupied
 * @param[in] idx pixel index
 * @param[in] ncol map width
 * @param[in,out] occupancy occupancy map to update
 * @param[in,out] log_odds log odds to update
 */
void MapBuilder::updatePointOccupancy(bool occupied, size_t idx, vector<int8_t>& occupancy, vector<double>& log_odds) const
{
  if (idx >= occupancy.size())
  {
    return;
  }

  if (occupancy.size() != log_odds.size())
  {
    ROS_ERROR("occupancy and count do not have the same number of elements");
    return;
  }

  // Update log_odds.
  double p;  // Probability of being occupied knowing current measurement.
  if (occupied)
  {
    p = p_occupied_when_laser_;
  }
  else
  {
    p = p_occupied_when_no_laser_;
  }
  // Original formula: Table 4.2, "Probabilistics robotics", Thrun et al., 2005:
  // log_odds[idx] = log_odds[idx] +
  //     std::log(p * (1 - p_occupancy) / (1 - p) / p_occupancy);
  // With p_occupancy = 0.5, this simplifies to:
  log_odds[idx] += std::log(p / (1 - p));
  if (log_odds[idx] < -large_log_odds_)
  {
    log_odds[idx] = -large_log_odds_;
  }
  else if(log_odds[idx] > large_log_odds_)
  {
    log_odds[idx] = large_log_odds_;
  }
  // Update occupancy.
  if (log_odds[idx] < -max_log_odds_for_belief_)
  {
    occupancy[idx] = 0;
  }
  else if (log_odds[idx] > max_log_odds_for_belief_)
  {
    occupancy[idx] = 100;
  }
  else
  {
    occupancy[idx] = static_cast<int8_t>(lround((1 - 1 / (1 + std::exp(log_odds[idx]))) * 100));
  }
}

/** Callback for the LaserScan subscriber.
 *
 * Update (geometrical transformation + probability update) the map with the current scan
 */
void MapBuilder::grow(const sensor_msgs::LaserScan& scan)
{
    if(first_scan_)
    {
        updateMap(scan, 0, 0, 0);
        correlative_scan_matcher_->updateMapLookUpTable(log_odds_);


        sensor_msgs::PointCloud2 cloud_scan;
        projector_.projectLaser(scan, cloud_scan);

        // convert from pointcloud2 to pcl
        pcl::PointCloud<pcl::PointXYZ>::Ptr pcl_cloud (new pcl::PointCloud<pcl::PointXYZ>);
        pcl::fromROSMsg(cloud_scan, *pcl_cloud);  // PointCloud2 to pcl::PointCloud


        tf::Transform current_tf;
        current_tf.setOrigin(tf::Vector3(1, 1, 0.0));
        tf::Quaternion rotation;
        rotation.setRPY(0, 0, 0.3);
        current_tf.setRotation(rotation);

        pcl::PointCloud<pcl::PointXYZ> new_cloud;
        pcl_ros::transformPointCloud(*pcl_cloud, new_cloud, current_tf);



        double x,y,theta;

        correlative_scan_matcher_->multiResolutionSearch(new_cloud, x, y, theta);

        first_scan_ = false;
    }
    return;


}

bool MapBuilder::updateMap(const sensor_msgs::LaserScan& scan, long int dx, long int dy, double theta)
{
    const bool has_moved = (dx != 0 || dy != 0);
    const int ncol = map_.info.width;


    // Update occupancy.
    const double xcenter = (map_width_ / 2) * resolution_;
    const double ycenter = (map_height_ / 2) * resolution_;
    int counter = 0;
    for (size_t i = 0; i < scan.ranges.size(); ++i)
    {
        if(scan.ranges[i] < 0.01 || scan.ranges[i] > 30)
        {
          continue;
        }
        const double angle = angles::normalize_angle(scan.angle_min + i * scan.angle_increment + theta);

        double x1 = scan.ranges[i] * std::cos(angle) + xcenter; // Can be negative
        double y1 = scan.ranges[i] * std::sin(angle) + ycenter;
        const long int xmap = lround(x1 / resolution_);
        const long int ymap = lround(y1 / resolution_);
        size_t point_index = (ymap * map_width_) + xmap;
        log_odds_[point_index] = 1;
        counter++;
        //   vector<size_t> pts;
        //     const bool obstacle_in_map = getRayCastToObstacle(map_, angle, scan.ranges[i], pts);
        //   if (pts.empty())
        //   {
        //     continue;
        //   }
        //   if (obstacle_in_map)
        //   {
        //     // The last point is the point with obstacle.
        //     const size_t last_pt = pts.back();
        //     updatePointOccupancy(true, last_pt, map_.data, log_odds_);
        //     pts.pop_back();
        //     count++;
        //   }
        //   // The remaining points are in free space.
        //   updatePointsOccupancy(false, pts, map_.data, log_odds_);
      }
      std::cout<<"total num:"<<counter<<std::endl;

      return has_moved;
}

/** Return the pixel list by ray casting from map origin to map border, first obstacle.
 *
 * Return the pixel list by ray casting from map origin to map border or first obstacle, whichever comes first.
 *
 * @param[in] map occupancy grid
 * @param[in] angle laser beam angle
 * @param[in] range laser beam range
 * @param[out] raycast list of pixel indexes touched by the laser beam
 * @return true if the last point of the pixel list is an obstacle (end of laser beam). 
 */
bool MapBuilder::getRayCastToObstacle(const nav_msgs::OccupancyGrid& map, double angle, double range, vector<size_t>& raycast)
{
  // Do not consider a 0-length range.
  if (range < 1e-10)
  {
    raycast.clear();
    return false;
  }

  const vector<size_t>& ray_to_map_border = ray_caster_.getRayCastToMapBorder(angle,
      map.info.height, map.info.width,  angle_resolution_);
  // range in pixel length. The ray length in pixels corresponds to the number
  // of pixels in the bresenham algorithm.
  const size_t pixel_range = lround(range * max(abs(std::cos(angle)), abs(std::sin(angle))) / map.info.resolution);
  size_t raycast_size;
  bool obstacle_in_map = pixel_range < ray_to_map_border.size();
  if (obstacle_in_map)
  {
    raycast_size = pixel_range;
  }
  else
  {
    raycast_size = ray_to_map_border.size();
  }
  raycast.clear();
  raycast.reserve(raycast_size);
  for (size_t i = 0; i < raycast_size; ++i)
  {
    raycast.push_back(ray_to_map_border[i]);
  }

  return obstacle_in_map;
}



} // namespace local_map

