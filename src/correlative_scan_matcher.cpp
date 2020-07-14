/** Defines a class with a local occupancy grid
*/

#include <cmath>
#include <fstream>
#include <math.h>  // for round(), std::round() is since C++11.

#include <correlative_scan_matcher.h>
#include <chrono>

namespace correlative_scan_math
{

correlativeScanMatcher::correlativeScanMatcher(int map_width, int map_height, float resolution)
{
  //ros::NodeHandle private_nh("~");
  map_width_ = map_width;
  map_height_ = map_height;
  resolution_ = resolution;
}

void correlativeScanMatcher::bruteForceSearch(const sensor_msgs::LaserScan& scan, double& x, double& y, double& theta)
{
  vector<Candidate> result_candidates;
  vector<Candidate> rotation_canditate_sets = generateRotationScan(scan, 0, 0, 0);

  int step_size = linear_search_window_ /  linear_search_step_;

  for(int i = -step_size; i < step_size; i++)
  {
    for (int j = -step_size; j < step_size; ++j)
    {
      vector<Candidate>::iterator it = rotation_canditate_sets.begin();
      for( ;it != rotation_canditate_sets.end(); ++it)
      {
        Candidate rotation_candidate = *it;
        Candidate new_candidate = rotation_candidate;
        new_candidate.x_offset = j;
        new_candidate.y_offset = i;
        new_candidate.score = LEAST_SCORE_NUMBER;
        scoreCandidate(new_candidate);
        result_candidates.push_back(new_candidate);
      }
    }
  }

  //const Candidate& best_candidate = *std::max(result_candidates.begin(), result_candidates.end());
  Candidate best_candidate;
  double max_score = LEAST_SCORE_NUMBER;
  for(int i = 0; i < result_candidates.size(); i++)
  {
    if(result_candidates[i].score > max_score)
    {
      max_score = result_candidates[i].score;
      best_candidate.x_offset = result_candidates[i].x_offset;
      best_candidate.y_offset = result_candidates[i].y_offset;
      best_candidate.orientation = result_candidates[i].orientation;
      best_candidate.score = result_candidates[i].score;
    }
  }
    
  x = best_candidate.x_offset * resolution_;
  y = best_candidate.y_offset * resolution_;
  theta = best_candidate.orientation;

  //std::cout<<"x: "<<best_candidate.x_offset<<"  y:"<<best_candidate.y_offset<<" theta:"<<best_candidate.orientation<<" score:"<<best_candidate.score<<std::endl;
  //std::cout<<"x: "<<x<<"  y:"<<y<<" theta:"<<theta<<" score:"<<best_candidate.score<<std::endl;
}


void correlativeScanMatcher::scoreCandidate(Candidate& candidate)
{
  double score = LEAST_SCORE_NUMBER;

  for(int i = 0; i < candidate.discretize_scan.size(); i++)
  {
    //TODO: Condider boundry!
    int offset = getPointIndex(candidate.x_offset, candidate.y_offset);
    size_t point_index = candidate.discretize_scan[i] + (size_t)offset;

    if(!pointInMap(point_index))
    {
      std::cout<<"Point out of boundry."<<std::endl;
      continue;
    }
    score += lookup_table_[point_index];
   
  }

  candidate.score = score;
}

int correlativeScanMatcher::getPointIndex(int x, int y)
{
  return  x + y * map_width_;
}

vector<Candidate> correlativeScanMatcher::generateRotationScan(const sensor_msgs::LaserScan& scan, float x, float y, float theta)
{
  vector<Candidate> rotation_canditate_sets;

  int step_size = 2 * angular_search_window_ /  angular_search_step_;
  
  const double xcenter = (map_width_ / 2) * resolution_;
  const double ycenter = (map_height_ / 2) * resolution_;

  for(int k = 0; k < step_size; k++)
  {
    float rotated_angle = -angular_search_window_ + angular_search_step_ * k;
    vector<size_t> pts;
    for (size_t i = 0; i < scan.ranges.size(); ++i)
    {
      const double angle = angles::normalize_angle(scan.angle_min + i * scan.angle_increment + rotated_angle);
      double range = scan.ranges[i];
      double x1 = x + range * std::cos(angle) + xcenter; // Can be negative
      double y1 = y + range * std::sin(angle) + ycenter;
      const long int xmap = lround(x1 / resolution_);
      const long int ymap = lround(y1 / resolution_);
      // TO DO: Consider boundry !
      size_t point_index = (ymap * map_width_) + xmap;

      if(pointInMap(point_index))
      {
        pts.push_back(point_index);
      }
    }
    Candidate new_candidate;
    new_candidate.orientation = rotated_angle;
    new_candidate.discretize_scan = pts;
    rotation_canditate_sets.push_back(new_candidate);
  }

  return rotation_canditate_sets;

}



} // namespace namespace correlative_scan_math


