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
  map_width_ = map_width;
  map_height_ = map_height;
  resolution_ = resolution;
  max_depth_ = 0;
}

void correlativeScanMatcher::bruteForceSearch(const sensor_msgs::LaserScan& scan, double& x, double& y, double& theta)
{
  vector<Candidate> result_candidates;


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


void correlativeScanMatcher::multiResolutionSearch(const sensor_msgs::LaserScan& scan, double& x, double& y, double& theta)
{
  vector<RotatedScan> rotated_scan_sets = generateRotationScan(scan, 0, 0, 0);
  vector<Candidate> candidates = generateLowestResolutionCell(rotated_scan_sets);
  
  for(int i = 0; i < candidates.size(); i++)
  {
    scoreCandidate(candidates[i], rotated_scan_sets);
  }
 
  Candidate best_candidate_in_low_resolution;
  double max_score = LEAST_SCORE_NUMBER;
  for(int i = 0; i < candidates.size(); i++)
  {
    if(candidates[i].score > max_score)
    {
      max_score = candidates[i].score;
      best_candidate_in_low_resolution.x_offset = candidates[i].x_offset;
      best_candidate_in_low_resolution.y_offset = candidates[i].y_offset;
      best_candidate_in_low_resolution.orientation = candidates[i].orientation;
      best_candidate_in_low_resolution.score = candidates[i].score;
    }
  }
  int start_x = best_candidate_in_low_resolution.x_offset * 2;
  int start_y = best_candidate_in_low_resolution.y_offset * 2;
  Candidate best_candidate = recursiveSearch(best_candidate_in_low_resolution, start_x, start_y, rotated_scan_sets);
   
  x = best_candidate.x_offset * resolution_;
  y = best_candidate.y_offset * resolution_;
  theta = best_candidate.orientation;
}




void correlativeScanMatcher::scoreCandidate(Candidate& candidate, const vector<RotatedScan>& rotated_scan_sets)
{
  double score = LEAST_SCORE_NUMBER;
  ScanInGrid scan_in_grid = rotated_scan_sets[candidate.scan_index].discretize_scan;
  for(int i = 0; i < scan_in_grid.size(); i++)
  {
    //TODO: Condider boundry!
    int offset = getPointIndex(candidate.x_offset, candidate.y_offset);
    size_t point_index = scan_in_grid[i] + (size_t)offset;

  
    std::vector<double>* lookup_table = getLayeredLookupTable(candidate.depth);
    if(lookup_table)
    {
      if(!pointInMap(candidate.depth, point_index))
      {
        std::cout<<"Point out of boundry."<<std::endl;
        continue;
      }
      score += (*lookup_table)[point_index];
    }
    
  }
  candidate.score = score;
}


Candidate correlativeScanMatcher::recursiveSearch(
  int current_depth, Candidate best_candidate, int start_x, int start_y, const vector<RotatedScan>& rotated_scan_sets)
{
  if(current_depth==0)
  {
    return false;
  }
  else{
    vector<Candidate> higher_resolution_candidates = generateLayeredCandidates(current_depth,  start_x,  start_y, rotated_scan_sets);
    scoreCandidate(higher_resolution_candidates, rotated_scan_sets);
  
    Candidate best_candidate;
    double max_score = LEAST_SCORE_NUMBER;
    for(int i = 0; i < candidates.size(); i++)
    {
      if(candidates[i].score > max_score)
      {
        max_score = candidates[i].score;
        best_candidate.scan_index = candidates[i].scan_index;
        best_candidate.depth = candidates[i].depth;
        best_candidate.x_offset = candidates[i].x_offset;
        best_candidate.y_offset = candidates[i].y_offset;
        best_candidate.orientation = candidates[i].orientation;
        best_candidate.score = candidates[i].score;
      }
    }

    start_x = best_candidate.x_offset * 2;
    start_y = best_candidate.y_offset * 2;
    best_candidate = recursiveSearch(current_depth-1, best_candidate, start_x, start_y, rotated_scan_sets);
    return best_candidate;
  }
}




int correlativeScanMatcher::getPointIndex(int x, int y)
{
  return  x + y * map_width_;
}

vector<RotatedScan> correlativeScanMatcher::generateRotationScan(const sensor_msgs::LaserScan& scan, float x, float y, float theta)
{
  vector<RotatedScan> rotation_scan_sets;

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

      if(pointInMap(0, point_index))
      {
        pts.push_back(point_index);
      }
    }
    RotatedScan rs;
    rs.index = k;
    rs.rotated_angle = rotated_angle;
    rs.discretize_scan = pts;
    rotation_scan_sets.push_back(rs);
  } 
  return rotation_scan_sets;
}

vector<Candidate> correlativeScanMatcher::generateAllCanditate(const vector<RotatedScan>& rotated_scan_sets)
{
  vector<Candidate> result_candidates;

  int step_size = linear_search_window_ /  linear_search_step_;

  for(int i = -step_size; i < step_size; i++)
  {
    for (int j = -step_size; j < step_size; ++j)
    {
      for(int k = 0; k < rotated_scan_sets.size(); ++k)
      {
        Candidate new_candidate;
        new_candidate.scan_index = rotated_scan_sets[i].index;
        new_candidate.x_offset = j;
        new_candidate.y_offset = i;
        new_candidate.orientation =  rotated_scan_sets[i].rotated_angle;
        new_candidate.score = LEAST_SCORE_NUMBER;
        result_candidates.push_back(new_candidate);
      }
    }
  }
  return result_candidates;
}


void correlativeScanMatcher::updateMapLookUpTable(const std::vector<double>& lookup_table)
{
      lookup_table_ = lookup_table;
      updateLayeredLookupTable();
}


void correlativeScanMatcher::updateLayeredLookupTable()
{
  layered_lookup_table_[0] = &lookup_table_;

  for(int i = 1; i <= max_depth_; i++)
  {
    vector<double>* intermediate_look_up_table = new vector<double>();
    int cell_length = 1 << max_depth_;
    int cell_number = map_width_ / cell_length;
    for(int i = 0; i < cell_number; i++)
    {
      int start_y = i * cell_length;
      for(int j = 0; j < cell_number; j++)
      {
        int start_x = j * cell_length;
        double max_prop = findMaxLogInCell(i, start_x, start_y, cell_length);
        intermediate_look_up_table->push_back(max_prop);
      }
    }
    layered_lookup_table_[i] = intermediate_look_up_table;
  }


}

vector<Candidate> correlativeScanMatcher::generateLowestResolutionCell(const vector<RotatedScan>& rotated_scan_sets)
{
  vector<Candidate> low_resolution_candidates;
  int step_size = linear_search_window_ /  linear_search_step_;
  int cell_length = 1 << max_depth_;
  int cell_number = map_width_ / cell_length;

  for(int i = -step_size; i < step_size; i = i+cell_length)
  {
    for(int j = -step_size; j < step_size; j = j+cell_length)
    {
      for(int k = 0; k < rotated_scan_sets.size(); ++k)
      {
        Candidate new_candidate;
        new_candidate.scan_index = rotated_scan_sets[k].index;
        new_candidate.depth = max_depth_;
        new_candidate.x_offset = j;
        new_candidate.y_offset = i;
        new_candidate.orientation = rotated_scan_sets[k].rotated_angle;
        new_candidate.score = LEAST_SCORE_NUMBER;
        low_resolution_candidates.push_back(new_candidate);
      }
    }
  }
  return low_resolution_candidates; 
}


vector<Candidate> correlativeScanMatcher::generateLayeredCandidates(
  int depth, int start_x, int start_y, const vector<RotatedScan>& rotated_scan_sets)
{
  vector<Candidate> layered_candidates;

  for(int i = 0; i < 2; i++)
  {
    for(int j = 0; j < 2; j++)
    {
      for(int k = 0; k < rotated_scan_sets.size(); ++k)
      {
        Candidate new_candidate;
        new_candidate.scan_index = rotated_scan_sets[k].index;
        new_candidate.depth = depth;
        new_candidate.x_offset = start_x + j;
        new_candidate.y_offset = start_y + i;
        new_candidate.orientation = rotated_scan_sets[k].rotated_angle;
        new_candidate.score = LEAST_SCORE_NUMBER;
        layered_candidates.push_back(new_candidate);
      }
    }
  }
  return layered_candidates; 
}


double correlativeScanMatcher::findMaxLogInCell(int depth, int start_x, int start_y, int cell_length)
{
  double max_prop = LEAST_SCORE_NUMBER;
  for(int i = 0; i < cell_length; i++)
  {
    for(int j = 0; j < cell_length; j++)
    {
      size_t point_index =  getPointIndex(start_x + j, start_y + i);
      if(pointInMap(depth, point_index))
      {
        if(lookup_table_[point_index] > max_prop)
        {
          max_prop = lookup_table_[point_index];
        }
      }
    }
  }
  return max_prop;
}

void correlativeScanMatcher::indexToXY(size_t index, int& x, int& y)
{
  y = index / map_width_;
  x = index % map_width_;
}


} // namespace namespace correlative_scan_math


