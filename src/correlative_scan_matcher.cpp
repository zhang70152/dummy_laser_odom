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
  last_x_ = 0;
  last_y_ = 0;
  last_theta_ = 0;
  fail_counter_ = 0;
}


bool correlativeScanMatcher::multiResolutionSearch(const sensor_msgs::LaserScan& scan, double& x, double& y, double& theta)
{


  vector<RotatedScan> rotated_scan_sets = generateRotationScan(scan, 0, 0, 0);
  vector<Candidate> candidates = generateLowestResolutionCell(rotated_scan_sets);
  
  for(int i = 0; i < candidates.size(); i++)
  {
    scoreCandidate(candidates[i], rotated_scan_sets);
  }
 
  Candidate best_candidate_in_low_resolution;
  int best_index=0;
  double max_score = LEAST_SCORE_NUMBER;
  for(int i = 0; i < candidates.size(); i++)
  {
    if(candidates[i].score > max_score)
    {
      max_score = candidates[i].score;
      best_candidate_in_low_resolution.scan_index = candidates[i].scan_index;
      best_candidate_in_low_resolution.x_offset = candidates[i].x_offset;
      best_candidate_in_low_resolution.y_offset = candidates[i].y_offset;
      best_candidate_in_low_resolution.orientation = candidates[i].orientation;
      best_candidate_in_low_resolution.score = candidates[i].score;
    }
  }
  //std::cout<<"search cells num:"<<candidates.size()/rotated_scan_sets.size()<<std::endl;
  
  int start_x = best_candidate_in_low_resolution.x_offset * 2;
  int start_y = best_candidate_in_low_resolution.y_offset * 2;
   std::cout<<"low resolution offset_x:"<<best_candidate_in_low_resolution.x_offset
            <<" offset_y:"<<best_candidate_in_low_resolution.y_offset
            <<" orientation:"<<best_candidate_in_low_resolution.orientation<<std::endl;
  Candidate best_candidate = recursiveSearch(max_depth_-1, best_candidate_in_low_resolution, start_x, start_y, rotated_scan_sets);

  x = best_candidate.x_offset * resolution_;
  y = best_candidate.y_offset * resolution_;
  theta = best_candidate.orientation;

  if(fabs(x - last_x_)>0.31 || fabs(y - last_y_)>0.31 || fabs(theta - last_theta_)>0.2 )
  {
    std::cout<<" large displacement! fail to match! "
             <<" delta x: "<<(double)fabs(x - last_x_)
             <<" delta y: "<<fabs(y - last_y_)
             <<" delta theta: "<<(double)fabs(theta - last_theta_)
             <<" max score:"<<max_score<<std::endl;
    //set zero increment for "failed search".
    x = 0;
    y = 0;
    theta = 0;
    fail_counter_++;
  }
  else{
    last_x_ = x;
    last_y_ = y;
    last_theta_ = theta;
    fail_counter_ = 0;
  }
  
  std::cout<<"final x: "<<x<<"final  y:"<<y<<" final theta:"<<theta<<" score:"<<best_candidate.score<<std::endl;
  if(fail_counter_>5)
  {
    return false;
  }
  else{
    return true;
  }
}



Candidate correlativeScanMatcher::recursiveSearch(
  int current_depth, Candidate best_candidate, int start_x, int start_y, const vector<RotatedScan>& rotated_scan_sets)
{
  if(current_depth<0)
  {
    return best_candidate;
  }
  else{
    vector<Candidate> higher_resolution_candidates = generateLayeredCandidates(current_depth,  start_x,  start_y, rotated_scan_sets);

    for(int i = 0; i < higher_resolution_candidates.size(); i++)
    {
      scoreCandidate(higher_resolution_candidates[i], rotated_scan_sets);
    }
    Candidate best_candidate;
    double max_score = LEAST_SCORE_NUMBER;
    for(int i = 0; i < higher_resolution_candidates.size(); i++)
    {
      if(higher_resolution_candidates[i].score > max_score)
      {
        max_score = higher_resolution_candidates[i].score;
        best_candidate = higher_resolution_candidates[i];
      }
    }
    std::cout<<"candidate num"<<higher_resolution_candidates.size()
    <<" depth:"<<current_depth<<" best score:"<<max_score
    <<" x_offset:"<<best_candidate.x_offset<<" y_offset:"<<best_candidate.y_offset
    <<" orientation:"<<best_candidate.orientation
    <<" max_score:"<< best_candidate.score<<std::endl;
    start_x = best_candidate.x_offset * 2;
    start_y = best_candidate.y_offset * 2;
    best_candidate = recursiveSearch(current_depth-1, best_candidate, start_x, start_y, rotated_scan_sets);
    return best_candidate;
  }
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

void correlativeScanMatcher::updateMapLookUpTable(const std::vector<double>& lookup_table)
{
      lookup_table_ = lookup_table;
      updateCellsLookupTable();
}


void correlativeScanMatcher::updateCellsLookupTable()
{
  layered_lookup_table_[0] = &lookup_table_;

  for(int k = 1; k <= max_depth_; k++)
  {
    vector<double>* intermediate_look_up_table = new vector<double>();
    int cell_length = 1 << k;
    int cell_number = map_width_ / cell_length;
    for(int i = 0; i < cell_number; i++)
    {
      int start_y = i * cell_length;
      if(i>0)
      {
        start_y--;
      }
      for(int j = 0; j < cell_number; j++)
      {
        int start_x = j * cell_length;
        if(j>0)
        {
          start_x--;
        }
        double max_prop = findMaxLogInCell(k, start_x, start_y, cell_length);

        //std::cout<<"START X:"<<start_x<<" START Y"<<start_y<<" cells length:"<<cell_length<<"   DEPTH"<<k<<" max prob"<<max_prop<<std::endl;
        intermediate_look_up_table->push_back(max_prop);
      }
    }
    layered_lookup_table_[k] = intermediate_look_up_table;
    //std::cout<<"depth:"<<i<<" look up table size"<<layered_lookup_table_[i]->size()<<"cells length:"<<cell_length<<std::endl;
  }


}

vector<Candidate> correlativeScanMatcher::generateLowestResolutionCell(const vector<RotatedScan>& rotated_scan_sets)
{
  vector<Candidate> low_resolution_candidates;
  int step_size = linear_search_window_ /  resolution_;
  int cell_length = 1 << max_depth_;
  int cell_number = step_size / cell_length;

  for(int i = -cell_number; i < cell_number; i++)
  {
    for(int j = -cell_number; j < cell_number; j++)
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
  int counter = 0;
  for(int i = 0; i < cell_length; i++)
  {
    for(int j = 0; j < cell_length; j++)
    {
      size_t point_index =  getPointIndex(start_x + j, start_y + i);
     
      if(pointInMap(0, point_index))
      {
        if(lookup_table_[point_index] > max_prop)
        {
          max_prop = lookup_table_[point_index];
        }
     
      }
      else{
         std::cout<<"invalid:"<<point_index<<std::endl;
      }
      counter++;
    }
  }

  return max_prop;
}

void correlativeScanMatcher::scoreCandidate(Candidate& candidate, const vector<RotatedScan>& rotated_scan_sets)
{
  double score = LEAST_SCORE_NUMBER;
  ScanInGrid scan_in_grid = rotated_scan_sets[candidate.scan_index].discretize_scan;
  int counter = 0;
  for(int i = 0; i < scan_in_grid.size(); i++)
  {
    //TODO: Condider boundry!
    int x,y;
    indexToXY(scan_in_grid[i], x, y);
    int cell_length = 1 << candidate.depth;
    int cell_x = x / cell_length;
    int cell_y = y / cell_length;

    size_t cell_index = getCellIndex(cell_x + candidate.x_offset, cell_y + candidate.y_offset, candidate.depth);

    int offset = getPointIndex(candidate.x_offset, candidate.y_offset);
    size_t point_index = scan_in_grid[i] + (size_t)offset;

    // if(candidate.depth == 0)
    // {
    //   std::cout<<"x:"<<x<<" y:"<<y<<" cell_index:"<<cell_index<<" point index"<<point_index<<" i:"<<scan_in_grid[i]<<std::endl;
    //   std::cout<<"offset x:"<<candidate.x_offset<<" offset y:"<<candidate.y_offset<<std::endl;
    // }
    std::vector<double>* lookup_table = getLayeredLookupTable(candidate.depth);
    if(lookup_table)
    {
      if(!pointInMap(candidate.depth, cell_index))
      {
        std::cout<<"Point out of boundry."<<"depth:"<<candidate.depth<<" index: "<<scan_in_grid[i]<<" cell:"<<cell_index<<" lk table size:"<<lookup_table->size()<<std::endl;
        counter++;
        continue;
      }
      double log_odds = (*lookup_table)[cell_index];
      double prob = (1 - 1 / (1 + std::exp(log_odds))) ;
      score = score + prob;
      if(candidate.depth == max_depth_)
      {
        // std::cout<<"x:"<<x<<" y:"<<y<<" cell_index:"<<cell_index<<" point index"<<point_index<<" i:"<<scan_in_grid[i]<<std::endl;
        // std::cout<<"offset x:"<<candidate.x_offset<<" offset y:"<<candidate.y_offset<<std::endl;
        // std::cout<<"cell_index:"<<cell_index<<" score"<<score<<" prob:"<<prob<<std::endl;
      }
      //std::cout<<"cell_index:"<<cell_index<<" score"<<score<<" once:"<<prob<<std::endl;
    }
    else{
      std::cout<<"Error: Cant find look up table!"<<std::endl;
    }

  }
  if(counter>50)
  {
    std::cout<<"invalid point number:"<<counter<<std::endl;
  }

  candidate.score = score;
  //std::cout<<"offset x:"<<candidate.x_offset<<" offset y:"<<candidate.y_offset<<" score"<<score<<std::endl;
  
}


} // namespace namespace correlative_scan_math


