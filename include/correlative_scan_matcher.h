#ifndef _CORRELATIVE_SCAN_MATCHER_H_
#define _CORRELATIVE_SCAN_MATCHER_H_

#include <cmath> // For std::exp.
#include <exception>
#include <string>
#include <vector>
#include <map>
#include <ros/ros.h>
#include <angles/angles.h>
#include <sensor_msgs/LaserScan.h>
#include <nav_msgs/OccupancyGrid.h>



#define LEAST_SCORE_NUMBER -9999.0
namespace correlative_scan_math
{

using std::vector;

typedef std::vector<size_t>  ScanInGrid;

struct RotatedScan {
  int index;
  double rotated_angle;
  ScanInGrid discretize_scan;
};

struct Candidate {
  int scan_index;
  int depth;
  int x_offset = 0;
  int y_offset = 0;
  double x = 0.;
  double y = 0.;
  double orientation = 0.;

  double score = LEAST_SCORE_NUMBER;

  bool operator<(const Candidate& other) const {return score < other.score;}
  bool operator>(const Candidate& other) const {return score > other.score;}
};


class correlativeScanMatcher
{
  public:
    correlativeScanMatcher(int map_width, int map_height, float resolution);

    ~correlativeScanMatcher()
    {
      
    }
    void updateMapLookUpTable(const std::vector<double>& lookup_table);

    void bruteForceSearch(const sensor_msgs::LaserScan& scan, double& x, double& y, double& theta);

    void multiResolutionSearch(const sensor_msgs::LaserScan& scan, double& x, double& y, double& theta);

    void setSearchParameters(float linear_search_window, float linear_search_step, float angular_search_window, float angular_search_step)
    {
      linear_search_window_ = linear_search_window;
      linear_search_step_ = linear_search_step;
      angular_search_window_ = angular_search_window; 
      angular_search_step_ = angular_search_step;
    }

  private:

    vector<RotatedScan> generateRotationScan(const sensor_msgs::LaserScan& scan, float x, float y, float theta);

    vector<Candidate> generateAllCanditate(const vector<RotatedScan>& rotated_scan_sets);

    void scoreCandidate(Candidate& candidate, const vector<RotatedScan>& rotated_scan_sets);

    Candidate recursiveSearch(int current_depth,  Candidate best_candidate, int start_x, int start_y, const vector<RotatedScan>& rotated_scan_sets);

    int getPointIndex(int x, int y);

    bool pointInMap(int depth, size_t index)
    {
      std::map<int, vector<double>*>::iterator it = layered_lookup_table_.find(depth);
      if(it != layered_lookup_table_.end())
      {
        vector<double>* lookup_table = it->second;
        return index < lookup_table->size();
      }
      else{
        return false;
      }
    }

    void updateLayeredLookupTable();

    vector<Candidate> generateLowestResolutionCell(const vector<RotatedScan>& rotated_scan_sets);

    vector<Candidate> generateLayeredCandidates(int depth, int start_x, int start_y, const vector<RotatedScan>& rotated_scan_sets);

    double findMaxLogInCell(int depth, int start_x, int start_y, int cell_length);

    void indexToXY(size_t index, int& x, int& y);

    std::vector<double>* getLayeredLookupTable(int depth)
    {
      std::map<int, vector<double>*>::iterator it = layered_lookup_table_.find(depth);
      if(it != layered_lookup_table_.end())
      {
        vector<double>* lookup_table = it->second;
        return lookup_table;
      }
      else{
        return 0;
      }
    }

    std::vector<double> lookup_table_;  //!< log odds ratios for the binary Bayes filter
                                    //!< log_odd = log(p(x) / (1 - p(x)))

    std::map<int, vector<double>*> layered_lookup_table_;

    int map_width_;
    int map_height_;
    float resolution_;


    float linear_search_window_;
    float linear_search_step_;
    float angular_search_window_; 
    float angular_search_step_;

    int max_depth_;
};


} 

#endif // # ifndef _CORRELATIVE_SCAN_MATCHER_H_
