#ifndef _CORRELATIVE_SCAN_MATCHER_H_
#define _CORRELATIVE_SCAN_MATCHER_H_

#include <cmath> // For std::exp.
#include <exception>
#include <string>
#include <vector>

#include <ros/ros.h>
#include <angles/angles.h>
#include <sensor_msgs/LaserScan.h>
#include <nav_msgs/OccupancyGrid.h>



#define LEAST_SCORE_NUMBER -9999.0
namespace correlative_scan_math
{

using std::vector;

typedef std::vector<size_t>  ScanInGrid;
//typedef std::map<float, ScanInGrid > CandidateSets;
struct Candidate {
  ScanInGrid discretize_scan;
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
    void updateMapLookUpTable(const std::vector<double>& lookup_table)
    {
      lookup_table_ = lookup_table;
    }

    void bruteForceSearch(const sensor_msgs::LaserScan& scan, double& x, double& y, double& theta);

    void setSearchParameters(float linear_search_window, float linear_search_step, float angular_search_window, float angular_search_step)
    {
      linear_search_window_ = linear_search_window;
      linear_search_step_ = linear_search_step;
      angular_search_window_ = angular_search_window; 
      angular_search_step_ = angular_search_step;
    }

  private:

    vector<Candidate> generateRotationScan(const sensor_msgs::LaserScan& scan, float x, float y, float theta);

    void scoreCandidate(Candidate& candidate);

    int getPointIndex(int x, int y);

    bool pointInMap(size_t index)
    {
      return index < map_width_ * map_height_;
    }

    std::vector<double> lookup_table_;  //!< log odds ratios for the binary Bayes filter
                                    //!< log_odd = log(p(x) / (1 - p(x)))

    int map_width_;
    int map_height_;
    float resolution_;


    float linear_search_window_;
    float linear_search_step_;
    float angular_search_window_; 
    float angular_search_step_;
};


} 

#endif // # ifndef _CORRELATIVE_SCAN_MATCHER_H_
