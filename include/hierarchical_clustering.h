#pragma once
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include "nn_cluster.h"
#include <memory>

typedef std::vector<double> point;
typedef unsigned int uint;
typedef std::pair<int, int> pair_int;

/**
* Hash function for the pair of integers.
*/
class hierarchical_clustering {

private:
  int dimension; // store the dimension of the data
  int size;      // store the size of the data
  int last_index; // to keep track of the last index
  double epsilon; // epsilion as it is stated in the main paper.
  double max_dist; // the maximum distance between any two points in the dataset
  double min_dist; // the minimum distance between any two points in the dataset
  double beta;     // size * max_dist / min_dist

  nnCluster nnc;   // a nearest neighbour data structure for clusters

  bool stop = false;
  int merged_weight;
  point merged_cluster;
  std::vector<std::pair<pair_int, pair_int>> merges;
  std::vector<bool> existed;
  std::unordered_set<int> lambda;
  std::unordered_set<int> unmerged_clusters;
  std::vector<std::tuple<pair_int, pair_int, pair_int>> output;
  void merge(point &mu_a, point &mu_b, int size_a, int size_b);
  std::unordered_set <int> unchecked;
  std::unordered_set<int> helper(std::unordered_set<int> &mp, double merge_value);

public:
  hierarchical_clustering(std::vector<point> &data, int n, int d, double epsilon_, int bucket, int bins, int run_time);
  std::vector<std::pair<pair_int, pair_int>> get_merges() const;
  void print_merges();
  void build_hierarchy();
  void print_file(const std::string filename);

};
