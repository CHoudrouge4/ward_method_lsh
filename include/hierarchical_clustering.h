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
  float epsilon; // epsilion as it is stated in the main paper.
  float gamma;   // essentially created to ensure the approximation of the NN
  double max_dist; // the maximum distance between any two points in the dataset
  double min_dist; // the minimum distance between any two points in the dataset
  double beta;     // size * max_dist / min_dist

  nnCluster nnc;   // a nearest neighbour data structure for clusters

  bool stop = false;

  std::vector<std::pair<pair_int, pair_int>> merges;
  std::unordered_map<pair_int, bool, pairhash> existed;
  std::unordered_set<pair_int> magic;
  std::unordered_set<std::pair<int, int>> unmerged_clusters;
  std::vector<std::tuple<pair_int, pair_int, pair_int>> output;
  std::vector<pair_int> to_erase;

  inline float * merge(float * mu_a, float * mu_b, int size_a, int size_b);

  std::unordered_set<pair_int> helper(std::unordered_set<pair_int> &mp, float merge_value);

public:
  hierarchical_clustering(float * data, int n, int d, float epsilon_, float gamma_, int tree_number, int visited_leaf);
  std::vector<std::pair<pair_int, pair_int>> get_merges() const;
  void print_merges();
  void build_hierarchy();
  void print_file(const std::string filename);
};
