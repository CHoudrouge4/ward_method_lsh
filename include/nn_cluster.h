#pragma once
#include "lsh.h"
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <memory>

typedef std::pair<int, int> pair_int;

/**
* A hash function for pair of integer.
*/
namespace std {
  template <>
  struct hash<pair_int> {
    const size_t num = 65537;
    inline size_t operator()(const std::pair<int, int> &x) const {
        return (x.first * num) ^ x.second;
    }
  };
}

/**
* Struct that implements a hash function for pair of integers
*
*/
struct pairhash {
private:
  const size_t num = 65537;
public:
  inline size_t operator()(const std::pair<int, int> &x) const {
    return (x.first * num) ^ x.second;
  }
};

class nnCluster {
public:
  // Constructor
  nnCluster (std::vector<std::vector<double>> * points_, int n, int d, double epsilon_, double gamma_, const size_t &tree_number, int visited_leaf_);

  /**
  * This function return a nearest neighbour cluster id, distance to query,
  * and the weight of the nearest neighbour.
  */
  std::tuple<int, double, int> query(const std::vector<double> &query, int query_size, bool itself=false);

  /**
  * This function adds a cluster to the data clusters of the clusters of size
  * cluster size
  *
  */
  int add_cluster(const std::vector<double> &cluster, const int cluster_size);

  /**
  * This function adds a cluster of a given size to the corresponding data structure
  */
  int add_cluster(const std::vector<double> &cluster, int cluster_size, int old_index, int new_index );


  std::tuple<int, double, int> add_new_cluster(const std::vector<double> &cluster, const int cluster_size);

  /**
  * This function deletes a cluseter with the id = idx, and weight = size.
  */
  void delete_cluster(int idx, int size);

  /**
  * This function returns the coordinates of the centroid with id = idx, and weight = size
  */
  double * get_point(int idx, int size);

  /**
  * This function returns the number of data structures.
  */
  int get_number_of_data_structures() const;

  /**
  * This function computes and returns the minimum distance between all the cluster.
  */
  double compute_min_dist(std::unordered_set<pair_int> &unmerged_clusters, std::unordered_map<pair_int, bool, pairhash> &existed);

  /**
  * This function computes and returns an approximation of the maximum distance between all the clusters.
  */
  double compute_max_dist(const double * points, const int n, const int d);

  pair_int get_index(int index, int weight);

  void update_dict(int new_idx, int new_weight, int old_idx, int old_weight);

  void update_size(int ds_index, int new_index, int size);

private:
  std::vector<std::vector<double>> points; // this stores the initial data points
  // size is the number of input points, number of data structures, number of visted leaves
  int size, dimension, number_of_data_structure, visited_leaf;
  // epsilon, gamma initially for appriximating nearest neighbour distance, maximum distance.
  double epsilon, gamma, max_distance;

  // this vectors stores the nearest neighbour data structures
  std::vector<LSHDataStructure> nn_data_structures;
  // entry i equals to true of the ith data structures is built
  std::vector<bool> build;
  // the number of points inside each data structures
  std::vector<int> sizes;
  // stores the indices returned from the NN query
  std::vector<std::vector<int>> indices;
  // stores the distance returned from the NN query
  std::vector<std::vector<double>> dists;

  // computes the ward's distance between two clusters of size_q, dise_b
  inline double distance(int size_a, int size_b, double dist);

  // maps the cluster to its weight, the cluster is uniquely determined by
  // the index of the data structure and its index.
  std::unordered_map<std::pair<int, int>, int> cluster_weight;
  // this maps helps in having a unique ID for each cluster.
  std::unordered_map<pair_int, pair_int> dict;
  // std::unordered_map<pair_int, int> idx_index;
};
