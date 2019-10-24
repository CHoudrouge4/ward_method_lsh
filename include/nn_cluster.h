#pragma once
#include "lsh.h"
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <vector>

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

class nnCluster {
public:
  // Constructor
  nnCluster (std::vector<std::vector<double>> &points_, int n, int d, double epsilon_,
                                              int bucket, int bins, int running_time);
  /**
  * This function return a nearest neighbour cluster id, distance to query,
  * and the weight of the nearest neighbour.
  */
  std::tuple<int, double, int> query(int q_id, const std::vector<double> &query, int query_size);

  /**
  * This function adds a cluster to the data clusters of the clusters of size
  * cluster size
  */
  void add_cluster(const std::vector<double> &cluster, const int cluster_size, const int id);
  void put_back(const std::vector<double> &cluster, const int id);
  void v_put_back(const int id);
  /**
  * This function deletes a cluseter with the id = idx, and weight = size.
  */
  void delete_cluster(int idx);
  void v_delete_cluster(int idx);
  /**
  * This function returns the coordinates of the centroid with id = idx, and weight = size
  */
  std::vector<double> get_point(int idx) {
    	return points[idx];
  }

  /**
  * This function returns the number of data structures.
  */
  int get_number_of_data_structures() const;

  /**
  * This function computes and returns the minimum distance between all the cluster.
  */
  double compute_min_dist(std::unordered_set<int> &unmerged_clusters, std::vector<bool> &existed);

  /**
  * This function computes and returns an approximation of the maximum distance between all the clusters.
  */
  double compute_max_dist(const std::vector<std::vector<double>> &points, const int n, const int d);

  pair_int get_index(int index, int weight);

  void update_size(int ds_index, int new_index, int size);

  int get_cluster_size(const int id) {
    return cluster_weight[id];
  }

private:
  std::vector<std::vector<double>> points; // this stores the initial data points
  // size is the number of input points, number of data structures, number of visted leaves
  int size, dimension, number_of_data_structure;
  int bucket_size, nb_bins, running_time;

  // epsilon, gamma initially for appriximating nearest neighbour distance, maximum distance.
  double epsilon, max_distance;

  // this vectors stores the nearest neighbour data structures
  std::vector<LSHDataStructure> nn_data_structures;
  // the number of points inside each data structures
  std::vector<int> sizes;

  // computes the ward's distance between two clusters of size_q, dise_b
  inline double distance(int size_a, int size_b, double dist);

  // maps the cluster to its weight, the cluster is uniquely determined by
  // the index of the data structure and its index.
  std::vector<int> cluster_weight;
  // std::unordered_map<pair_int, int> idx_index;
  std::unordered_map<int, int> id_ds;
};
