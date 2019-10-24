#include "nn_cluster.h"
#include "lsh.h"
#include "utilities.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>
#include <sstream>
//#include <assert.h>
#include <vector>

// maybe it should be named delta_ess
inline double nnCluster::distance(int size_a, int size_b, double dist) {
  return (size_a * size_b * dist) / ((double)size_a + size_b);
}

//constructor
nnCluster::nnCluster(std::vector<std::vector<double>> &points_, int n, int d, double epsilon_, int bucket, int bins, int run_time):
				size(n), dimension(d), epsilon(epsilon_), bucket_size(bucket), nb_bins(bins), running_time(run_time) {

  cluster_weight = std::vector<int>(n * 4);
	int nb_ds = (int) ceil(log_base_(n, 1 + epsilon));
	number_of_data_structure = (std::max(nb_ds, 1))  + 5;
	points = std::vector<std::vector<double>>(points_);

	std::cout << "epsilon " << epsilon << ' ' << number_of_data_structure << std::endl;
	nn_data_structures.reserve(number_of_data_structure);

  sizes = std::vector<int> (number_of_data_structure);
  //assert(sizes.size() == number_of_data_structure);

	for (int i = 0; i < number_of_data_structure; ++i) {
		LSHDataStructure idx(bucket , bins, d);
	  nn_data_structures.push_back(idx);
    sizes[i] = 0;
	}

  for (size_t i = 0; i < points.size(); ++i) add_cluster(points[i], 1, i);
}

/**
* output : id, distance, and weight (the true one)
*/
std::tuple<int, double, int> nnCluster::query(int q_id, const std::vector<double> &query, const int query_size) {
  double min_distance = std::numeric_limits<double>::max();
  int res = -1;
  int res_index = 0;
  for (int i = 0; i < number_of_data_structure; ++i) {
    if (sizes[i] <= 0) continue;

	 	auto p = nn_data_structures[i].QueryPoint(q_id, query, running_time);
    //std::cout << i << ' ' << p.first << std::endl;
    if (p.first == -1) continue;
    int tmp_index = p.first;
		int tmp_size = cluster_weight[p.first];
  //  assert(tmp_size > 0);
    double tmp_dist = distance(query_size, tmp_size, p.second);
    if (tmp_dist <= min_distance) {
      min_distance = tmp_dist;
      res = tmp_index;
      res_index = i;
	  }
  }
  return std::make_tuple(res, min_distance, cluster_weight[res]);
}

// I should add the cluster to the array of points
void nnCluster::add_cluster(const std::vector<double> &cluster, const int cluster_size, const int id) {
	    int ds = floor(log_base_(cluster_size, 1 + epsilon));
      id_ds.insert({id, ds});
      nn_data_structures[ds].InsertPoint(id, cluster);
			sizes[ds]++;
      cluster_weight[id] = cluster_size;
    //  assert(sizes[ds] <= size);
      if(id >= points.size()) points.push_back(cluster);
}

void nnCluster::put_back(const std::vector<double> &cluster, const int id) {
  int ds = id_ds[id];
  nn_data_structures[ds].InsertPoint(id, cluster);
  if(id >= points.size()) points.push_back(cluster);
}

void nnCluster::v_put_back(const int id) {
  int ds = id_ds[id];
  sizes[ds]++;
}

/**
* TO CHANGE
*/
void nnCluster::delete_cluster(int idx) {
  //assert(id_ds.find(idx) != id_ds.end());
  int ds = id_ds[idx];
  nn_data_structures[ds].RemovePoint(idx);
  std::cout << "sizes "<< ds << ' ' << sizes[ds] << std::endl;
  //assert(sizes[ds] >= 0);
}

void nnCluster::v_delete_cluster(int idx) {
  //assert(id_ds.find(idx) != id_ds.end());
  int ds = id_ds[idx];
  // std::cout << "ds index : " << ds << " ind " << idx  << " w: " << cluster_weight[{ds, idx}] << std::endl;
  //nn_data_structures[ds].RemovePoint(idx);
  sizes[ds]--;
  // std::cout << "sizes "<< ds << ' ' << sizes[ds] << std::endl;
  //assert(sizes[ds] >= 0);
}

double nnCluster::compute_min_dist(std::unordered_set<int> &unmerged_clusters, std::vector<bool> &existed) {
  double min_dis = std::numeric_limits<double>::max();
	for (int i = 0; i < size; ++i) {
  	auto res = get_point(i);

    auto t = query(i, res, 1);
    if(std::get<1>(t) > 0.0)
		min_dis = std::min(2 * std::get<1>(t), min_dis);

    //assert(min_dis > 0);

    cluster_weight[i] = 1;

    unmerged_clusters.insert(i);
		existed[i] = true;
  }
  return min_dis;
}

// GOOD
int nnCluster::get_number_of_data_structures() const {
  return nn_data_structures.size();
}

void nnCluster::update_size(int ds_index, int new_index, int size) {
	cluster_weight[new_index] = size;
}

double nnCluster::compute_max_dist(const std::vector<std::vector<double>> &points, const int n, const int d) {
	std::vector<double> max_pt(d);
	std::vector<double> min_pt(d);

	for (int i = 0; i < d; ++i) {
		max_pt[i] = std::numeric_limits<double>::min();
		min_pt[i] = std::numeric_limits<double>::max();
 	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < d; j++) {
			max_pt[j] = std::max(points[i][j], max_pt[j]);
			min_pt[j] = std::min(points[i][j], min_pt[j]);
		}
	}

	double dist = 0.0;
	double x, y;
	for (int i = 0; i < d; ++i) {
		x = max_pt[i];
		y = min_pt[i];
		dist += (x - y) * (x - y);
	}
	return dist;
}
