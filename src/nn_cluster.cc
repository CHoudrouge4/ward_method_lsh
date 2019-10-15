#include "nn_cluster.h"
#include "lsh.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>
#include <sstream>

/**
* TODO: add point by pushing back to the std::vector
* Question: is the distance of LSH squared ?
*/
inline double log_base(double num, double base) { return std::log(num) / std::log(base); }

// maybe it should be named delta_ess
inline double nnCluster::distance(int size_a, int size_b, double dist) {
  return (size_a * size_b * dist) / (size_a + size_b);
}

//constructor
nnCluster::nnCluster(std::vector<std::vector<double>> &points_, int n, int d, double epsilon_):
				size(n), dimension(d), epsilon(epsilon_) {

	int nb_ds = (int) ceil(log_base(n, 1 + epsilon));
	number_of_data_structure = (std::max(nb_ds, 1)) * 2 + 5;
	points = std::vector<std::vector<double>>(points_);

	std::cout << "epsilon " << epsilon << ' ' << number_of_data_structure << std::endl;

	LSHDataStructure index(1000 , 1, d); // change the number of bucket later // parametrise it
	// insert the points
	for (size_t i = 0; i < points_.size(); ++i) {
		index.InsertPoint(i, points_[i]);
		cluster_weight[{0, i}] = 1;
	}

	build = std::vector<bool>(number_of_data_structure, false);
	sizes = std::vector<int> (number_of_data_structure, 0);
	nn_data_structures.reserve(number_of_data_structure);
	sizes[0] = n;
	nn_data_structures.push_back(index);
	for (int i = 1; i < number_of_data_structure; ++i) {
		LSHDataStructure idx(1000 , 1, d);
	  nn_data_structures.push_back(idx);
	}
}

/**
* output : id, distance, and weight (the true one)
*/
std::tuple<int, double, int> nnCluster::query(const std::vector<double> &query, const int query_size, bool itself) {
  double min_distance = std::numeric_limits<double>::max();
  int res = -1;
  int res_index = 0;
  for (int i = 0; i < number_of_data_structure; ++i) {
    if (sizes[i] <= 0) continue;

	 	auto p = nn_data_structures[i].QueryPoint(query, 10);

    int tmp_index = p.first;
		int tmp_size = cluster_weight[{i, p.first}];
	  double tmp_dist;
    if (itself)
      tmp_dist = p.second; //recheck
    else
      tmp_dist = distance(query_size, tmp_size, p.second);

    if (tmp_dist <= min_distance) {
      min_distance = tmp_dist;
      res = tmp_index;
      res_index = i;
	  }
  }
  return std::make_tuple(res, min_distance, cluster_weight[{res_index, res}]);
}

// I should add the cluster to the array of points
int nnCluster::add_cluster(const std::vector<double> &query, int cluster_size, int id) {
	    int idx = floor(log_base(cluster_size, 1 + epsilon));
      nn_data_structures[idx].InsertPoint(id, query);
			sizes[idx] = sizes[idx] + 1;
      cluster_weight[{idx, id}] = cluster_size ;
      std::cout << "cluster " << id << " and weight " << cluster_size << " is added to DS " << idx << '\n';
      if(id >= points.size()) points.push_back(query);
			return idx;
}

/**
* TO CHANGE
*/
void nnCluster::delete_cluster(int idx, int size) {
  int i = (int) floor(log_base(size, 1 + epsilon));// I think it is wrong
  nn_data_structures[i].RemovePoint(idx);
	sizes[i] = sizes[i] - 1;
}

std::vector<double> nnCluster::get_point(int idx) {
  std::cout << "getting point " << idx << " from total number of points " << points.size() << std::endl;
	return points[idx];
}

// good + need review
double nnCluster::compute_min_dist(std::unordered_set<pair_int> &unmerged_clusters, std::unordered_map<pair_int, bool, pairhash> &existed) {
  double min_dis = std::numeric_limits<double>::max();
	for (int i = 0; i < size; ++i) {
  	auto res = get_point(i);
    delete_cluster(i, 1);
    auto t = query(res, 1);
    std::cout << "The nearest neighbor of cluster " << i << " is cluster " << std::get<0>(t) << std::endl;
		min_dis = std::min(std::get<1>(t), min_dis);

    iszero(min_dis < 0);

		add_cluster(res, 1, i);
		dict[{i, 1}] = {i, 1};


		dict[{i, 1}] = {i, 1};
		cluster_weight[{0, i}] = 1;

		unmerged_clusters.insert({i, 1});
		existed[{i, 1}] = true; // check it
  }
  std::cout << "the minimum distance is: " << min_dis + min_dis << std::endl;
  return min_dis + min_dis;
}

// GOOD
int nnCluster::get_number_of_data_structures() const {
  return nn_data_structures.size();
}

// to use while building the representation
pair_int nnCluster::get_index(int index, int weight) {
	return dict[std::make_pair(index, weight)];
}

// we can later use one weight
void nnCluster::update_dict(int new_idx, int new_weight, int old_idx, int old_weight) {
	dict[{new_idx, new_weight}] = dict[{old_idx, old_weight}];
}

void nnCluster::update_size(int ds_index, int new_index, int size) {
	cluster_weight[{ds_index, new_index}] = size;
}

double nnCluster::compute_max_dist(const std::vector<std::vector<double>> points, const int n, const int d) {
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
  std::cout << "The maximum distance is: " << 2 * dist << std::endl;
	return 2 * dist;
}
