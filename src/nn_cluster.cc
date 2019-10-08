#include "nn_cluster.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>
#include <sstream>

/**
* TODO:
* 	create a function just to query the point itself.
*  	find a way to get the right number of data structures
*   Add memo for the size to index
*/

template <class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
	if(v.size() == 0) {
		out << "()";
		return out;
	} else  {
		out << '(';
		out << ' ' << v[0] ;
		for(size_t i = 1; i < v.size(); ++i)
			out << " , " << v[i];
		out << ')';
		return out;
	}
}

inline double log_base(double num, double base) { return std::log(num) / std::log(base); }

inline double nnCluster::distance(int size_a, int size_b, double dist) {
  return (size_a * size_b * dist) /(size_a + size_b);
}

nnCluster::nnCluster(std::vector<vector<double>> &points_, int n, int d, double epsilon_, double gamma_, const size_t &tree_number, int visited_leaf_):
      points(points_, n, d), size(n), dimension(d), epsilon(epsilon_) , gamma(gamma_), visited_leaf(visited_leaf_) {

	//std::cout << "n " << n << std::endl;
	//assert(visited_leaf == visited_leaf_);
	int nb_ds = (int) ceil(log_base(n, 1 + epsilon));
	number_of_data_structure = (std::max(nb_ds, 1) )* 2 + 5;
	std::cout << "epsilon " << epsilon << ' ' << number_of_data_structure << std::endl;
	LSHDataStructure index(1000 , 1, d); // change the number of bucket later
	// insert the points
	for (int i = 0; i < points_.size(); ++i) {
		index.InsertPoint(i, points_[i]);
	}

	build = std::vector<bool>(number_of_data_structure, false);
	sizes = std::vector<int> (number_of_data_structure, 0);
	nn_data_structures.reserve(number_of_data_structure);
	sizes[0] = n;
	nn_data_structures.push_back(index);
	nn_data_structures[0].buildIndex();
	build[0] = true;
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
    if (!build[i] || sizes[i] <= 0) continue;
		nn_data_structures[i].knnSearch(query, indices, dists, 1, flann::SearchParams(visited_leaf));
		if(indices.size() == 0) continue;
		if(indices[0].size() == 0) continue;
    int tmp_index = indices[0][0];
		int tmp_size = cluster_weight[{i, tmp_index}];
	  double tmp_dist;
    if (itself)
      tmp_dist = dists[0][0];
    else
      tmp_dist = distance(query_size, tmp_size, dists[0][0]);

    if (tmp_dist <= min_distance) {
      min_distance = tmp_dist;
      res = tmp_index;
      res_index = i;
	  }

    indices[0].clear();
    dists[0].clear();
  }
	//
  return std::make_tuple(res, min_distance, cluster_weight[{res_index, res}]);
}

int nnCluster::add_cluster(const std::vector<double> &query, int cluster_size) {
      int idx = floor(log_base(cluster_size, 1 + epsilon));
			//assert(idx < number_of_data_structure);
			//assert(idx >= 0);
		  if (!build[idx]) {
        nn_data_structures[idx].buildIndex(cluster);
        build[idx] = true;
      } else {
        nn_data_structures[idx].addPoints(cluster);
      }
			sizes[idx] = sizes[idx] + 1;
			return idx;
}

int nnCluster::add_cluster(const std::vector<double> &cluster, int cluster_size, int old_index, int new_index) {
	int idx = add_cluster(cluster, cluster_size);
	dict[{new_index, cluster_size}] = dict[{old_index, cluster_size}];
//	idx_index[{idx, new_index}] = idx_index[{idx, old_index}];
	return idx;
}

std::tuple<int, double, int> nnCluster::add_new_cluster(const std::vector<double> &cluster, const int cluster_size) {
	int idx = add_cluster(cluster, cluster_size);
	auto t = query(cluster, cluster_size, true);
	dict[{std::get<0>(t), cluster_size}] = {std::get<0>(t), cluster_size};
	cluster_weight[{idx, std::get<0>(t)}] = cluster_size;
	return {std::get<0>(t), std::get<1>(t), cluster_size};
}

void nnCluster::delete_cluster(int idx, int size) {
  int i = (int) floor(log_base(size, 1 + epsilon));// I think it is wrong
	//assert(i >= 0);
	//assert((size_t)i < nn_data_structures.size());
  nn_data_structures[i].removePoint(idx);
	sizes[i] = sizes[i] - 1;
}

double * nnCluster::get_point(int idx, int size) {
	int i = (int) floor(log_base(size, 1 + epsilon));
//	std::cout << idx << ' ' << size << std::endl;
	//assert(i < number_of_data_structure);
	//assert(i >= 0);
  if(build[i]) {
    return nn_data_structures[i].getPoint(idx);
  } else {
	  return nullptr;
	}
}

double nnCluster::compute_min_dist(std::unordered_set<pair_int> &unmerged_clusters, std::unordered_map<pair_int, bool, pairhash> &existed) {
  double min_dis = std::numeric_limits<double>::max();
  for (int i = 0; i < size; ++i) {
      double * res_ = get_point(i, 1);
// edit
			flann::Matrix<double> res(res_, 1, dimension);
      delete_cluster(i, 1);
      auto t = query(res, 1);
			min_dis = std::min(std::get<1>(t), min_dis);
			add_cluster(res, 1);
			cluster_weight[{0, i}] = 1;
			dict[{i, 1}] = {i, 1};

			t = query(res, 1, true);
      dict[{std::get<0>(t), 1}] = {i, 1};
			cluster_weight[{0, std::get<0>(t)}] = 1;

			unmerged_clusters.insert({std::get<0>(t), 1});
			existed[{std::get<0>(t), cluster_weight[{0, std::get<0>(t)}]}] = true;

		//	std::cout << "( " << std::get<0>(t) << ' ' <<  cluster_weight[{0, std::get<0>(t)}] << " ) " << std::endl;
  }
  return min_dis;
}

int nnCluster::get_number_of_data_structures() const {
  return nn_data_structures.size();
}

// to use while building the representation
pair_int nnCluster::get_index(int index, int weight) {
	return dict[std::make_pair(index, weight)];
}

// we can later use one weight
void nnCluster::update_dict(int new_idx, int new_weight, int old_idx, int old_weight) {
	//std::cout << "new_idx " << new_idx << " new_weight " << new_weight << " old_idx " << old_idx << " old_weight " << old_weight << std::endl;
	dict[{new_idx, new_weight}] = dict[{old_idx, old_weight}];
}

void nnCluster::update_size(int ds_index, int new_index, int size) {
	cluster_weight[{ds_index, new_index}] = size;
}

nnCluster::~nnCluster() {
  delete [] points.ptr();
	// for(size_t i = 0; i < nn_data_structures.size(); ++i)
	// 	delete nn_data_structures[i];
}

double nnCluster::compute_max_dist(const double * points, const int n, const int d) {
	double * max_pt = (double *) malloc(d * sizeof(double));
	double * min_pt = (double *) malloc(d * sizeof(double));

	for (int i = 0; i < d; ++i) {
		* (max_pt + i) = std::numeric_limits<double>::min();
		* (min_pt + i) = std::numeric_limits<double>::max();
 	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < d; j++) {
			* (max_pt + j) = std::max(*(points + i * d + j), *(max_pt + j));
			* (min_pt + j) = std::min(*(points + i * d + j), *(min_pt + j));
		}
	}

	double dist = 0.0;
	double x, y;
	for (int i = 0; i < d; ++i) {
		x = * (max_pt + i);
		y = * (min_pt + i);
		dist += (x - y) * (x - y) + (x - y) * (x - y);
	}

	free(max_pt);
	free(min_pt);
	return dist;
}
