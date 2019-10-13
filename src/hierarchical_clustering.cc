#include "hierarchical_clustering.h"
#include "utilities.h"
#include <fstream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <tuple>
#include <utility>
#include <sstream>
//
/**
* TODO:
*   Remove the log_base_ funstion
*   check the precision of fpa
*   make unchecked a field
*   don't use tuples if possible
*   make the removal of elments from the unordered_map faster.
*/

#define id first
#define w second

typedef std::pair<int, int> pair_int;

inline std::string toString(const std::pair< size_t, size_t> & data) {
     std::ostringstream str;
     str << data.first << "," << data.second << ";";
     return str.str();
}

hierarchical_clustering::hierarchical_clustering(std::vector<std::vector<double>> &data, int n, int d, double epsilon_, double gamma_):
                                                               nnc(data, n, d, epsilon_, gamma_),
                                                               dimension(d),
                                                               size(n),
                                                               epsilon(epsilon_),
                                                               gamma(gamma_) {

  max_dist = nnc.compute_max_dist(data, n, d);
  unmerged_clusters.max_load_factor(std::numeric_limits<double>::infinity());
  min_dist = nnc.compute_min_dist(unmerged_clusters, existed);
  beta = ceil(log_base_((max_dist/min_dist) * n, 1 + epsilon)); // be carefull four the double / double
  output.reserve(n);
  to_erase.reserve(n);
}

/**
* this function takes the centroid of two clusters and computes the new centroid
* based on the following formula.
* new_centroid = ((size_a * mu_a) + (size_b * mu_b))/(size_a + size_b).
*
*/
inline std::vector<double> hierarchical_clustering::merge(std::vector<double> &mu_a, std::vector<double> &mu_b, int size_a, int size_b) {
  double den = size_a + size_b;
  double coeff_a = size_a / den; // make sure if this is efficient
  double coeff_b = size_b / den;
  std::vector<double> res(dimension);
  for (int i = 0; i < dimension; ++i)
    res[i] = coeff_a * mu_a[i] + coeff_b * mu_b[i];
  return res;
}

std::unordered_set<pair_int>  hierarchical_clustering::helper(std::unordered_set<pair_int> &mp, double merge_value) {
  std::unordered_set <pair_int> unchecked; // this one should be placed maybe in different place, maybe it should be in the fields
  std::vector<double> merged_cluster;
  merged_cluster.reserve(dimension);
  int merged_weight;
  for (auto&& p : mp) {
     if(existed[p]) {
     bool ok = false;
     bool flag = false;
     int u = p.id;
     int u_weight = p.w;
     if(u_weight == 0) continue;
     to_erase.push_back({u, u_weight});

     //double * res_ = (double *) malloc(dimension * sizeof(double));
      auto res = nnc.get_point(u);
      nnc.delete_cluster(u, u_weight); // we delete the cluster because it is the nn of itself

      auto t = nnc.query(res, u_weight);
     //std::cout << "t: " << std::get<0>(t) << ' ' << std::get<1>(t) << ' ' << std::get<2>(t) << std::endl;

     double dist = std::get<1>(t); // getting the distance
     int t_weight = std::get<2>(t);


     while (dist < merge_value) {
       ok = true;

       auto nn_pt = nnc.get_point(std::get<0>(t));// to check

       // merging phase
       merged_cluster = merge(res, nn_pt, u_weight, t_weight);
       merged_weight = u_weight + t_weight;

       {
           int u_tmp = u;
           int weight_tmp = u_weight;
           auto p = nnc.add_new_cluster(merged_cluster, merged_weight, last_index);

           u = std::get<0>(p);
           u_weight = merged_weight;
           nnc.delete_cluster(u, u_weight);

           output.push_back(std::make_tuple(nnc.get_index(u, u_weight), nnc.get_index(u_tmp, weight_tmp), nnc.get_index(std::get<0>(t), std::get<2>(t))));
       }
//
        existed[p] = false;
        to_erase.push_back({u, u_weight});

        existed[{std::get<0>(t), t_weight}] = false;
        to_erase.push_back({std::get<0>(t), t_weight});
        nnc.delete_cluster(std::get<0>(t), t_weight);

        if(merged_weight == size) {
          stop = true;
          break;
        }

        t        = nnc.query(merged_cluster, merged_weight);
        dist     = std::get<1>(t);
        t_weight = std::get<2>(t);

        if(t_weight <= 0 || std::get<0>(t) < 0) break;
        auto nnn_pt = nnc.get_point(std::get<0>(t));
//          if(nnn_pt == nullptr) break;
        if(dist < merge_value) {
          res  = merged_cluster;
          flag = true;
        }
    }

    if(u_weight == size) break;
    if(u_weight == 0) break;
    if(!ok) {
     //assert(res_ != nullptr);
    int idx = nnc.add_cluster(res, u_weight, last_index);
    t = nnc.query(res, u_weight, true);
    nnc.update_size(idx, std::get<0>(t), u_weight);
    t_weight = std::get<2>(t);

    nnc.update_dict(std::get<0>(t), u_weight, u, u_weight);

    existed[{std::get<0>(t), u_weight}] = true;
    magic.insert({std::get<0>(t), u_weight});
   } else {
     // fix the id
     int idx = nnc.add_cluster(merged_cluster, merged_weight, last_index);
     auto tt = nnc.query(merged_cluster, merged_weight, true);
     nnc.update_size(idx, std::get<0>(tt), merged_weight);

     tt = nnc.query(merged_cluster, merged_weight, true);

     existed[{std::get<0>(tt), merged_weight}] = true;
     unchecked.insert({std::get<0>(tt), merged_weight});
     nnc.update_dict(std::get<0>(tt), merged_weight, u, u_weight);
     }
    }
   }

   for(size_t i = 0; i < to_erase.size(); ++i) mp.erase(to_erase[i]);
   to_erase.clear();
   return unchecked;
}
//
void hierarchical_clustering::build_hierarchy() {
  double merge_value;
  for (int i = 0; i < beta; ++i) {
    merge_value = pow(1 + epsilon, i); // find an efficient one
    //std::cout << "merge value " << merge_value << std::endl;

    auto ss = helper(this->unmerged_clusters, merge_value); // these are the merges
    while (ss.size() > 1) {
      auto tmp = helper(ss, merge_value);
      ss.clear();
      for(auto&& p: tmp) {
        existed[p] = true;
        ss.insert(p);
      }
     if(stop) break;
   }

   if(stop) return;
   if(ss.size() == 1) unmerged_clusters.insert(*ss.begin());
   for(auto m: magic) unmerged_clusters.insert(m);

     magic.clear();
     if(unmerged_clusters.size() <= 1) break;
  }
}
//
std::vector<std::pair<pair_int, pair_int>> hierarchical_clustering::get_merges() const {
   return merges;
}

void hierarchical_clustering::print_merges() {
   for (auto&& t: output) {
     std::cout << toString(std::get<0>(t)) << toString(std::get<1>(t)) << toString(std::get<2>(t)) << std::endl;
   }
}

void hierarchical_clustering::print_file(const std::string filename) {
  std::ofstream out(filename);
  for (auto&& t: output)
    out << toString(std::get<0>(t)) << toString(std::get<1>(t)) << toString(std::get<2>(t)) << std::endl;
  out.close();
}
