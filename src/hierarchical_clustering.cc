#include "hierarchical_clustering.h"
#include "utilities.h"
#include <fstream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <tuple>
#include <utility>
#include <sstream>
//#include <assert.h>

#define id first
#define w second

typedef std::pair<int, int> pair_int;

inline std::string toString(const std::pair< size_t, size_t> & data) {
     std::ostringstream str;
     str << data.first << "," << data.second << ";";
     return str.str();
}

hierarchical_clustering::hierarchical_clustering(std::vector<point> &data, int n, int d, double epsilon_, int buckets, int bins, int run_time):
                                                               nnc(data, n, d, epsilon_, buckets, bins, run_time),
                                                               dimension(d),
                                                               size(n),
                                                               epsilon(epsilon_) {
  merged_cluster = point(dimension);
  last_index = data.size();
  max_dist = 10 * nnc.compute_max_dist(data, n, d);
  unmerged_clusters.max_load_factor(std::numeric_limits<double>::infinity());
  existed = std::vector<bool>(size * 2, false);
  //std::cout << "first min dist " << std::endl;
  min_dist =  nnc.compute_min_dist(unmerged_clusters, existed);
  //min_dist = 1000;
  //std::cout << "MINIMUM distance: " << min_dist << std::endl;

  for (auto && p: unmerged_clusters) lambda.insert(p);
  beta = ceil(log_base_((max_dist/min_dist) * n, 1 + epsilon)); // be carefull four the double / double
  output.reserve(n);
}

/**
* this function takes the centroid of two clusters and computes the new centroid
* based on the following formula.
* new_centroid = ((size_a * mu_a) + (size_b * mu_b))/(size_a + size_b).
*/
void hierarchical_clustering::merge(point &mu_a, point &mu_b, int size_a, int size_b) {
  double den = size_a + size_b;
  double coeff_a = size_a / den; // make sure if this is efficient
  double coeff_b = size_b / den;
  for (int i = 0; i < dimension; ++i)
    merged_cluster[i] = coeff_a * mu_a[i] + coeff_b * mu_b[i];
}

// create a vector for the next round
std::unordered_set<int>  hierarchical_clustering::helper(std::unordered_set<int> &to_merge, double merge_value) {
  unchecked.clear();
  for (auto&& u : to_merge) {
     if(!existed[u]) {
       continue;
     }

     int u_weight = nnc.get_cluster_size(u);
    // assert(u_weight <= size);
     if(u >= 0 && existed[u] && u_weight < size) {
       existed[u] = false;
       //assert(u_weight > 0);

       auto res = nnc.get_point(u);
       nnc.v_delete_cluster(u);
       nnc.delete_cluster(u);
       //assert(u >= 0 && u < 2 * (size + 1) + 1);
       auto t = nnc.query(u, res, u_weight);
       int nn_id    = std::get<0>(t);
       double dist  = std::get<1>(t);
       int t_weight = std::get<2>(t);
    //   std::cout << "nn_id " << nn_id << " u " << u << std::endl;
       if (nn_id >= 0 && dist <= merge_value) {
         auto nn_pt = nnc.get_point(nn_id);
         merge(res, nn_pt, u_weight, t_weight);
         merged_weight = u_weight + t_weight;

         nnc.add_cluster(merged_cluster, merged_weight, last_index); // it should be added to the points as well
         unchecked.insert(last_index);
         lambda.insert(last_index);
         existed[last_index] = true;
         // register the merge operation
         output.push_back(std::make_tuple(std::make_pair(last_index, merged_weight),
                          std::make_pair(u, u_weight), std::make_pair(nn_id, t_weight)));

         last_index = last_index + 1;

         existed[nn_id] = false;
         nnc.v_delete_cluster(nn_id);
         nnc.delete_cluster(nn_id);
         nnc.delete_cluster(u);

         lambda.erase(u);
         lambda.erase(nn_id);
      //   std::cout << "lambda size " << lambda.size() << std::endl;
         if (lambda.size() == 2) {
           int first = 0;
           int sec = 0;
           int i = 0;
           for(auto&& e : lambda) {
             if(i == 0) first = e;
             if(i == 1) sec = e;
             i++;
           }
        //   std::cout << "first " << first << " second " << sec << std::endl;
           int f_weight = nnc.get_cluster_size(first);
           int s_weight = nnc.get_cluster_size(sec);
           existed[first] = false;
           existed[sec] = false;
           merged_weight = f_weight + s_weight;
           output.push_back(std::make_tuple(std::make_pair(last_index, merged_weight),
                            std::make_pair(first, f_weight), std::make_pair(sec, s_weight)));

            lambda.erase(first);
            lambda.erase(sec);
            lambda.insert(last_index);
            unchecked.insert(last_index);
            last_index++;
         }


         if(merged_weight == size) {
           stop = true;
           break;
         }
        } else {
          nnc.v_put_back(u);
          nnc.put_back(res, u);
          existed[u] = true;
        }
        //assert(u_weight <= size);
        if(u_weight >= size) break;
        if(u_weight == 0) break;
      }
   }

   return unchecked;
}
//
void hierarchical_clustering::build_hierarchy() {
  double merge_value;
  //std::cout << "beta " << beta << std::endl;
  for (int i = 0; i <= 2 * beta; ++i) {
    if(lambda.size() < 2) return;
    merge_value = pow(1 + epsilon, i); // find an efficient one
    //std::cout <<"merge value " << merge_value << std::endl;
    auto the_merged_cluster = helper(this->unmerged_clusters, merge_value); // these are the merges
    while (the_merged_cluster.size() > 0) {
      auto tmp = helper(the_merged_cluster, merge_value);
      if(stop) return;
      the_merged_cluster.clear();
      for(auto&& p: tmp) {
        if(existed[p]) the_merged_cluster.insert(p);
      }
   }

   if (stop) return;
   if (the_merged_cluster.size() == 1) unmerged_clusters.insert(*the_merged_cluster.begin());
   for (auto&& m: lambda) unmerged_clusters.insert(m);
  }
}

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
