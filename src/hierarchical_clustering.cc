#include "hierarchical_clustering.h"
#include "utilities.h"
#include <fstream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <tuple>
#include <utility>
#include <sstream>

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

hierarchical_clustering::hierarchical_clustering(std::vector<std::vector<double>> &data, int n, int d, double epsilon_):
                                                               nnc(data, n, d, epsilon_),
                                                               dimension(d),
                                                               size(n),
                                                               epsilon(epsilon_) {
  last_index = data.size();
  max_dist = nnc.compute_max_dist(data, n, d);
  unmerged_clusters.max_load_factor(std::numeric_limits<double>::infinity());
  min_dist = nnc.compute_min_dist(unmerged_clusters, existed);
  std::cout << "the number of unmerged clusters at t0: " << unmerged_clusters.size() << std::endl;
  for (auto && p: unmerged_clusters) lambda.insert(p);
  beta = ceil(log_base_((max_dist/min_dist) * n, 1 + epsilon)); // be carefull four the double / double
  std::cout << "beta (the maximim merge value is) " << beta << std::endl;
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

// create a vector for the next round
std::unordered_set<pair_int>  hierarchical_clustering::helper(std::unordered_set<pair_int> &to_merge, double merge_value) {
  std::unordered_set <pair_int> unchecked; // this one should be placed maybe in different place, maybe it should be in the fields
  std::vector<double> merged_cluster;
//  merged_cluster.reserve(dimension);
  int merged_weight;
  for (auto&& p : to_merge) {
     if(existed[p] && p.w < size) {
       bool entered = false;
       int u = p.id;
       int u_weight = p.w;
       std::cout << "u - u_weight " <<  u << ' ' << u_weight << std::endl;
       if(u_weight == 0) continue;
       to_erase.push_back({u, u_weight});

       auto res = nnc.get_point(u);
       nnc.delete_cluster(u, u_weight); // we delete the cluster because it is the nn of itself

       auto t = nnc.query(res, u_weight);
       double dist  = std::get<1>(t); // getting the distance
       int t_weight = std::get<2>(t);
       std::cout << "comparing distance and merge values " << dist << " ? " << merge_value << std::endl;
       if (dist <= merge_value) {
         entered = true;
         auto nn_pt = nnc.get_point(std::get<0>(t));

         merged_cluster = merge(res, nn_pt, u_weight, t_weight);
         merged_weight = u_weight + t_weight;
         std::cout << ">>>> merged_weight: " << merged_weight << ' ' << u_weight << ' ' << t_weight << std::endl;

        // insert the merged cluster in lambda to check it again in this round, (easy)
        // or we can do while then remove it
        // register the operation


        // add the cluster to the data strcture -- check
        std::cout << "insert the merged cluster with the index: " << last_index << std::endl;
        pair_int mc = std::make_pair(last_index, merged_weight);
        nnc.add_cluster(merged_cluster, merged_weight, last_index); // it should be added to the points as well
        unchecked.insert(mc);
        lambda.insert(mc);
        existed[mc] = true;
        // register the merge operation
        output.push_back(std::make_tuple(nnc.get_index(last_index, merged_weight), nnc.get_index(u, u_weight), nnc.get_index(std::get<0>(t), std::get<2>(t))));
        last_index++;

        existed[p] = false;
        existed[{std::get<0>(t), std::get<2>(t)}] = false;

        to_erase.push_back({std::get<0>(t), std::get<2>(t)});
        nnc.delete_cluster(std::get<0>(t), std::get<2>(t));

        lambda.erase({u, u_weight});
        lambda.erase({std::get<0>(t), std::get<2>(t)});

        if(merged_weight == size) {
          stop = true;
          break;
        }
      } else {
        /**
        * If we did not have a merge in the during the process what we should
        * do?
        * -> reinsert the coordinates of p in the data structure
        * -> is existed still true ?
        * -> should we inserted in lambda ?
        */
       // assert(res_ != nullptr);
       // since we did not mege the cluster we have to take some action

       nnc.add_cluster(res, u_weight, u);
       existed[{u, u_weight}] = true;
       //lambda.insert({u, u_weight});
      }

      // should we have these in this place ?
      if(u_weight >= size) break;
      if(u_weight == 0) break;

      //if(entered) {
        /**
        * add the merged to the next round, I can do it above or here ?!!
        */
        // I already did these actions
        //auto tt = nnc.query(merged_cluster, merged_weight, true);
        //nnc.update_size(idx, std::get<0>(tt), merged_weight);

        //tt = nnc.query(merged_cluster, merged_weight, true);

        //existed[{std::get<0>(tt), merged_weight}] = true;
        // do I need this?
        //unchecked.insert({std::get<0>(tt), merged_weight});
        //nnc.update_dict(std::get<0>(tt), merged_weight, u, u_weight);
      //}

      }
   }

   // deleting the merged clusters
   for(size_t i = 0; i < to_erase.size(); ++i) {
     to_merge.erase(to_erase[i]);
     //lambda.erase(to_erase[i]);
   }
   to_erase.clear();
   return unchecked;
}
//
void hierarchical_clustering::build_hierarchy() {
  /**
  * What is the purpose of this function and how it interact with the helper?
  *
  */
  double merge_value;
  for (int i = 0; i < beta; ++i) {
    merge_value = pow(1 + epsilon, i); // find an efficient one
    auto the_merged_cluster = helper(this->unmerged_clusters, merge_value); // these are the merges
    while (the_merged_cluster.size() > 0) {
      auto tmp = helper(the_merged_cluster, merge_value);
      the_merged_cluster.clear();
      for(auto&& p: tmp) {
        existed[p] = true;
        the_merged_cluster.insert(p);
      }
     if(stop) break;
   }

   if(stop) return;
   if(the_merged_cluster.size() == 1) unmerged_clusters.insert(*the_merged_cluster.begin());
   for(auto m: lambda) unmerged_clusters.insert(m);

  //  lambda.clear();
    if(unmerged_clusters.size() <= 1) break;
  }
}

// return the merges
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
  std::cout << output.size() << std::endl;
  for (auto&& t: output)
    out << toString(std::get<0>(t)) << toString(std::get<1>(t)) << toString(std::get<2>(t)) << std::endl;
  out.close();
}
