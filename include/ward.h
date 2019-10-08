#ifndef WARD_H_
#define WARD_H_


#include "flann/flann.hpp"
#include "hierarchical_clustering.h"
#include <string>
#include <vector>
#include <functional>
#include <unordered_map>
// 
// typedef std::pair<std::string, std::string> pr;
//
//  struct pairhash {
// private:
//    const std::size_t num = 65537;
//    std::hash<std::string> hash_fn;
// public:
//    std::size_t operator()(const std::pair<std::string, std::string> &x) const {
//      return (hash_fn (x.first ) * num) ^ (hash_fn(x.second));
//    }
// };
//
// class ward : public hierarchical_clustering {
// public:
//     ward(const std::string);
//     flann::Matrix<float> dataset;
//     void run();
// private:
//   std::vector<point> M; // dissimilarity matrix,
//   std::unordered_map<pr, double, pairhash> mp;
// };

#endif
