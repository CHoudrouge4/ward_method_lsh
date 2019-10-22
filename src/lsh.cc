

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

#include <unordered_set>
#include <unordered_map>
#include "lsh.h"


using std::unordered_set;
using std::unordered_map;
using std::max;
using std::min;
using std::pair;
using std::vector;

double LSHDataStructure::SqrDist(const vector<double>& p1,
                                 const vector<double>& p2) {
  double d = 0;
  for (int i = 0; i < p1.size(); i++) {
    d += (p1[i] - p2[i]) * (p1[i] - p2[i]);
  }
  return d;
}

vector<int> LSHDataStructure::Project(const vector<double>& coordinates) {
  vector<int> projections;

  for (int i = 0; i < nb_bins_; i++) {
    int b = projectors_[i].first;
    double c = 0;
    for (int j = 0; j < projectors_[i].second.size(); j++) {
      c += projectors_[i].second[j] * coordinates[j];
    }
    projections.push_back(static_cast<int>((c + b) / r_));
  }

  return projections;
}

void LSHDataStructure::InsertPoint(int id, const vector<double> &coordinates) {
  points_.insert(pair<int, vector<double>>(id, coordinates));

  vector<int> proj = Project(coordinates);
  vector < pair < int, int >> proj_with_id;

  for (int i = 0; i < nb_bins_; i++) {
    unordered_map<int, unordered_set<int>>::iterator it_bin;
    int pos = 0;
    it_bin = bins_collection_[i].find(proj[i]);
    if (it_bin != bins_collection_[i].end()) {
      (it_bin->second).insert(id);
      pos = (it_bin->second).size()-1;
    } else {
	unordered_set<int> new_bin;
	//vector<int> new_bin;
	new_bin.insert(id);
	bins_collection_[i].insert(pair<int, unordered_set<int>>(proj[i], new_bin));
    }
    proj_with_id.push_back(pair<int, int > (proj[i], pos));
  }
  points_to_bins_.insert(pair<int,
			 vector<pair<int, int>>> (id,
						  proj_with_id));
}


// void LSHDataStructure::RemovePoint(int id){
//     unordered_map<int,vector<double>>::iterator it;
//     it = points_.find(id);
//     if(it == points_.end()) {
//       std::cout << "#id " << id << std::endl;
//       return;
//     }
//     points_.erase(id);
//
//     unordered_map<int,vector<pair<int,int>>>::iterator it_to_bins;
//     it_to_bins = points_to_bins_.find(id);
//
//     vector < pair<int,int> > proj =  it_to_bins->second;
//
//     for(int i = 0; i < bins_collection_.size(); i++){
//     	int bin = (proj[i]).first;
//     	int pos = (proj[i]).second;
//     	unordered_map<int,unordered_set<int>>::iterator it_bins;
//     	it_bins = bins_collection_[i].find(bin);
//     	(it_bins->second).erase(id);
//     }
// }


void LSHDataStructure::RemovePoint(int id){
     unordered_map<int,vector<double>>::iterator it;
     it = points_.find(id);
     if(it == points_.end()){
       std::cout << "#id " << id << std::endl;
       return;
     } 
     points_.erase(it);

     unordered_map<int,vector<pair<int,int>>>::iterator it_to_bins;
     it_to_bins = points_to_bins_.find(id);
     if(it_to_bins == points_to_bins_.end()) {
       std::cout << "#id " << id << std::endl;
       return;
     }
     //

     vector < pair<int,int> > proj =  it_to_bins->second;

     for(int i = 0; i < bins_collection_.size(); i++){
        int bin = (proj[i]).first;
      //  std::cout << "Erasing: " << id << " in bin " << bin << std::endl;
        // int pos = (proj[i]).second;
        unordered_map<int,unordered_set<int>>::iterator it_bins;
        it_bins = bins_collection_[i].find(bin);
        (it_bins->second).erase(id);
        // this is slow
//        for (int j = 0; j < (it_bins->second).size(); j++){
//            if((it_bins->second)[j] == id){
//                (it_bins->second).erase((it_bins->second).begin()+j);
//                break;/
//            }
//    }
     }
     points_to_bins_.erase(it_to_bins);
}


pair<int, double> LSHDataStructure::QueryPoint(const vector<double>& coordinates,
                                    int running_time) {
  // 1. Get the projection: i.e. a list of bins b_1,...,b_nb_bins
  // 2. Consider the elements in the bins b_1,.., b_nb_bins up to a fixed budget
  // 3. Output the closest one.
//  if(points_.size() == 0) return {-1, 0.0};
  vector<int> proj = Project(coordinates);
  int nb_comparisons = 0;
  int id = points_.begin()->first;
  double min_dist = SqrDist(points_.begin()->second, coordinates);

  unordered_map<int, unordered_set<int>>::iterator it_bin;
  for (int i = 0; i < nb_bins_; i++) {
    it_bin = bins_collection_[i].find(proj[i]);
    if (it_bin == bins_collection_[i].end()) continue;
    unordered_set<int> myset = (it_bin->second);
    for ( auto it = myset.begin(); it != myset.end(); ++it ){
      unordered_map<int, vector<double>>::iterator p;
      p = points_.find(*it);
      if(p == points_.end()) continue;
      double d = SqrDist(coordinates, p->second);
      if (d < min_dist) {
    	  min_dist = d;
    	  id = p->first;
    	  nb_comparisons++;
      }
      if (nb_comparisons > running_time) return std::make_pair(id, min_dist);
    }
  }
  return std::make_pair(id, min_dist);
}

// TODO(cohenaddad): Adding deletion

void LSHDataStructure::Print() {
  for (int i = 0; i < nb_bins_; i++) {
    std::cout << "Hash fun #" << i << std::endl;
    for (std::unordered_map<int, unordered_set<int>>::iterator it = bins_collection_[i].begin();
         it != bins_collection_[i].end(); ++it) {
      std::cout << "Bin #" << it->first << ":  ";
      for ( auto it2 = (it->second).begin(); it2 != (it->second).end();
	    ++it2 ){
        std::cout << *it2 << "   ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
  }
}

LSHDataStructure::LSHDataStructure(int bucket_size, int nb_bins1,
                                   int dimension) {
  nb_bins_ = nb_bins1;
  r_ = bucket_size;

  std::normal_distribution<double> distrib{0, 1};

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);

  for (int i = 0; i < nb_bins_; i++) {
      srand (time(NULL));
      int offset = rand() % r_ +1;
    vector<double> coordinates;
    coordinates.reserve(dimension);
    for (int j = 0; j < dimension; j++) {
      coordinates.push_back(distrib(generator));
    }
    projectors_.push_back(pair<int, vector<double>>(offset, coordinates));
    unordered_map<int, unordered_set<int>> new_unordered_map;
    bins_collection_.push_back(new_unordered_map);
  }
}


// int main (){
//
//     LSHDataStructure L = LSHDataStructure(5,5,3);
//
//     vector < double > a = vector < double > (3);
//     vector < double > b = vector < double > (3);
//     vector < double > c = vector < double > (3);
//     vector < double > d = vector < double > (3);
//     vector < double > e = vector < double > (3);
//
//     for (int i = 0; i < 3 ; i++){
// 	a[i] = i;
// 	b[i] = i*i;
// 	c[i] = i*i*i;
// 	d[i] = i*i*i*i;
// 	e[i] = i*i*i*i*i;
//     }
//
//     L.InsertPoint(0, a);
//     L.InsertPoint(1, b);
//     L.InsertPoint(2, c);
//     L.InsertPoint(3, d);
//     L.InsertPoint(4, e);
//
//     L.Print();
//
//     L.QueryPoint(b,100);
//
//     L.RemovePoint(3);
//     L.RemovePoint(2);
//     L.RemovePoint(2);
//     L.RemovePoint(4);
//     L.RemovePoint(4);
//     L.RemovePoint(2);
//     L.RemovePoint(4);
//     L.Print();
//
// }
