#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

#include "lsh.h"

using std::map;
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
  points_to_bins_.insert(pair<int, vector<int>>(id, proj));

  for (int i = 0; i < nb_bins_; i++) {
    map<int, vector<int>>::iterator it_bin;
    it_bin = bins_collection_[i].find(proj[i]);
    if (it_bin != bins_collection_[i].end()) {
      (it_bin->second).push_back(id);
    } else {
      vector<int> new_bin;
      new_bin.push_back(id);
      bins_collection_[i].insert(pair<int, vector<int>>(proj[i], new_bin));
    }
  }
}

void LSHDataStructure::RemovePoint(int id){
    map<int,vector<double>>::iterator it;
    it = points_.find(id);
    if(it == points_.end()) return;
    points_.erase(it);

    map<int,vector<int>>::iterator it_to_bins;
    it_to_bins = points_to_bins_.find(id);

    vector < int > proj =  it_to_bins->second;

    for(int i = 0; i < bins_collection_.size(); i++){
	int bin = proj[i];
	map<int,vector<int>>::iterator it_bins;
	it_bins = bins_collection_[i].find(bin);

	// this is slow
	for (int j = 0; j < (it_bins->second).size(); j++){
	    if((it_bins->second)[j] == id){
		(it_bins->second).erase((it_bins->second).begin()+j);
		break;
	    }
	}
    }
}


pair<int, double> LSHDataStructure::QueryPoint(const vector<double>& coordinates,
                                    int running_time) {
  // 1. Get the projection: i.e. a list of bins b_1,...,b_nb_bins
  // 2. Consider the elements in the bins b_1,.., b_nb_bins up to a fixed budget
  // 3. Output the closest one.
  vector<int> proj = Project(coordinates);
  int nb_comparisons = 0;
  int id = points_.begin()->first;
  double min_dist = SqrDist(points_.begin()->second, coordinates);

  for (int i = 0; i < nb_bins_; i++) {
    map<int, vector<int>>::iterator it_bin;
    it_bin = bins_collection_[i].find(proj[i]);
    if (it_bin == bins_collection_[i].end()) continue;
    for (int j = 0; j < (it_bin->second).size(); j++) {
      map<int, vector<double>>::iterator p;
      p = points_.find((it_bin->second)[j]);
      double d = SqrDist(coordinates, p->second);
      if (d < min_dist) {
        min_dist = d;
        id = p->first;
        nb_comparisons++;
      }
      if (nb_comparisons > running_time) return pair<int,double>(id,min_dist);
    }
  }
  return pair<int,double>(id,min_dist);
}

// TODO(cohenaddad): Adding deletion

void LSHDataStructure::Print() {
  for (int i = 0; i < nb_bins_; i++) {
    std::cout << "Hash fun #" << i << std::endl;
    for (std::map<int, vector<int>>::iterator it = bins_collection_[i].begin();
         it != bins_collection_[i].end(); ++it) {
      std::cout << "Bin #" << it->first << ":  ";
      for (int j = 0; j < (it->second).size(); j++) {
        std::cout << (it->second)[j] << "   ";
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
    map<int, vector<int>> new_map;
    bins_collection_.push_back(new_map);
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
//     L.Print();
//
// }
