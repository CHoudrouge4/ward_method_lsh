#include "nn_cluster.h"
#include "hierarchical_clustering.h"
#include <vector>
#include <fstream>
#include <cmath>
#include <time.h>

std::vector<std::vector<double>> read_file(const std::string file_name, int &n, int &m) {
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0.1, 0.5);
   std::ifstream in(file_name);
   int k;
   in >> n >> m;
   std::vector<std::vector<double>> points(n, std::vector<double>(m));
   for (int i = 0; i < n; ++i) {
     for(int j = 0; j < m; ++j) {
       in >> points[i][j];
       //points[i][j] *= 10000;// dis(gen);
     }
   }
   in.close();
   return points;
}

void test_nn_cluster() {
  int n, m;
  std::string file_name = "data.in";
  auto p = read_file(file_name, n, m);
  nnCluster nnc(p, n, m, 0.9, 1, 1, 1);

  std::unordered_set<int> unmerged_clusters;
  std::vector<bool> existed(n * 2, false);
  double dd = nnc.compute_min_dist(unmerged_clusters, existed);
}

double distance(std::vector<double> &p, std::vector<double> &q) {
  double sum = 0.0;
  for (size_t i = 0; i < p.size() && i < q.size(); ++i) {
    sum += (p[i] - q[i]) * (p[i] - q[i]);
  }
  return sum;
}

void compute_matrix_distance(std::vector<std::vector<double>> &data) {
 // for (size_t i = 0; i < data.size(); ++i) {
  //  std::cout << '\t' << i;
 // }
//  std::cout << std::endl;
  double min_dist = 1000;
  for (size_t i = 0; i < data.size(); ++i) {
  //  std::cout << i << '\t';
    for (size_t j = i + 1; j < data.size(); ++j) {
     	if(distance(data[i], data[j]) > 0)
	    min_dist = std::min(distance(data[i], data[j]), min_dist);
    }
//    std::cout << std::endl;
  }
  std::cout << "min distance bf " << min_dist << std::endl;
}

void test_HC(std::string input_file, std::string output_file, double epsilon) {
  std::cout << "test HC" << std::endl;
  int n, d;
  auto data = read_file(input_file, n, d);
  std::cout << "data dimension " << data.size() << ' ' << data[0].size() << std::endl;
  compute_matrix_distance(data);
  int bucket = 2;
  int bins = std::max(2, (int)ceil(std::pow(n, 1/4.0)));
  int run_time = 3 * bins;
  hierarchical_clustering hc(data, n, d, epsilon, bucket, bins, run_time);
  std::cout << "start building" << std::endl;
  clock_t start = clock();
  hc.build_hierarchy();
  clock_t end = clock();
  std::cout << (float)(end - start)/CLOCKS_PER_SEC << std::endl;
  hc.print_file(output_file);
}

int main() {
//  std::vector<std::string> data_names = {"1000", "2000", "5000", "11314"};
  std::vector<std::string> data_names = {"iris", "cancer", "digits", "boston"};
  for (auto&& name: data_names) {
  //  std::string name = "news_" + n + "_100";
  //std::string name = "iris";
    test_HC(name + ".in", name + ".out", 10);
  }
  return 0;
}
