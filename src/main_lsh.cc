#include "nn_cluster.h"
#include "hierarchical_clustering.h"
#include <vector>
#include <fstream>
#include <cmath>
#include <time.h>

std::vector<std::vector<double>> read_file(const std::string file_name, int &n, int &m) {
   std::ifstream in(file_name);
   int k;
   in >> n >> m;
   std::vector<std::vector<double>> points(n, std::vector<double>(m));
   for (int i = 0; i < n; ++i) {
     for(int j = 0; j < m; ++j) {
       in >> points[i][j];
     }
   }
   in.close();
   return points;
}

void test_nn_cluster() {
  int n, m;
  std::string file_name = "data.in";
  auto p = read_file(file_name, n, m);
  nnCluster nnc(p, n, m, 0.9,1,1,1);

  std::unordered_set<pair_int> unmerged_clusters;
  std::unordered_map<pair_int, bool, pairhash> existed;
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
  for (size_t i = 0; i < data.size(); ++i) {
    std::cout << '\t' << i;
  }
  std::cout << std::endl;
  for (size_t i = 0; i < data.size(); ++i) {
    std::cout << i << '\t';
    for (size_t j = 0; j < data.size(); ++j) {
      std::cout << distance(data[i], data[j]) << '\t';
    }
    std::cout << std::endl;
  }
}

void test_HC(std::string input_file, std::string output_file, double epsilon) {
  std::cout << "test HC" << std::endl;
  int n, d;
  auto data = read_file(input_file, n, d);
  std::cout << "data dimension " << data.size() << ' ' << data[0].size() << std::endl;
  //compute_matrix_distance(data);
  int bucket = 3;
  int bins = 10;//(int)floor(std::pow(n, 1/2.0));
  int run_time = 5 * bins;
  hierarchical_clustering hc(data, n, d, epsilon, bucket, bins, run_time);
  clock_t start = clock();
  hc.build_hierarchy();
  clock_t end = clock();
  std::cout << (float)(end - start)/CLOCKS_PER_SEC << std::endl;
  hc.print_file(output_file);
}

int main() {

  //test_nn_cluster();
  std::vector<std::string> data_names = {"digits"};
  for (auto&& name : data_names) {
    test_HC(name + ".in", name + ".out", 2);
  }
  return 0;
}
