#include "nn_cluster.h"
#include <vector>
#include <fstream>

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
  nnCluster nnc(p, n, m, 0.9);

  std::unordered_set<pair_int> unmerged_clusters;
  std::unordered_map<pair_int, bool, pairhash> existed;
  double dd = nnc.compute_min_dist(unmerged_clusters, existed);
}

int main() {

  test_nn_cluster();
  return 0;
}
