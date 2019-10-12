// /**
// * This file is for testing
// *
// */
// #include <iostream>
// #include <fstream>
// #include <stdlib.h>
// #include <vector>
// #include <random>
// #include <dlfcn.h>
// #include "nn_cluster.h"
// #include "hierarchical_clustering.h"
// #include <time.h>
// #include <limits>
// #include <cmath>
//
// clock_t start;
// clock_t current;
// clock_t elapsed = 0;
//
// inline void begin_record_time() { start = clock(); current = start; }
// inline void get_current_time() { current = clock(); }
// inline clock_t elapsed_time() { return current - start;}
//
// float power(float x, int y) {
//     if (y == 0)
//         return 1;
//     else if (y % 2 == 0)
//         return power(x, y / 2) * power(x, y / 2);
//     else
//         return x * power(x, y / 2) * power(x, y / 2);
// }
//
// float * read_file(const std::string file_name, int &n, int &m) {
//   std::ifstream in(file_name);
//   int k;
//   in >> n >> m;
//   float * array = (float *) malloc(n * m * sizeof(float));
//   for (int i = 0; i < n; ++i) {
//     for(int j = 0; j < m; ++j) {
//       in >> *(array + i*m + j);
//     }
//   }
//   in.close();
//   return array;
// }
//
// float compute_min_distance(float * array, int n, int m) {
//   float result = std::numeric_limits<float>::max();
//   for (int i = 0; i < n; ++i) {
//     for (int j = i + 1; j < n; ++j) {
//       float dist = 0.0;
//       for (int k = 0; k < m; ++k) {
//         float tmp = (*(array + i * m + k)) - (* (array + j * m + k));
//         dist += tmp * tmp;
//       }
//       if(dist < result) result = dist;
//     }
//   }
//   return result;
// }
//
// void print_array (float * array, int n, int m) {
//   std::cout << "[" << ' ';
//   for (int i = 0; i < n; ++i) {
//     std::cout << '[';
//     for(int j = 0; j < m; ++j) {
//       std::cout << *(array + i*m + j) << ' ';
//     }
//     std::cout << ']' << '\n';
//   }
//   std::cout << "]" << '\n';
// }
//
// float * generate_random_matrix(int n, int m) {
//   std::random_device rd;
//   std::mt19937 gen(rd());
//   std::uniform_real_distribution<> dis(-10.0, 10.0);
//   float * array = (float *) malloc(n * m * sizeof(float));
//   for (int i = 0; i < n; ++i) {
//     for (int j = 0; j < m; ++j) {
//       * (array + i * m + j) = dis(gen);
//     }
//   }
//   return array;
// }
//
// void printing_result (const std::vector<std::vector<int>> &indices, const std::vector<std::vector<float>> &dists) {
//   std::cout << "Indices" << '\n';
//   for (auto e: indices) {
//       for (auto ee: e) std::cout << ee << '\n';
//   }
//   std::cout << "Distance" << '\n';
//   for (auto d_: dists) {
//       for (auto d: d_) std::cout << d << '\n';
//   }
// }
//
// void print_matrix (const flann::Matrix<float> &dataset, const int n, const int d) {
//   for (int i = 0; i < n; ++i) {
//      for (int j = 0; j < d; ++j) {
//        std::cout << dataset[i][j] << ' ';
//      }
//      std::cout << '\n';
//   }
// }
//
// extern "C" typedef double (*func_t)(int n, int d, void * array);
//
// void test_add_delete_cluster(nnCluster &index, const int n, const int d) {
//
//   std::random_device rd;
//   std::mt19937 gen(rd());
//   int n_d_s = index.get_number_of_data_structures();
//   std::uniform_int_distribution<> dis(1, n_d_s);
//   begin_record_time();
//   while(elapsed_time() < 100000000) {
//     float * new_cluster_ = generate_random_matrix(1, d);
//     int random_size = dis(gen);
//     flann::Matrix<float> new_cluster(new_cluster_, 1, d);
//
//     std::cout << ">>>> the random size is " << random_size << std::endl;;
//     index.add_new_cluster(new_cluster, random_size);
//
//
//     auto t = index.query(new_cluster, random_size, true);
//     std::cout << "query result " << std::get<0>(t) << ' ' << std::get<1>(t) << ' ' << std::get<2>(t) << '\n';
//     float * res = index.get_point(std::get<0>(t), std::get<2>(t));
//     float dist = 0;
//     for (int i = 0; i < d; ++i) {
//       std::cout << new_cluster[0][i] << ' ' << *(new_cluster_ + i) << ' ' << *(res + i) << std::endl;
//       float tmp = (*(new_cluster_ + i) -  (*(res + i)));
//       dist += tmp * tmp;
//     }
//     std::cout << "distance "<< dist << '\n';
//     for(int i = 0; i < d; ++i) {
//       assert(*(new_cluster_ + i) == *(res + i));
//     }
//     get_current_time();
//   }
// }
//
// void test_data_structure() {
//   std::cout << "testing the data structure" << '\n';
//   int n = 5;
//   int d = 2;
//   float * points = generate_random_matrix(n, d);
//   double epsilon = 0.5;
//   double gamma = 0.9;
//   nnCluster index(points, n, d, epsilon, gamma, 16, 4);
//   //index.compute_min_dist();
//   //test_add_delete_cluster(index, n, d);
// }
//
// void test_HC() {
//   int n;
//   int d;
//   //float * points = generate_random_matrix( n, d);
//   std::vector<std::string> data = {"./data/iris", "./data/cancer", "./data/digits", "./data/boston"};
//
//   const int trees = 2;
//   const int leaves = 10;
//   //std::vector<std::string> data = {"./data/boston.in"};
//   float epsilon;
//   std::cin >> epsilon;
//   std::cout << "epsilon read" << std::endl;
//   for(auto&& data_name : data) {
//   //  std::string data_name = "./data/data10000_0_100";
//     float * points = read_file(data_name + ".in", n , d);
//     std::cout << "done reading" << std::endl;
//
//     hierarchical_clustering hc(points, n, d, epsilon, 0.9, trees, leaves);
//     std::cout << "done initializing" << std::endl;
//     clock_t start = clock();
//     hc.build_hierarchy();
//     clock_t end = clock();
//     std::cout << (float)(end - start)/CLOCKS_PER_SEC << std::endl;
//     epsilon = epsilon * 100;
//     std::string output_file = data_name + std::to_string((int)floor((epsilon))) + "_" + std::to_string(trees) + "_"  + std::to_string(leaves) + ".out";
//     std::cout << output_file << std::endl;
//     hc.print_file(output_file);
//     epsilon /= 100;
//   }
// }
//
// void the_big_exp() {
//   int n;
//   int d = 20;
//   int k = 10;
//   std::ofstream out("perfs7.txt", std::ios_base::app);
//   std::vector<int> trees = {2};
//   std::vector<int> leaves = {10};
//   std::vector<float> epsilons = {8};
//   for (auto&& e: epsilons) {
//     for (auto&& tr: trees) {
//       for(auto&& l: leaves) {
//          for (int i = 10000; i < 20000; i += 1000) {
//             out << e << ' ' << tr << ' ' << l << ' ' << i << ' ' << d;
//             for (int j = 0; j < 1; ++j) {
//
//         //      if(i == 18000) {j = 1;}
//               std::string data_name = "data" + std::to_string(i) + '_' + std::to_string(j) + '_' + std::to_string(d) + '_' + std::to_string(k);
//               std::string file_name = "./data/" + data_name + ".in";
//               float epsilon = e * 100;
//               std::string output_file = data_name + '_' + std::to_string((int)floor((epsilon))) + "_" + std::to_string(tr) + "_"  + std::to_string(l) + ".out";
//               std::cout << output_file << std::endl;
//
//               float * points = read_file(file_name, n, d);
//               //std::cout << "done reading" << std::endl;
//               hierarchical_clustering hc(points, i, d, e, 0.9, tr, l);
//               //std::cout << "done initializing" << std::endl;
//               clock_t start = clock();
//               hc.build_hierarchy();
//               clock_t end = clock();
//               out << ' ' << (float)(end - start)/CLOCKS_PER_SEC;
//               std::cout << (float)(end - start)/CLOCKS_PER_SEC << std::endl;
//               hc.print_file(output_file);
//           }
//           out << std::endl;
//         }
//       }
//     }
//   }
//
//   out.close();
// }
//
// float logb(float num, float base) {
//   return std::log10(num)/ std::log10(base);
// }
//
// void test_power_log(const int n, const float epsilon) {
//   int number_of_ds = floor(logb(n , 1 + epsilon)) + 1;
//   for (int i = 1; i < n; ++i)  {
//     int pow_res = power(1 + epsilon, i);
//     int log_res = floor(logb(i , 1 + epsilon));
//     std::cout << n << ' ' << i << ' ' << log_res << '/' << number_of_ds << '\n';
//     assert(log_res >= 0);
//     assert(log_res < number_of_ds);
//   }
// }
//
// void test_power_log_n(const float epsilon) {
//   for (int n = 1; n < 10000; ++n) test_power_log(n, epsilon);
// }
//
// void test_news_group(std::string file_name) {
//   std::cout << "testing news group" << std::endl;
//   int n;
//   int d;
//
//   // std::string file_name = "./news_1000.in";
//   std::ofstream out("epsilons_newsgroup_perfs.txt", std::ios_base::app);
//   std::vector<int> trees = {16};
//   std::vector<int> leaves = {32};
//   std::vector<float> epsilons = {10};
//   for (auto&& e: epsilons) {
//     for (auto&& tr: trees) {
//       for (auto&& l: leaves) {
//         out << e << ' ' << tr << ' ' << l << ' ' << n << ' ' << d;
//
// 	float * points = read_file(file_name, n, d);
// 	std::cout << n << ' ' << d << std::endl;
//         hierarchical_clustering hc(points, n, d, e, 0.9, tr, l);
//
//         // building HC
//         clock_t start = clock();
//         hc.build_hierarchy();
//         clock_t end = clock();
//
//         out << ' ' << (float)(end - start)/CLOCKS_PER_SEC;
//         std::cout << (float)(end - start)/CLOCKS_PER_SEC << std::endl;
//
//         std::string data_name = "news" + std::to_string(n) + '_' + std::to_string(d);
//         float epsilon = e * 100;
//         std::string output_file = data_name + '_' + std::to_string((int)floor((epsilon))) + "_" + std::to_string(tr) + "_"  + std::to_string(l) + ".out";
//         std::cout << output_file << std::endl;
// 		    hc.print_file(output_file);
// 		    out << std::endl;
//       }
//     }
//   }
//   out.close();
// }
//
// void test_news_group_several_sizes() {
// 	int sizes = 100;
// 	while (sizes < 11314) {
// 		std::string file_name = "./news_" + std::to_string(sizes) + ".in";
// 		test_news_group(file_name);
// 		sizes *= 2;
// 	}
// }
//
//int main () {
  //test_data_structure();
  //test_HC();
  //the_big_exp();
  //int n = 100;
  // float epsilon = 0.5f;
  // //test_power_log(n, epsilon);
  // for (int i = 0; i < 100; ++i) {
  //   epsilon *= 2;
  //   test_power_log_n(epsilon);
  // }
  //
  //std::string file_name = "./news_11314.in";
  //test_news_group(file_name);
  //test_news_group_several_sizes();
//   std::cout << "yeyy" << std::endl;
// 	return 0;
// }
