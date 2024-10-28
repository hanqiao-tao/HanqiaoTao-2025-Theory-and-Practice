#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cstdio>
#include <emmintrin.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <stack>
#include <array>
#include <cassert>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ilcplex/ilocplex.h>
#include <numeric>
#include <math.h>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

using namespace std;
#define  BIG_INT 2147483647
#pragma region Defination
typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<NumMatrix>   Num3Matrix;

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<NumVarMatrix>   NumVar3Matrix;
#pragma endregion


#pragma region Global variables
int seed = 922;								 // seed from time for generating             

double intTolerance = 1e-5;                  // solve floating-point precision error
string file_name;                            // instance file name
clock_t starting_time;                       // starting time for the solution procedure
clock_t ending_time;                         // ending time for the solution procedure
int iter_1 = 15;							 // max iteration of (P)LNS 
int iter_2 = 999;                            // max iteration of KP-LS 
int iter_3 = 999;                            // max iteration of BPP-LS 
int upperbound;                              // global upperbound
double constructing_time;					 // time consumption
double heuristic_time;						 // time consumption
int n_pool = 10;							 // number of solutions in the solution pool
int init_ub;
int z_cpx = 0;
int cpx_ub = 999999;
//Input
// Warning: m_flight / K_config must be definitly an integer!!

int num_test;                                // number of test points
int K_config = 2;                            // number of prototype configurations
int m_flight = 36;                            // total number of flights
vector<int> m_k_flight;                      // number of flights per configuration
vector<int> profit;                          // profits of each test point

int T_cap;									 // maximum flight time of each flight
int D_cap;									 // maximum oil level of each flight
vector<int> t_test;                          // time consumption of each test point
vector<int> d_test;                          // oil consumption of each test point
double density;								 // density of the conflict graph

vector<vector<bool>> adj_matr;               // adjacent precedence matrix of the test points
vector<vector<bool>> trans_matr;             // transitive precedence matrix of the test points
vector<pair<int, int>> pre_rela;             // direct precedence relations in pairs, according to adj_matr
vector<vector<int>> pre_sucs;                // all precedence successors of each test point
vector<vector<int>> pre_preds;               // all precedence predecessors of each test point
vector<int> test_early_flight;               // earliest flight of each test point
vector<int> test_late_flight;                // latest flight of each test point

vector<vector<bool>> conflict_matr;          // conflict precedence matrix of the test points
vector<pair<int, int>> conflict_rela;        // direct conflict relations in pairs, according to conflict_matr
vector<vector<int>> conflict_sucs;           // all conflict successors of each test point
vector<vector<int>> conflict_preds;          // all conflict predecessors of each test point

vector<vector<int>> test_config;             // feasible configuration of each test point

int delete_point_num;						 // number of test points deleted in pertubation
// Preprocess
vector<bool> fea_set_map;                    // whether test point can be allocated to the flights existed
vector<int> fea_point;                       // test points can be allocated to the flights existed

// Output
double total_time;							 // total solution time for the heuristic algorithm
double dev;									 // absolute deviation between the incumbent soluiton and upper bound
double gap;									 // percentage gap between the incumbent soluiton and upper bound
vector<int> sol_test;						 // the incumbent soluiton (final solution): store the number of the test point allocated to, or whether the test point is executed
vector<int> sol_flight;						 // the incumbent soluiton (final solution): store the configuration of the flight, or whether the flight is used
string output_file = "yyy";  // NOTE: "yyy" is the path of output file!
//Parallel
vector<vector<vector<int>>> sol_pool;
vector<int> sort_tests;
vector<int> random_list;
int p_min;
#pragma endregion

#pragma region Functions
// functions in heuristic.cpp
double ceil_t(double aa);
double max(double a, double b);
double min(double a, double b);
double dff(double x, int p, double eps);
bool instance_reader_otto(string filename, int & num_test, int & T_cap, vector<int> & t_test, vector<vector<bool>> & pre_matr);
void generate_instance(int & D_cap, vector<int> & d_test, int & m_flight, vector<int> & m_k_flight, vector<vector<bool>> & conflict_matr, vector<vector<int>> & test_config, vector<int> & profit);
void FTSP_closure_warshall(vector<vector<bool>> & trans_matr);
void FTSP_matrix_to_graph_vec(vector<pair<int, int>> & pre_rela, vector<vector<int>> & pre_sucs, vector<vector<int>> & pre_preds, vector<pair<int, int>> & conflict_rela, vector<vector<int>> & conflict_sucs, vector<vector<int>> & conflict_preds);

void FTSP_early_late_flight(vector<int> & test_early_flight, vector<int> & test_late_flight);
void FTSP_preprocessing(vector<bool> & fea_set_map, vector<int> & fea_point);

vector<int> bubble_sort(vector<double> index);
vector<int> greedy_gen_load(int config, vector<int> sort_tests, vector<int> partial_sol_test);
vector<vector<vector<int>>>  BDP(int alpha, int beta);

vector<vector<int>> sub_BPP(vector<vector<int>> init_sol);
vector<vector<int>> local_search(vector<vector<int>> init_sol, int p_min, vector<int> sort_tests, int rand_seed);
vector<vector<int>> SOS(vector<vector<int>> init_sol);
vector<vector<int>> pertubation(vector<vector<int>> sol, int rho, int rand_seed);
void LNS_1(vector<vector<int>>& best_sol);
vector<vector<int>> sol_no1;
vector<vector<int>> sol_no2;

int cpxSolve(int T_solve);
void check(vector<vector<int>> sol);

void WriteLogFile(const string& szString, const string& filename);
void WriteLogFile_byEntry(const string& szString, const string& filename);
void output();
void free_memory();
#pragma endregion

