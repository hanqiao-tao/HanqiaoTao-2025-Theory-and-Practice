#include "FTSP-P.h"

int main(void)
{
    vector<int> instance_number;
    for (int i = 0; i <= 2; i++) 
    {
        for (int j = 51; j <= 75; j++) 
        {
            instance_number.push_back(150 * i + j);
        }
    }
    vector<double> dens_vec = { 0.1,0.4};
    for (int den_i = 0; den_i < dens_vec.size(); den_i++)
    {
        density = dens_vec[den_i];
        for (int kk = 0; kk <= 2; kk++)
        {
            //FTSP preparation
            for (int num_id = 0; num_id < instance_number.size(); num_id++)
            {
                K_config = 2 + kk;
                int num = instance_number[num_id];

                bool read_ins = instance_reader_otto("xxx" + to_string(num) + ".alb", num_test, T_cap, t_test, adj_matr);
                // NOTE: "xxx" is the path that stores the data files !!

                cout << "1000_" + to_string(num) + ".alb" << endl;
                generate_instance(D_cap, d_test, m_flight, m_k_flight, conflict_matr, test_config, profit);
                FTSP_closure_warshall(trans_matr);
                FTSP_matrix_to_graph_vec(pre_rela, pre_sucs, pre_preds, conflict_rela, conflict_sucs, conflict_preds);

                //Preprocessing
                FTSP_early_late_flight(test_early_flight, test_late_flight);
                FTSP_preprocessing(fea_set_map, fea_point);

                //Initial solution
                starting_time = clock();
                sol_pool = BDP(3, 2000);
                ending_time = clock();
                constructing_time = (double) (ending_time - starting_time) / CLOCKS_PER_SEC;
                init_ub = sol_pool[0][0][0];
                cout << "the best objective value of the initial solution is: " << sol_pool[0][0][0] << endl;
                printf("the time of constructing solution pool is: %.4f seconds\n", constructing_time);
                vector<vector<int>> best_sol = sol_pool[0];
                upperbound = best_sol[0][0];

                //Local Search
                //First: solve a bin packing problem, to reduce the utilization of the flights and add new test points
                //Second: do neighborhoods, try to increase the objective value
                //Last: Pertubation, remove rho maximum profits tests
        
                //parallel
                for (int i = 0; i < n_pool; i++) 
                {
                    random_list.push_back(i);
                }
                random_shuffle(random_list.begin(), random_list.end());
                random_list[0] = 0;  // for singal thread

                starting_time = clock();
                vector<double> profit_unit_1(num_test + 1);
                for (int i = 1; i <= num_test; i++) 
                {
                    profit_unit_1[i] = min((double)(profit[i]) / (double) (t_test[i]), (double)(profit[i]) / (double) (d_test[i]));
                }
                sort_tests = bubble_sort(profit_unit_1);
                p_min = 100;
                for (int i = 1; i <= num_test; i++) 
                {
                    if (p_min > profit[i]) 
                    {
                        p_min = profit[i];
                    }
                }

                LNS_1(best_sol);

                ending_time = clock();

                //find best solution
                upperbound = best_sol[0][0];

                heuristic_time = (double) (ending_time - starting_time) / CLOCKS_PER_SEC;
                printf("the time of heuristic is: %.4f seconds\n", heuristic_time);
                printf("the objective value is: %d\n", upperbound);
                for (int i = 0; i < best_sol[0].size(); i++) 
                {
                  //cout << best_sol[0][i] << endl;
                }
                for (int i = 0; i < best_sol[1].size(); i++) 
                {
                  //cout << best_sol[1][i] << endl;
                }
                ////check(best_solution);

                // MIP solver
                //z_cpx = cpxSolve(600);
                //cout << "the objective value of the solution solved by cplex is: " << z_cpx << endl;
                //cout << "the upperbound of the instance reported by cplex is: " << cpx_ub << endl;
                //cout << endl;
                output();
                free_memory();

            }
        }
        cout << endl;
    }
    return 0;
}

#pragma region FTSP-P Preparation
double ceil_t(double aa)
{
    double a_d = round(aa);

    if (abs(aa - a_d) <= intTolerance)
    {
        int(a_d);
    }

    double a_ceil;
    if (a_d >= aa - intTolerance)
    {
        a_ceil = double(a_d);
    }
    else
    {
        a_ceil = double(a_d + 1);
    }
    return a_ceil;
}

double max(double a, double b)
{
    if (a > b)
    {
        return a;
    }
    else
    {
        return b;
    }
}

double min(double a, double b)
{
    if (a > b)
    {
        return b;
    }
    else
    {
        return a;
    }
}

// dual feasible function
double dff(double x, int p, double eps)
{
    double u;
    if (x > 1 - eps)
    {
        x = 1.0;
    }
    if (x < eps)
    {
        x = 0.0;
    }
    if (fabs((p + 1) * x - round((p + 1) * x)) < intTolerance)
    {
        u = x;
    }
    else
    {
        u = floor((p + 1) * x) / p;
    }
    return u;
}

bool instance_reader_otto(string filename, int& num_test, int& T_cap, vector<int>& t_test, vector<vector<bool>>& pre_matr) //pre_matr size nj+1 / nj+1
{
    ifstream fin;
    fin.open(filename);
    if (fin.is_open() == true)
    {
        string line_str;
        getline(fin, line_str);
        assert(line_str == "<number of tasks>");
        getline(fin, line_str);
        num_test = stoi(line_str);
        //cout << "number of tasks: " << num_test << endl;

        pre_matr.resize(num_test + 1);
        for (int i = 0; i < pre_matr.size(); i++)
        {
            pre_matr[i].resize(num_test + 1);
        }

        getline(fin, line_str);//emptyline
        getline(fin, line_str);
        assert(line_str == "<cycle time>");
        getline(fin, line_str);
        T_cap = stoi(line_str);
        //cout << "Cycle time: " << T_cap << endl;

        getline(fin, line_str);//emptyline
        getline(fin, line_str);
        assert(line_str == "<order strength>");
        getline(fin, line_str);
        double O_S = stod(line_str);
        //cout << "Order strength of the instance: " << O_S << endl;

        getline(fin, line_str);//emptyline
        getline(fin, line_str);//emptyline
        getline(fin, line_str);
        assert(line_str == "<task times>");

        //read in processing times
        t_test.resize(num_test + 1);
        t_test[0] = 0;
        for (int j = 1; j <= num_test; j++)
        {
            string no_j, dur_j;
            getline(fin, line_str);
            stringstream ss(line_str);
            getline(ss, no_j, ' ');
            int cur_j = stoi(no_j);
            assert(cur_j == j);

            getline(ss, dur_j, ' ');
            t_test[j] = stoi(dur_j);
        }

        //read in precedence relations into the matrix
        getline(fin, line_str);//emptyline
        getline(fin, line_str);
        assert(line_str == "<precedence relations>");
        while (getline(fin, line_str) && line_str.empty() == false)
        {
            string pre_s, suc_s;
            stringstream ss(line_str);
            getline(ss, pre_s, ',');
            int pre = stoi(pre_s);
            getline(ss, suc_s, ',');
            int suc = stoi(suc_s);
            assert(pre < suc);
            pre_matr[pre][suc] = true;
        }
        getline(fin, line_str);
        assert(line_str == "<end>");
        int sum_lower_tri = 0;
        for (int r = 1; r <= num_test; r++)
        {
            for (int c = 0; c <= r; c++)
                sum_lower_tri += pre_matr[r][c];
        }
        assert(sum_lower_tri == 0);

        fin.close();
        return true;
    }
    else
    {
        std::cerr << "Instance reading Error!!!" << endl;
        return false;
    }
}

void generate_instance(int& D_cap, vector<int>& d_test, int& m_flight, vector<int>& m_k_flight, vector<vector<bool>>& conflict_matr, vector<vector<int>>& test_config, vector<int>& profit)
{
    srand(seed);
    // generate oil consumption for each test point i: randomly choose an integer in range [t_test[i] - int(0.5* t_test[i]), t_test[i] + int(0.5* t_test[i]))
    D_cap = T_cap;
    d_test.resize(num_test + 1);
    for (int i = 1; i <= num_test; i++)
    {
        d_test[i] = rand() % t_test[i] + t_test[i] - t_test[i] / 2;
        if (d_test[i] > D_cap)
        {
            d_test[i] = D_cap;
        }
    }

    // detemine number of feasible flights for each configuration
    m_k_flight.resize(K_config + 1);

    for (int k = 1; k <= K_config; k++)
    {
        m_k_flight[k] = m_flight / K_config;
    }

    // generate conflict graph
    int conf_edge = 0;
    conflict_matr.resize(num_test + 1);
    for (int i = 0; i <= num_test; i++)
    {
        conflict_matr[i].resize(num_test + 1);
    }
    while (conf_edge <= (int)(num_test * (num_test - 1)) / 2 * density)
    {
        int j_1 = (rand() % num_test);
        int j_2 = (rand() % (num_test - j_1)) + j_1 + 1;
        if (adj_matr[j_1][j_2])
        {
            continue;
        }
        else
        {
            conflict_matr[j_1][j_2] = true;
            conflict_matr[j_2][j_1] = true;
            conf_edge++;
        }
    }

    // detemine feasible configuration for each test point
    test_config.resize(K_config + 1);
    if (K_config == 1)
    {
        test_config[1].resize(num_test + 1, 1);
    }
    else
    {
        for (int k = 1; k <= K_config; k++)
        {
            test_config[k].resize(num_test + 1);
            for (int j = 1; j <= num_test; j++)
            {
                double r = (double)rand() / RAND_MAX;
                if (r <= 0.8)
                {
                    test_config[k][j] = 1;
                }
            }
        }
        // check if there is a point still have no feasible configuration
        for (int i = 1; i <= num_test; i++)
        {
            int flag_sum = 0;
            for (int k = 1; k <= K_config; k++)
            {
                flag_sum = flag_sum + test_config[k][i];
            }
            if (flag_sum == 0)
            {
                int kkk = rand() % K_config + 1;
                test_config[kkk][i] = true;
            }
        }
    }

    // detemine profit for each test point
    profit.resize(num_test + 1);
    for (int i = 1; i <= num_test; i++)
    {
        profit[i] = (rand() % 51) + 50;
    }
}

//based on the adjacency, find the trans closure, both input and output are 2D matrices
void FTSP_closure_warshall(vector<vector<bool>>& trans_matr)
{
    trans_matr = adj_matr;
    //each step, find with k as intermediate node, can we get from i to j
    for (int k = 0; k <= num_test; k++)
    {
        // Pick all vertices as source one by one 
        for (int i = 0; i <= num_test; i++)
        {
            // Pick all vertices as destination for the 
            // above picked source 
            for (int j = 0; j <= num_test; j++)
            {
                // If vertex k is on a path from i to j, 
                // then make sure that the value of reach[i][j] is 1 
                trans_matr[i][j] = trans_matr[i][j] || (trans_matr[i][k] && trans_matr[k][j]);
            }
        }
    }
}

//transform the matrices into vectors for the algorithms
void FTSP_matrix_to_graph_vec(vector<pair<int, int>>& pre_rela, vector<vector<int>>& pre_sucs, vector<vector<int>>& pre_preds, vector<pair<int, int>>& conflict_rela, vector<vector<int>>& conflict_sucs, vector<vector<int>>& conflict_preds)
{
    // transform precedence matrix
    for (int i = 1; i < num_test; i++)
    {
        for (int j = i + 1; j <= num_test; j++)
        {
            if (adj_matr[i][j] == true)
            {
                pre_rela.push_back(make_pair(i, j));
            }
        }
    }

    pre_sucs.resize(num_test + 1);
    pre_preds.resize(num_test + 1);
    for (int i = 1; i <= num_test; i++)
    {
        for (int j = 1; j <= num_test; j++)
        {
            if (trans_matr[i][j])
            {
                pre_sucs[i].push_back(j);
            }
            if (trans_matr[j][i])
            {
                pre_preds[i].push_back(j);
            }
        }
    }

    // transform conflict matrix
    for (int i = 1; i < num_test; i++)
    {
        for (int j = i + 1; j <= num_test; j++)
        {
            if (conflict_matr[i][j] == true)
            {
                conflict_rela.push_back(make_pair(i, j));
            }
        }
    }

    conflict_sucs.resize(num_test + 1);
    conflict_preds.resize(num_test + 1);
    for (int i = 1; i <= num_test; i++)
    {
        for (int j = 1; j <= num_test; j++)
        {
            if (conflict_matr[i][j])
            {
                conflict_sucs[i].push_back(j);
            }
            if (conflict_matr[j][i])
            {
                conflict_preds[i].push_back(j);
            }
        }
    }
}
#pragma endregion

#pragma region Preprocessing
void FTSP_early_late_flight(vector<int>& test_early_flight, vector<int>& test_late_flight)
{
    test_early_flight.clear();
    test_early_flight.shrink_to_fit();
    test_early_flight.resize(num_test + 1);
    test_late_flight.clear();
    test_late_flight.shrink_to_fit();
    test_late_flight.resize(num_test + 1);

    // time capacity
    for (int j = 1; j <= num_test; j++)
    {
        test_late_flight[j] = m_flight;
        vector<double> pred_job = { double(t_test[j]) / double(T_cap) };
        vector<double> suc_job = { double(t_test[j]) / double(T_cap) };
        vector<double> Eps = { 0.0 };
        //calculate earliest station 

        // dual feasible function
        int sum_t = t_test[j];
        for (int pre = 0; pre < pre_preds[j].size(); pre++)
        {
            pred_job.push_back(double(t_test[pre_preds[j][pre]]) / double(T_cap));
            if (double(t_test[pre_preds[j][pre]]) / double(T_cap) < 0.5 - intTolerance && num_test < 1000)
            {
                Eps.push_back(double(t_test[pre_preds[j][pre]]) / double(T_cap));
            }
        }
        for (int r = 0; r < Eps.size(); r++)
        {
            double eps = Eps[r];
            if (eps < 0.5)
            {
                for (int p = 1; p <= 100; p++)
                {
                    double b = 0;
                    for (int i = 0; i < pred_job.size(); i++)
                    {
                        b = b + dff(pred_job[i], p, eps);
                    }
                    b = ceil_t(b);
                    if (b > test_early_flight[j])
                    {
                        test_early_flight[j] = b;
                    }
                }
            }
        }

        //calculate latest station
        // dual feasible function
        vector<double> Eps_2 = { 0.0 };
        for (int suc = 0; suc < pre_sucs[j].size(); suc++)
        {
            suc_job.push_back(double(t_test[pre_sucs[j][suc]]) / double(T_cap));
            if (double(t_test[pre_sucs[j][suc]]) / double(T_cap) < 0.5 - intTolerance && num_test < 1000)
            {
                Eps_2.push_back(double(t_test[pre_sucs[j][suc]]) / double(T_cap));
            }
        }
        for (int r = 0; r < Eps_2.size(); r++)
        {
            double eps = Eps_2[r];
            if (eps < 0.5)
            {
                for (int p = 1; p <= 100; p++)
                {
                    double b = 0;
                    for (int i = 0; i < suc_job.size(); i++)
                    {
                        b = b + dff(suc_job[i], p, eps);
                    }
                    b = ceil_t(b);
                    if (m_flight + 1 - test_late_flight[j] < b)
                    {
                        test_late_flight[j] = m_flight + 1 - int(ceil_t(b));
                    }
                }
            }
        }
    }
    // oil capacity 
    for (int j = 1; j <= num_test; j++)
    {
        vector<double> pred_job = { double(d_test[j]) / double(D_cap) };
        vector<double> suc_job = { double(d_test[j]) / double(D_cap) };
        vector<double> Eps = { 0.0 };

        // dual feasible function
        int sum_t = d_test[j];
        for (int pre = 0; pre < pre_preds[j].size(); pre++)
        {
            pred_job.push_back(double(d_test[pre_preds[j][pre]]) / double(D_cap));
            if (double(d_test[pre_preds[j][pre]]) / double(D_cap) < 0.5 - intTolerance && num_test < 1000)
            {
                Eps.push_back(double(d_test[pre_preds[j][pre]]) / double(D_cap));
            }
        }
        for (int r = 0; r < Eps.size(); r++)
        {
            double eps = Eps[r];
            if (eps < 0.5)
            {
                for (int p = 1; p <= 100; p++)
                {
                    double b = 0;
                    for (int i = 0; i < pred_job.size(); i++)
                    {
                        b = b + dff(pred_job[i], p, eps);
                    }
                    b = ceil_t(b);
                    if (b > test_early_flight[j])
                    {
                        test_early_flight[j] = b;
                    }
                }
            }
        }

        //calculate latest station
        // dual feasible function
        vector<double> Eps_2 = { 0.0 };
        for (int suc = 0; suc < pre_sucs[j].size(); suc++)
        {
            suc_job.push_back(double(d_test[pre_sucs[j][suc]]) / double(D_cap));
            if (double(d_test[pre_sucs[j][suc]]) / double(D_cap) < 0.5 - intTolerance && num_test < 1000)
            {
                Eps_2.push_back(double(d_test[pre_sucs[j][suc]]) / double(D_cap));
            }
        }
        for (int r = 0; r < Eps_2.size(); r++)
        {
            double eps = Eps_2[r];
            if (eps < 0.5)
            {
                for (int p = 1; p <= 100; p++)
                {
                    double b = 0;
                    for (int i = 0; i < suc_job.size(); i++)
                    {
                        b = b + dff(suc_job[i], p, eps);
                    }
                    b = ceil_t(b);
                    if (m_flight + 1 - test_late_flight[j] < b)
                    {
                        test_late_flight[j] = m_flight + 1 - int(ceil_t(b));
                    }
                }
            }
        }
    }

    // longest path (precedence)
    vector<int> head;         // from 0 to each point i
    head.resize(num_test + 1);
    vector<int> tail;         // from each point i to n + 1
    tail.resize(num_test + 1);
    for (int i = 1; i <= num_test; i++)
    {
        for (int j = 1; j < i; j++)
        {
            if (adj_matr[j][i] && head[i] <= head[j])
            {
                head[i] = head[j] + 1;
            }
        }
    }

    for (int i = num_test; i >= 1; i--)
    {
        for (int j = num_test; j > i; j--)
        {
            if (adj_matr[i][j] && tail[i] <= tail[j])
            {
                tail[i] = tail[j] + 1;
            }
        }
    }

    for (int i = 1; i <= num_test; i++)
    {
        if (head[i] + 1 > test_early_flight[i])
        {
            test_early_flight[i] = head[i] + 1;
        }
        if (m_flight - tail[i] < test_late_flight[i])
        {
            test_late_flight[i] = m_flight - tail[i];
        }
    }
}

void FTSP_preprocessing(vector<bool>& fea_set_map, vector<int>& fea_point)
{
    fea_set_map.resize(num_test + 1);
    for (int i = 1; i <= num_test; i++)
    {
        if (m_flight >= test_early_flight[i])
        {
            fea_set_map[i] = true;
            fea_point.push_back(i);
        }
    }
    delete_point_num = (fea_point.size() / 50 - 1) * 10;
}
#pragma endregion

#pragma region Initial solution
//Bubble sort
vector<int> bubble_sort(vector<double> index)
{
    vector<int> sort_test;
    sort_test.resize(num_test + 1);
    for (int i = 1; i <= num_test; i++)
    {
        sort_test[i] = i;
    }
    int tmp;
    for (int i = 1; i <= num_test; i++)
    {
        for (int j = 1; j <= num_test - i; j++)
        {
            if (index[sort_test[j]] < index[sort_test[j + 1]])
            {
                tmp = sort_test[j];
                sort_test[j] = sort_test[j + 1];
                sort_test[j + 1] = tmp;
            }
        }
    }
    return sort_test;
}

//greedy generate a load for a flight with a certain configuration
vector<int> greedy_gen_load(int config, vector<int> sort_tests, vector<int> partial_sol_test)
{
    vector<bool> feasible_map;
    feasible_map.resize(num_test + 1);
    for (int i = 1; i <= num_test; i++)
    {
        if (partial_sol_test[i] > 0 || fea_set_map[i] == false)
        {
            feasible_map[i] = false;
        }
        else
        {
            feasible_map[i] = true;
        }
    }

    vector<int> tmp_partial_sol;
    int T_c = 0;
    int D_c = 0;

    for (int i_test = 1; i_test <= num_test; i_test++)
    {
        int test = sort_tests[i_test];
        if (feasible_map[test])
        {
            bool flag = true;
            // check configuration constraints
            flag = test_config[config][test];
            // check capacity constraints
            if (flag)
            {
                if (T_c + t_test[test] > T_cap || D_c + d_test[test] > D_cap)
                {
                    flag = false;
                    break;
                }
            }
            // check precedence constraints
            if (flag)
            {
                for (int j = 0; j < pre_preds[test].size(); j++)
                {
                    int test_pred = pre_preds[test][j];
                    if (partial_sol_test[test_pred] <= 0)       // there is a predecessor of current test not executed
                    {
                        flag = false;
                        break;
                    }
                }
            }
            // check conflict constraints
            if (flag)
            {
                for (int j = 0; j < tmp_partial_sol.size(); j++)
                {
                    int test_conflict = tmp_partial_sol[j];
                    if (conflict_matr[test][test_conflict])
                    {
                        flag = false;
                        break;
                    }
                }
            }
            // no violated constraints
            if (flag)
            {
                T_c = T_c + t_test[test];
                D_c = D_c + d_test[test];
                tmp_partial_sol.push_back(test);
            }
        }
    }
    return tmp_partial_sol;
}

//bounded dynamic programming:
//there are two parameters controlling the performance of the dynamic programming method
//alpha: control the number of new state from a current state
//beta:  control the maximum number of states in the current stage
vector<vector<vector<int>>>  BDP(int alpha, int beta)
{
    vector<vector<int>> current_stage_test;
    vector<vector<int>> current_stage_flight;
    vector<vector<int>> next_stage_test;
    vector<vector<int>> next_stage_flight;

    vector<vector<int>> sort_index_list;
    vector<double> profit_max;
    profit_max.resize(num_test + 1);
    vector<double> profit_min;
    profit_min.resize(num_test + 1);
    for (int i = 1; i <= num_test; i++)
    {
        profit_max[i] = max((double)(profit[i]) / (double)(t_test[i]), (double)(profit[i]) / (double)(d_test[i]));
        profit_min[i] = min((double)(profit[i]) / (double)(t_test[i]), (double)(profit[i]) / (double)(d_test[i]));
    }
    sort_index_list.push_back(bubble_sort(profit_max));
    sort_index_list.push_back(bubble_sort(profit_min));

    for (int k = 0; k <= alpha; k++)
    {
        vector<double> unit_profit;
        unit_profit.resize(num_test + 1);
        for (int i = 1; i <= num_test; i++)
        {
            unit_profit[i] = (double)(alpha * profit[i]) / (double)(k * t_test[i] + (alpha - k) * d_test[i]);
        }
        vector<int> sort_tests = bubble_sort(unit_profit);
        sort_index_list.push_back(sort_tests);
    }

    //initialize BDP
    current_stage_test.push_back(vector<int>(num_test + 1));
    current_stage_flight.push_back(vector<int>(m_flight + 1));

    //BDP
    bool if_over = true;
    int num_flight = 1;
    while (if_over)
    {
        if (num_flight > m_flight)
        {
            if_over = false;
            break;
        }
        for (int i = 0; i < current_stage_test.size(); i++)
        {
            vector<int> k_sort;
            vector<int> k_value;
            for (int k = 0; k <= K_config; k++)
            {
                k_sort.push_back(k);
            }
            k_value.resize(K_config + 1);
            for (int j = 1; j <= m_flight; j++)
            {
                if (current_stage_flight[i][j] > 0)
                {
                    k_value[current_stage_flight[i][j]] = k_value[current_stage_flight[i][j]] + 1;
                }
            }
            int tmp;
            for (int j = 1; j <= K_config; j++)
            {
                for (int k = 1; k <= K_config - j; k++)
                {
                    if (k_value[k] > k_value[k + 1])
                    {
                        tmp = k_sort[k];
                        k_sort[k] = k_sort[k + 1];
                        k_sort[k + 1] = tmp;
                    }
                }
            }
            for (int kk = 1; kk <= K_config; kk++)
            {
                int k = k_sort[kk];
                //check whther there is a spare flight with configuration k
                bool if_config = true;
                int config_k_num = 0;
                for (int j = 1; j <= num_flight; j++)
                {
                    if (current_stage_flight[i][j] == k)
                    {
                        config_k_num = config_k_num + 1;
                        if (config_k_num >= m_flight / K_config)
                        {
                            if_config = false;
                            break;
                        }
                    }
                }
                if (if_config)
                {
                    vector<vector<int>> load_set;
                    for (int id_index = 0; id_index < sort_index_list.size(); id_index++)
                    {
                        vector<int> load = greedy_gen_load(k, sort_index_list[id_index], current_stage_test[i]);
                        bool if_same = false;
                        for (int load_id = 0; load_id < load_set.size(); load_id++)
                        {
                            if (load_set[load_id].size() == load.size())
                            {
                                for (int load_test = 0; load_test < load.size(); load_test++)
                                {
                                    if (load[load_test] != load_set[load_id][load_test])
                                    {
                                        break;
                                    }
                                    if_same = true;
                                }
                            }
                            if (if_same)
                            {
                                break;
                            }
                        }
                        if (!if_same)
                        {
                            load_set.push_back(load);
                            if (load.size() == 0)
                            {
                                next_stage_test.push_back(current_stage_test[i]);
                                next_stage_flight.push_back(current_stage_flight[i]);
                            }
                            vector<int> state_test = current_stage_test[i];
                            for (int test_id = 0; test_id < load.size(); test_id++)
                            {
                                state_test[load[test_id]] = num_flight;
                                state_test[0] = state_test[0] + profit[load[test_id]];
                            }
                            next_stage_test.push_back(state_test);
                            vector<int> state_flight = current_stage_flight[i];
                            state_flight[num_flight] = k;
                            next_stage_flight.push_back(state_flight);
                        }

                    }
                }
            }
        }
        // compute and sort objective value of states in next_stage, and delete redundant states
        current_stage_test.clear();
        current_stage_test.shrink_to_fit();
        current_stage_flight.clear();
        current_stage_flight.shrink_to_fit();

        vector<int> obj_map;
        for (int i = 0; i < next_stage_test.size(); i++)
        {
            obj_map.push_back(i);
        }
        // bubble sort
        int tmp;
        for (int i = 0; i < next_stage_test.size(); i++)
        {
            for (int j = 0; j < next_stage_test.size() - i - 1; j++)
            {
                if (next_stage_test[obj_map[j]][0] < next_stage_test[obj_map[j + 1]][0])
                {
                    tmp = obj_map[j];
                    obj_map[j] = obj_map[j + 1];
                    obj_map[j + 1] = tmp;
                }
            }
        }
        // update new stage
        for (int i = 0; i < obj_map.size(); i++)
        {
            if (beta <= i)
            {
                break;
            }
            current_stage_test.push_back(next_stage_test[obj_map[i]]);
            current_stage_flight.push_back(next_stage_flight[obj_map[i]]);
        }
        next_stage_test.clear();
        next_stage_test.shrink_to_fit();
        next_stage_flight.clear();
        next_stage_flight.shrink_to_fit();
        obj_map.clear();
        obj_map.shrink_to_fit();
        num_flight = num_flight + 1;
    }
    vector<vector<vector<int>>> sol_pool;
    vector<vector<int>> init_sol;
    for (int i = 1; i <= n_pool; i++)
    {
        init_sol.push_back(current_stage_test[0]);
        init_sol.push_back(current_stage_flight[0]);
        sol_pool.push_back(init_sol);
        vector<vector<int>>().swap(init_sol);
    }

    return sol_pool;
}
#pragma endregion

#pragma region Local search & Pertubation
vector<vector<int>> sub_BPP(vector<vector<int>> init_sol)
{
    bool improve = true;
    int loop = 0;

    while (improve == true)
    {
        improve = false;
        //find executed test points
        vector<int> test_set;
        vector<vector<int>> flight_cap(m_flight + 1, vector<int>(2));
        for (int i = 1; i <= num_test; i++)
        {
            if (init_sol[0][i] > 0)
            {
                test_set.push_back(i);
                flight_cap[init_sol[0][i]][0] = flight_cap[init_sol[0][i]][0] + t_test[i];
                flight_cap[init_sol[0][i]][1] = flight_cap[init_sol[0][i]][1] + d_test[i];
            }
        }
        // Relocate neighborhood
        if (improve == false)
        {
            for (int i_test = 0; i_test < test_set.size(); i_test++)
            {
                int i = test_set[i_test];
                for (int bin_2 = 1; bin_2 <= m_flight; bin_2++)
                {
                    int bin_1 = init_sol[0][i];
                    int t_1 = flight_cap[bin_1][0];
                    int d_1 = flight_cap[bin_1][1];
                    int t_2 = flight_cap[bin_2][0];
                    int d_2 = flight_cap[bin_2][1];
                    // capacity constraints
                    if (bin_1 != bin_2 && init_sol[1][bin_2] > 0 && t_2 + t_test[i] <= T_cap && d_2 + d_test[i] <= D_cap && t_2 + d_2 + t_test[i] + d_test[i] > t_1 + d_1)
                    {
                        bool check_fea = true;
                        // check precedence constraints
                        if (check_fea)
                        {
                            for (int l = 0; l < pre_sucs[i].size(); l++)
                            {
                                if (init_sol[0][pre_sucs[i][l]] > 0 && init_sol[0][pre_sucs[i][l]] - bin_2 < 1)
                                {
                                    check_fea = false;
                                    break;
                                }
                            }
                            for (int l = 0; l < pre_preds[i].size(); l++)
                            {
                                if (init_sol[0][pre_preds[i][l]] > 0 && bin_2 - init_sol[0][pre_preds[i][l]] < 1)
                                {
                                    check_fea = false;
                                    break;
                                }
                            }
                        }
                        // check conflict constraints
                        if (check_fea)
                        {
                            for (int l = 0; l < conflict_sucs[i].size(); l++)
                            {
                                if (init_sol[0][conflict_sucs[i][l]] == bin_2)
                                {
                                    check_fea = false;
                                    break;
                                }
                            }
                            for (int l = 0; l < conflict_preds[i].size(); l++)
                            {
                                if (bin_2 == init_sol[0][conflict_preds[i][l]])
                                {
                                    check_fea = false;
                                    break;
                                }
                            }
                        }
                        // check configuration constraints
                        check_fea = check_fea && test_config[init_sol[1][bin_2]][i];
                        // all constraints satisfied
                        if (check_fea)
                        {
                            init_sol[0][i] = bin_2;
                            flight_cap[bin_1][0] = flight_cap[bin_1][0] - t_test[i];
                            flight_cap[bin_1][1] = flight_cap[bin_1][1] - d_test[i];
                            flight_cap[bin_2][0] = flight_cap[bin_2][0] + t_test[i];
                            flight_cap[bin_2][1] = flight_cap[bin_2][1] + d_test[i];
                            improve = true;
                        }
                    }
                }
            }
        }
        //check(init_sol);
        // Swap neighborhood
        if (improve == false)
        {
            for (int i_test_1 = 0; i_test_1 < test_set.size(); i_test_1++)
            {
                for (int i_test_2 = i_test_1 + 1; i_test_2 < test_set.size(); i_test_2++)
                {
                    int i_1 = test_set[i_test_1];
                    int i_2 = test_set[i_test_2];
                    int flight_1 = init_sol[0][i_1];
                    int flight_2 = init_sol[0][i_2];
                    int t_1 = flight_cap[flight_1][0];
                    int d_1 = flight_cap[flight_1][1];
                    int t_2 = flight_cap[flight_2][0];
                    int d_2 = flight_cap[flight_2][1];
                    // capacity constraints
                    if (T_cap >= t_1 - t_test[i_1] + t_test[i_2] && T_cap >= t_2 + t_test[i_1] - t_test[i_2] &&
                        D_cap >= d_1 - d_test[i_1] + d_test[i_2] && D_cap >= d_2 + d_test[i_1] - d_test[i_2] &&
                        flight_1 != flight_2)
                    {
                        bool check_fea = true;
                        // check precedence constraints
                        if (check_fea)
                        {
                            for (int l = 0; l < pre_sucs[i_1].size(); l++)
                            {
                                if (init_sol[0][pre_sucs[i_1][l]] > 0 && init_sol[0][pre_sucs[i_1][l]] - flight_2 < 1)
                                {
                                    check_fea = false;
                                    break;
                                }
                            }
                            for (int l = 0; l < pre_preds[i_1].size(); l++)
                            {
                                if (init_sol[0][pre_preds[i_1][l]] > 0 && flight_2 - init_sol[0][pre_preds[i_1][l]] < 1)
                                {
                                    check_fea = false;
                                    break;
                                }
                            }
                            for (int l = 0; l < pre_sucs[i_2].size(); l++)
                            {
                                if (init_sol[0][pre_sucs[i_2][l]] > 0 && init_sol[0][pre_sucs[i_2][l]] - flight_1 < 1)
                                {
                                    check_fea = false;
                                    break;
                                }
                            }
                            for (int l = 0; l < pre_preds[i_2].size(); l++)
                            {
                                if (init_sol[0][pre_preds[i_2][l]] > 0 && flight_1 - init_sol[0][pre_preds[i_2][l]] < 1)
                                {
                                    check_fea = false;
                                    break;
                                }
                            }
                        }
                        // check conflict constraints
                        if (check_fea)
                        {
                            for (int l = 0; l < conflict_sucs[i_1].size(); l++)
                            {
                                if (init_sol[0][conflict_sucs[i_1][l]] == flight_2)
                                {
                                    check_fea = false;
                                    break;
                                }
                            }
                            for (int l = 0; l < conflict_preds[i_1].size(); l++)
                            {
                                if (flight_2 == init_sol[0][conflict_preds[i_1][l]])
                                {
                                    check_fea = false;
                                    break;
                                }
                            }
                            for (int l = 0; l < conflict_sucs[i_2].size(); l++)
                            {
                                if (init_sol[0][conflict_sucs[i_2][l]] == flight_1)
                                {
                                    check_fea = false;
                                    break;
                                }
                            }
                            for (int l = 0; l < conflict_preds[i_2].size(); l++)
                            {
                                if (flight_1 == init_sol[0][conflict_preds[i_2][l]])
                                {
                                    check_fea = false;
                                    break;
                                }
                            }
                        }
                        // check configuration constraints
                        check_fea = check_fea && test_config[init_sol[1][flight_2]][i_1] && test_config[init_sol[1][flight_1]][i_2];
                        // all constraints satisfied
                        if (check_fea)
                        {
                            init_sol[0][i_1] = flight_2;
                            init_sol[0][i_2] = flight_1;
                            flight_cap[flight_1][0] = flight_cap[flight_1][0] - t_test[i_1] + t_test[i_2];
                            flight_cap[flight_1][1] = flight_cap[flight_1][1] - d_test[i_1] + d_test[i_2];
                            flight_cap[flight_2][0] = flight_cap[flight_2][0] - t_test[i_2] + t_test[i_1];
                            flight_cap[flight_2][1] = flight_cap[flight_2][1] - d_test[i_2] + d_test[i_1];
                            improve = true;
                        }
                    }
                }
            }
        }
        //check(init_sol);

        //for (int j = 1; j <= m_flight; j++) 
        //{
        //    if (flights_test[j] == 0) 
        //    {
        //        empty_flights.push_back(j);
        //        init_sol[1][j] = 0;
        //    }
        //}
        //if (empty_flights.size() > 0) 
        //{
        //    for (int i = 1; i <= num_test; i++) 
        //    {
        //        for (int id = 0; id < empty_flights.size(); id++) 
        //        {
        //            if (id < empty_flights.size() - 1 && empty_flights[id] <= init_sol[0][i] && empty_flights[id + 1] > init_sol[0][i])
        //            {
        //                init_sol[0][i] = init_sol[0][i] - 1 - id;
        //            }
        //            if (id == empty_flights.size() - 1 && empty_flights[id] <= init_sol[0][i])
        //            {
        //                init_sol[0][i] = init_sol[0][i] - 1 - id;
        //            }
        //        }
        //    }
        //    for (int id = empty_flights.size() - 1; id >= 0; id--) 
        //    {
        //        for (int j = 1; j < m_flight; j++) 
        //        {
        //            if (empty_flights[id] <= j)
        //            {
        //                init_sol[1][j] = init_sol[1][j + 1];
        //            }
        //        }
        //        init_sol[1][m_flight] = 0;
        //    }
        //}

        if (loop > iter_3)
        {
            break;
        }
        loop = loop + 1;
    }
    return init_sol;
}

vector<vector<int>> local_search(vector<vector<int>> init_sol, int p_min, vector<int> sort_tests, int rand_seed)
{
    srand(rand_seed);
    int best_obj = init_sol[0][0];
    vector<vector<int>> best_sol = init_sol;
    int threshold = 0 - rand() % 21 - p_min;
    vector<int> set_fea(num_test + 1);
    for (int i = 1; i <= num_test; i++)
    {
        if (init_sol[0][i] > 0)
        {
            set_fea[i] = 1;
        }
        if (init_sol[0][i] == 0 && fea_set_map[i])
        {
            set_fea[i] = 2;
        }
    }
    bool improve = true;
    int loop = 0;
    while (improve)
    {
        vector<vector<int>> flight_load(m_flight + 1, vector<int>(2));
        for (int i = 1; i <= num_test; i++)
        {
            if (init_sol[0][i] > 0)
            {
                flight_load[init_sol[0][i]][0] = flight_load[init_sol[0][i]][0] + t_test[i];
                flight_load[init_sol[0][i]][1] = flight_load[init_sol[0][i]][1] + d_test[i];
            }
        }
        improve = false;
        //Insert neighborhood
        if (improve == false)
        {
            for (int i = 1; i <= num_test; i++)
            {
                if (init_sol[0][sort_tests[i]] == 0 && fea_set_map[sort_tests[i]])
                {
                    for (int j = 1; j <= m_flight; j++)
                    {
                        bool flag = true;
                        //capacity constraints
                        if (flight_load[j][0] + t_test[sort_tests[i]] > T_cap || flight_load[j][1] + d_test[sort_tests[i]] > D_cap)
                        {
                            flag = false;
                        }

                        //precedence constraints
                        if (flag)
                        {
                            for (int test_id = 0; test_id < pre_sucs[sort_tests[i]].size(); test_id++)
                            {
                                int suc = pre_sucs[sort_tests[i]][test_id];
                                if (init_sol[0][suc] > 0 && init_sol[0][suc] - j < 1)
                                {
                                    flag = false;
                                    break;
                                }
                            }
                            for (int test_id = 0; test_id < pre_preds[sort_tests[i]].size(); test_id++)
                            {
                                int pred = pre_preds[sort_tests[i]][test_id];
                                if (init_sol[0][pred] > 0 && j - init_sol[0][pred] < 1)
                                {
                                    flag = false;
                                    break;
                                }
                                if (init_sol[0][pred] == 0)
                                {
                                    flag = false;
                                    break;
                                }
                            }
                        }

                        //conflict constraints
                        if (flag)
                        {
                            for (int test_id = 0; test_id < conflict_sucs[sort_tests[i]].size(); test_id++)
                            {
                                if (j == init_sol[0][conflict_sucs[sort_tests[i]][test_id]])
                                {
                                    flag = false;
                                    break;
                                }
                            }
                            for (int test_id = 0; test_id < conflict_preds[sort_tests[i]].size(); test_id++)
                            {
                                if (init_sol[0][conflict_preds[sort_tests[i]][test_id]] == j)
                                {
                                    flag = false;
                                    break;
                                }
                            }
                        }

                        //configuration constraints
                        if (flag && init_sol[1][j] > 0)
                        {
                            flag = flag && test_config[init_sol[1][j]][sort_tests[i]];
                        }

                        //satisfy all constraints, insert test sort_tests[i] into flight j
                        if (flag)
                        {
                            improve = true;
                            init_sol[0][sort_tests[i]] = j;
                            set_fea[sort_tests[i]] = 1;
                            init_sol[0][0] = init_sol[0][0] + profit[sort_tests[i]];
                            best_obj = init_sol[0][0];
                            best_sol = init_sol;
                            flight_load[j][0] = flight_load[j][0] + t_test[sort_tests[i]];
                            flight_load[j][1] = flight_load[j][1] + d_test[sort_tests[i]];
                            if (init_sol[1][j] == 0)           // determine a configuration
                            {
                                // choose a configutaion with least utilization
                                vector<int> config_use(K_config + 1);
                                for (int flight_id = 1; flight_id <= m_flight; flight_id++)
                                {
                                    if (init_sol[1][flight_id] > 0)
                                    {
                                        config_use[init_sol[1][flight_id]] = config_use[init_sol[1][flight_id]] + 1;
                                    }
                                }
                                int config_id;
                                for (int k = 1; k <= K_config; k++)
                                {
                                    if (test_config[k][sort_tests[i]])
                                    {
                                        config_id = k;
                                        break;
                                    }
                                }
                                for (int k = 1; k <= K_config; k++)
                                {
                                    if (config_use[k] < config_use[config_id] && test_config[k][sort_tests[i]])
                                    {
                                        config_id = k;
                                    }
                                }
                                init_sol[1][j] = config_id;
                            }
                            break;
                        }
                    }
                }
            }
        }

        //2-exchange neighborhood
        if (improve == false)
        {
            for (int i_1 = 1; i_1 <= num_test; i_1++)
            {
                if (set_fea[i_1] == 1)
                {
                    vector<vector<int>> moves;
                    int j = init_sol[0][i_1];
                    for (int i_2 = i_1 + 1; i_2 <= num_test; i_2++)
                    {
                        if (set_fea[i_2] == 2 && profit[i_2] - profit[i_1] >= threshold)
                        {
                            bool flag = true;
                            //capacity constraints
                            if (flight_load[j][0] + t_test[i_2] - t_test[i_1] > T_cap || flight_load[j][1] + d_test[i_2] - d_test[i_1] > D_cap)
                            {
                                flag = false;
                            }

                            //precedence constraints
                            if (flag)
                            {
                                for (int test_id = 0; test_id < pre_sucs[i_2].size(); test_id++)
                                {
                                    int suc = pre_sucs[i_2][test_id];
                                    if (init_sol[0][suc] > 0 && init_sol[0][suc] - j < 1)
                                    {
                                        flag = false;
                                        break;
                                    }
                                }
                                for (int test_id = 0; flag && test_id < pre_preds[i_2].size(); test_id++)
                                {
                                    int pred = pre_preds[i_2][test_id];
                                    if (init_sol[0][pred] > 0 && j - init_sol[0][pred] < 1)
                                    {
                                        flag = false;
                                        break;
                                    }
                                    if (init_sol[0][pred] == 0)
                                    {
                                        flag = false;
                                        break;
                                    }
                                }
                                for (int test_id = 0; flag && test_id < pre_sucs[i_1].size(); test_id++)
                                {
                                    if (init_sol[0][pre_sucs[i_1][test_id]] > 0)
                                    {
                                        flag = false;
                                        break;
                                    }
                                }
                            }

                            //conflict constraints
                            if (flag)
                            {
                                for (int test_id = 0; test_id < conflict_sucs[i_2].size(); test_id++)
                                {
                                    int conf = conflict_sucs[i_2][test_id];
                                    if (j == init_sol[0][conf])
                                    {
                                        flag = false;
                                        break;
                                    }
                                }
                                for (int test_id = 0; test_id < conflict_preds[i_2].size(); test_id++)
                                {
                                    int conf = conflict_preds[i_2][test_id];
                                    if (j == init_sol[0][conf])
                                    {
                                        flag = false;
                                        break;
                                    }
                                }
                            }

                            //configuration constraints
                            flag = flag && test_config[init_sol[1][j]][i_2];

                            //satisfy all constraints
                            if (flag)
                            {
                                vector<int> move = { profit[i_2] - profit[i_1], i_2 };
                                moves.push_back(move);
                            }
                        }
                    }
                    if (moves.size() > 0)
                    {
                        int id = 0;
                        int best_p = moves[0][0];
                        for (int ii = 0; ii < moves.size(); ii++)
                        {
                            if (moves[ii][0] > best_p)
                            {
                                best_p = moves[ii][0];
                                id = ii;
                            }
                        }
                        improve = true;
                        int ii_2 = moves[id][1];
                        init_sol[0][i_1] = 0;
                        init_sol[0][ii_2] = j;
                        init_sol[0][0] = init_sol[0][0] + profit[ii_2] - profit[i_1];
                        flight_load[j][0] = flight_load[j][0] + t_test[ii_2] - t_test[i_1];
                        flight_load[j][1] = flight_load[j][1] + d_test[ii_2] - d_test[i_1];
                        set_fea[i_1] = 2;
                        set_fea[ii_2] = 1;
                        if (best_obj < init_sol[0][0])
                        {
                            best_obj = init_sol[0][0];
                            best_sol = init_sol;
                        }
                    }

                }

            }

        }

        //3-exchange neighborhood
        if (improve == false)
        {
            //i_3 to i_1, i_1 to i_2, i_2 drop
            for (int i_1 = 1; i_1 <= num_test; i_1++)
            {
                if (set_fea[i_1] == 1)
                {
                    for (int i_2 = i_1 + 1; i_2 <= num_test; i_2++)
                    {
                        if (set_fea[i_2] == 1)
                        {
                            vector<vector<int>> moves;
                            for (int i_3 = 1; i_3 <= num_test; i_3++)
                            {
                                if (init_sol[0][i_1] != init_sol[0][i_2] && set_fea[i_3] == 2)
                                {
                                    int t_1 = t_test[i_1];
                                    int t_2 = t_test[i_2];
                                    int t_3 = t_test[i_3];
                                    int d_1 = d_test[i_1];
                                    int d_2 = d_test[i_2];
                                    int d_3 = d_test[i_3];
                                    int j_1 = init_sol[0][i_1];
                                    int j_2 = init_sol[0][i_2];
                                    // check capacity
                                    if (profit[i_3] - profit[i_2] >= threshold && flight_load[j_1][0] - t_1 + t_3 <= T_cap && flight_load[j_2][0] - t_2 + t_1 <= T_cap
                                        && flight_load[j_1][1] - d_1 + d_3 <= D_cap && flight_load[j_2][1] - d_2 + d_1 <= D_cap)
                                    {
                                        bool flag = true;
                                        // check configuration
                                        flag = flag && test_config[init_sol[1][j_1]][i_3] && test_config[init_sol[1][j_2]][i_1];
                                        // check conflict
                                        if (flag)
                                        {
                                            for (int test_id = 0; test_id < conflict_sucs[i_1].size(); test_id++)
                                            {
                                                int conf = conflict_sucs[i_1][test_id];
                                                if (j_2 == init_sol[0][conf])
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                            for (int test_id = 0; test_id < conflict_preds[i_1].size(); test_id++)
                                            {
                                                int conf = conflict_preds[i_1][test_id];
                                                if (j_2 == init_sol[0][conf])
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                            for (int test_id = 0; test_id < conflict_sucs[i_3].size(); test_id++)
                                            {
                                                int conf = conflict_sucs[i_3][test_id];
                                                if (j_1 == init_sol[0][conf])
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                            for (int test_id = 0; test_id < conflict_preds[i_3].size(); test_id++)
                                            {
                                                int conf = conflict_preds[i_3][test_id];
                                                if (j_1 == init_sol[0][conf])
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                        }
                                        //check precedence
                                        if (flag)
                                        {
                                            for (int test_id = 0; test_id < pre_sucs[i_1].size(); test_id++)
                                            {
                                                int suc = pre_sucs[i_1][test_id];
                                                if (init_sol[0][suc] > 0 && init_sol[0][suc] - j_2 < 1)
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                            for (int test_id = 0; flag && test_id < pre_preds[i_1].size(); test_id++)
                                            {
                                                int pred = pre_preds[i_1][test_id];
                                                if (init_sol[0][pred] > 0 && j_2 - init_sol[0][pred] < 1)
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                                if (init_sol[0][pred] == 0)
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                            for (int test_id = 0; test_id < pre_sucs[i_3].size(); test_id++)
                                            {
                                                int suc = pre_sucs[i_3][test_id];
                                                if (init_sol[0][suc] > 0 && init_sol[0][suc] - j_1 < 1)
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                            for (int test_id = 0; flag && test_id < pre_preds[i_3].size(); test_id++)
                                            {
                                                int pred = pre_preds[i_3][test_id];
                                                if (init_sol[0][pred] > 0 && j_1 - init_sol[0][pred] < 1)
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                                if (init_sol[0][pred] == 0)
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                            for (int test_id = 0; flag && test_id < pre_sucs[i_2].size(); test_id++)
                                            {
                                                if (init_sol[0][pre_sucs[i_2][test_id]] > 0)
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                        }
                                        // all constraints satisfied
                                        if (flag)
                                        {
                                            vector<int> move = { profit[i_3] - profit[i_2], i_3 };
                                            moves.push_back(move);
                                        }
                                    }

                                }
                            }
                            // find best i_3
                            if (moves.size() > 0)
                            {
                                int best_p = moves[0][0];
                                int id = 0;
                                for (int i = 0; i < moves.size(); i++)
                                {
                                    if (moves[i][0] > best_p)
                                    {
                                        id = i;
                                        best_p = moves[i][0];
                                    }
                                }
                                int ii_3 = moves[id][1];
                                init_sol[0][ii_3] = init_sol[0][i_1];
                                init_sol[0][i_1] = init_sol[0][i_2];
                                init_sol[0][i_2] = 0;
                                init_sol[0][0] = init_sol[0][0] + profit[ii_3] - profit[i_2];
                                flight_load[init_sol[0][ii_3]][0] = flight_load[init_sol[0][ii_3]][0] + t_test[ii_3] - t_test[i_1];
                                flight_load[init_sol[0][ii_3]][1] = flight_load[init_sol[0][ii_3]][1] + d_test[ii_3] - d_test[i_1];
                                flight_load[init_sol[0][i_1]][0] = flight_load[init_sol[0][i_1]][0] + t_test[i_1] - t_test[i_2];
                                flight_load[init_sol[0][i_1]][1] = flight_load[init_sol[0][i_1]][1] + d_test[i_1] - d_test[i_2];
                                set_fea[i_2] = 2;
                                set_fea[ii_3] = 1;
                                if (best_obj < init_sol[0][0])
                                {
                                    best_obj = init_sol[0][0];
                                    best_sol = init_sol;
                                }
                            }
                        }
                    }
                }
            }
        }

        //Drop neighborhood
        if (improve == false)
        {
            for (int id = num_test; id > 0; id--)
            {
                int i = sort_tests[id];
                if (init_sol[0][i] > 0 && init_sol[0][0] - profit[i] >= threshold + best_obj)
                {
                    bool flag = true;
                    for (int j = 0; j < pre_sucs[i].size(); j++)
                    {
                        if (init_sol[0][pre_sucs[i][j]] > 0)
                        {
                            flag = false;
                            break;
                        }
                    }
                    if (flag)
                    {
                        init_sol[0][i] = 0;
                        set_fea[i] = 2;
                        init_sol[0][0] = init_sol[0][0] - profit[i];
                        improve = true;
                    }
                }
            }
        }

        loop = loop + 1;
        if (loop > iter_2)
        {
            break;
        }
    }
    return best_sol;
}

vector<vector<int>> local_search_nothreshold(vector<vector<int>> init_sol, int p_min, vector<int> sort_tests, int rand_seed)
{
    srand(rand_seed);
    int best_obj = init_sol[0][0];
    vector<vector<int>> best_sol = init_sol;
    int threshold = 0 - rand() % 21 - p_min;
    threshold = 0;
    vector<int> set_fea(num_test + 1);
    for (int i = 1; i <= num_test; i++)
    {
        if (init_sol[0][i] > 0)
        {
            set_fea[i] = 1;
        }
        if (init_sol[0][i] == 0 && fea_set_map[i])
        {
            set_fea[i] = 2;
        }
    }
    bool improve = true;
    int loop = 0;
    while (improve)
    {
        vector<vector<int>> flight_load(m_flight + 1, vector<int>(2));
        for (int i = 1; i <= num_test; i++)
        {
            if (init_sol[0][i] > 0)
            {
                flight_load[init_sol[0][i]][0] = flight_load[init_sol[0][i]][0] + t_test[i];
                flight_load[init_sol[0][i]][1] = flight_load[init_sol[0][i]][1] + d_test[i];
            }
        }
        improve = false;
        //Insert neighborhood
        if (improve == false)
        {
            for (int i = 1; i <= num_test; i++)
            {
                if (init_sol[0][sort_tests[i]] == 0 && fea_set_map[sort_tests[i]])
                {
                    for (int j = 1; j <= m_flight; j++)
                    {
                        bool flag = true;
                        //capacity constraints
                        if (flight_load[j][0] + t_test[sort_tests[i]] > T_cap || flight_load[j][1] + d_test[sort_tests[i]] > D_cap)
                        {
                            flag = false;
                        }

                        //precedence constraints
                        if (flag)
                        {
                            for (int test_id = 0; test_id < pre_sucs[sort_tests[i]].size(); test_id++)
                            {
                                int suc = pre_sucs[sort_tests[i]][test_id];
                                if (init_sol[0][suc] > 0 && init_sol[0][suc] - j < 1)
                                {
                                    flag = false;
                                    break;
                                }
                            }
                            for (int test_id = 0; test_id < pre_preds[sort_tests[i]].size(); test_id++)
                            {
                                int pred = pre_preds[sort_tests[i]][test_id];
                                if (init_sol[0][pred] > 0 && j - init_sol[0][pred] < 1)
                                {
                                    flag = false;
                                    break;
                                }
                                if (init_sol[0][pred] == 0)
                                {
                                    flag = false;
                                    break;
                                }
                            }
                        }

                        //conflict constraints
                        if (flag)
                        {
                            for (int test_id = 0; test_id < conflict_sucs[sort_tests[i]].size(); test_id++)
                            {
                                if (j == init_sol[0][conflict_sucs[sort_tests[i]][test_id]])
                                {
                                    flag = false;
                                    break;
                                }
                            }
                            for (int test_id = 0; test_id < conflict_preds[sort_tests[i]].size(); test_id++)
                            {
                                if (init_sol[0][conflict_preds[sort_tests[i]][test_id]] == j)
                                {
                                    flag = false;
                                    break;
                                }
                            }
                        }

                        //configuration constraints
                        if (flag && init_sol[1][j] > 0)
                        {
                            flag = flag && test_config[init_sol[1][j]][sort_tests[i]];
                        }

                        //satisfy all constraints, insert test sort_tests[i] into flight j
                        if (flag)
                        {
                            improve = true;
                            init_sol[0][sort_tests[i]] = j;
                            set_fea[sort_tests[i]] = 1;
                            init_sol[0][0] = init_sol[0][0] + profit[sort_tests[i]];
                            best_obj = init_sol[0][0];
                            best_sol = init_sol;
                            flight_load[j][0] = flight_load[j][0] + t_test[sort_tests[i]];
                            flight_load[j][1] = flight_load[j][1] + d_test[sort_tests[i]];
                            if (init_sol[1][j] == 0)           // determine a configuration
                            {
                                // choose a configutaion with least utilization
                                vector<int> config_use(K_config + 1);
                                for (int flight_id = 1; flight_id <= m_flight; flight_id++)
                                {
                                    if (init_sol[1][flight_id] > 0)
                                    {
                                        config_use[init_sol[1][flight_id]] = config_use[init_sol[1][flight_id]] + 1;
                                    }
                                }
                                int config_id;
                                for (int k = 1; k <= K_config; k++)
                                {
                                    if (test_config[k][sort_tests[i]])
                                    {
                                        config_id = k;
                                        break;
                                    }
                                }
                                for (int k = 1; k <= K_config; k++)
                                {
                                    if (config_use[k] < config_use[config_id] && test_config[k][sort_tests[i]])
                                    {
                                        config_id = k;
                                    }
                                }
                                init_sol[1][j] = config_id;
                            }
                            break;
                        }
                    }
                }
            }
        }

        //2-exchange neighborhood
        if (improve == false)
        {
            for (int i_1 = 1; i_1 <= num_test; i_1++)
            {
                if (set_fea[i_1] == 1)
                {
                    vector<vector<int>> moves;
                    int j = init_sol[0][i_1];
                    for (int i_2 = i_1 + 1; i_2 <= num_test; i_2++)
                    {
                        if (set_fea[i_2] == 2 && profit[i_2] - profit[i_1] >= threshold)
                        {
                            bool flag = true;
                            //capacity constraints
                            if (flight_load[j][0] + t_test[i_2] - t_test[i_1] > T_cap || flight_load[j][1] + d_test[i_2] - d_test[i_1] > D_cap)
                            {
                                flag = false;
                            }

                            //precedence constraints
                            if (flag)
                            {
                                for (int test_id = 0; test_id < pre_sucs[i_2].size(); test_id++)
                                {
                                    int suc = pre_sucs[i_2][test_id];
                                    if (init_sol[0][suc] > 0 && init_sol[0][suc] - j < 1)
                                    {
                                        flag = false;
                                        break;
                                    }
                                }
                                for (int test_id = 0; flag && test_id < pre_preds[i_2].size(); test_id++)
                                {
                                    int pred = pre_preds[i_2][test_id];
                                    if (init_sol[0][pred] > 0 && j - init_sol[0][pred] < 1)
                                    {
                                        flag = false;
                                        break;
                                    }
                                    if (init_sol[0][pred] == 0)
                                    {
                                        flag = false;
                                        break;
                                    }
                                }
                                for (int test_id = 0; flag && test_id < pre_sucs[i_1].size(); test_id++)
                                {
                                    if (init_sol[0][pre_sucs[i_1][test_id]] > 0)
                                    {
                                        flag = false;
                                        break;
                                    }
                                }
                            }

                            //conflict constraints
                            if (flag)
                            {
                                for (int test_id = 0; test_id < conflict_sucs[i_2].size(); test_id++)
                                {
                                    int conf = conflict_sucs[i_2][test_id];
                                    if (j == init_sol[0][conf])
                                    {
                                        flag = false;
                                        break;
                                    }
                                }
                                for (int test_id = 0; test_id < conflict_preds[i_2].size(); test_id++)
                                {
                                    int conf = conflict_preds[i_2][test_id];
                                    if (j == init_sol[0][conf])
                                    {
                                        flag = false;
                                        break;
                                    }
                                }
                            }

                            //configuration constraints
                            flag = flag && test_config[init_sol[1][j]][i_2];

                            //satisfy all constraints
                            if (flag)
                            {
                                vector<int> move = { profit[i_2] - profit[i_1], i_2 };
                                moves.push_back(move);
                            }
                        }
                    }
                    if (moves.size() > 0)
                    {
                        int id = 0;
                        int best_p = moves[0][0];
                        for (int ii = 0; ii < moves.size(); ii++)
                        {
                            if (moves[ii][0] > best_p)
                            {
                                best_p = moves[ii][0];
                                id = ii;
                            }
                        }
                        improve = true;
                        int ii_2 = moves[id][1];
                        init_sol[0][i_1] = 0;
                        init_sol[0][ii_2] = j;
                        init_sol[0][0] = init_sol[0][0] + profit[ii_2] - profit[i_1];
                        flight_load[j][0] = flight_load[j][0] + t_test[ii_2] - t_test[i_1];
                        flight_load[j][1] = flight_load[j][1] + d_test[ii_2] - d_test[i_1];
                        set_fea[i_1] = 2;
                        set_fea[ii_2] = 1;
                        if (best_obj < init_sol[0][0])
                        {
                            best_obj = init_sol[0][0];
                            best_sol = init_sol;
                        }
                    }

                }

            }

        }

        //3-exchange neighborhood
        if (improve == false)
        {
            //i_3 to i_1, i_1 to i_2, i_2 drop
            for (int i_1 = 1; i_1 <= num_test; i_1++)
            {
                if (set_fea[i_1] == 1)
                {
                    for (int i_2 = i_1 + 1; i_2 <= num_test; i_2++)
                    {
                        if (set_fea[i_2] == 1)
                        {
                            vector<vector<int>> moves;
                            for (int i_3 = 1; i_3 <= num_test; i_3++)
                            {
                                if (init_sol[0][i_1] != init_sol[0][i_2] && set_fea[i_3] == 2)
                                {
                                    int t_1 = t_test[i_1];
                                    int t_2 = t_test[i_2];
                                    int t_3 = t_test[i_3];
                                    int d_1 = d_test[i_1];
                                    int d_2 = d_test[i_2];
                                    int d_3 = d_test[i_3];
                                    int j_1 = init_sol[0][i_1];
                                    int j_2 = init_sol[0][i_2];
                                    // check capacity
                                    if (profit[i_3] - profit[i_2] >= threshold && flight_load[j_1][0] - t_1 + t_3 <= T_cap && flight_load[j_2][0] - t_2 + t_1 <= T_cap
                                        && flight_load[j_1][1] - d_1 + d_3 <= D_cap && flight_load[j_2][1] - d_2 + d_1 <= D_cap)
                                    {
                                        bool flag = true;
                                        // check configuration
                                        flag = flag && test_config[init_sol[1][j_1]][i_3] && test_config[init_sol[1][j_2]][i_1];
                                        // check conflict
                                        if (flag)
                                        {
                                            for (int test_id = 0; test_id < conflict_sucs[i_1].size(); test_id++)
                                            {
                                                int conf = conflict_sucs[i_1][test_id];
                                                if (j_2 == init_sol[0][conf])
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                            for (int test_id = 0; test_id < conflict_preds[i_1].size(); test_id++)
                                            {
                                                int conf = conflict_preds[i_1][test_id];
                                                if (j_2 == init_sol[0][conf])
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                            for (int test_id = 0; test_id < conflict_sucs[i_3].size(); test_id++)
                                            {
                                                int conf = conflict_sucs[i_3][test_id];
                                                if (j_1 == init_sol[0][conf])
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                            for (int test_id = 0; test_id < conflict_preds[i_3].size(); test_id++)
                                            {
                                                int conf = conflict_preds[i_3][test_id];
                                                if (j_1 == init_sol[0][conf])
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                        }
                                        //check precedence
                                        if (flag)
                                        {
                                            for (int test_id = 0; test_id < pre_sucs[i_1].size(); test_id++)
                                            {
                                                int suc = pre_sucs[i_1][test_id];
                                                if (init_sol[0][suc] > 0 && init_sol[0][suc] - j_2 < 1)
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                            for (int test_id = 0; flag && test_id < pre_preds[i_1].size(); test_id++)
                                            {
                                                int pred = pre_preds[i_1][test_id];
                                                if (init_sol[0][pred] > 0 && j_2 - init_sol[0][pred] < 1)
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                                if (init_sol[0][pred] == 0)
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                            for (int test_id = 0; test_id < pre_sucs[i_3].size(); test_id++)
                                            {
                                                int suc = pre_sucs[i_3][test_id];
                                                if (init_sol[0][suc] > 0 && init_sol[0][suc] - j_1 < 1)
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                            for (int test_id = 0; flag && test_id < pre_preds[i_3].size(); test_id++)
                                            {
                                                int pred = pre_preds[i_3][test_id];
                                                if (init_sol[0][pred] > 0 && j_1 - init_sol[0][pred] < 1)
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                                if (init_sol[0][pred] == 0)
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                            for (int test_id = 0; flag && test_id < pre_sucs[i_2].size(); test_id++)
                                            {
                                                if (init_sol[0][pre_sucs[i_2][test_id]] > 0)
                                                {
                                                    flag = false;
                                                    break;
                                                }
                                            }
                                        }
                                        // all constraints satisfied
                                        if (flag)
                                        {
                                            vector<int> move = { profit[i_3] - profit[i_2], i_3 };
                                            moves.push_back(move);
                                        }
                                    }

                                }
                            }
                            // find best i_3
                            if (moves.size() > 0)
                            {
                                int best_p = moves[0][0];
                                int id = 0;
                                for (int i = 0; i < moves.size(); i++)
                                {
                                    if (moves[i][0] > best_p)
                                    {
                                        id = i;
                                        best_p = moves[i][0];
                                    }
                                }
                                int ii_3 = moves[id][1];
                                init_sol[0][ii_3] = init_sol[0][i_1];
                                init_sol[0][i_1] = init_sol[0][i_2];
                                init_sol[0][i_2] = 0;
                                init_sol[0][0] = init_sol[0][0] + profit[ii_3] - profit[i_2];
                                flight_load[init_sol[0][ii_3]][0] = flight_load[init_sol[0][ii_3]][0] + t_test[ii_3] - t_test[i_1];
                                flight_load[init_sol[0][ii_3]][1] = flight_load[init_sol[0][ii_3]][1] + d_test[ii_3] - d_test[i_1];
                                flight_load[init_sol[0][i_1]][0] = flight_load[init_sol[0][i_1]][0] + t_test[i_1] - t_test[i_2];
                                flight_load[init_sol[0][i_1]][1] = flight_load[init_sol[0][i_1]][1] + d_test[i_1] - d_test[i_2];
                                set_fea[i_2] = 2;
                                set_fea[ii_3] = 1;
                                if (best_obj < init_sol[0][0])
                                {
                                    best_obj = init_sol[0][0];
                                    best_sol = init_sol;
                                }
                            }
                        }
                    }
                }
            }
        }

        ////Drop neighborhood
        //if (improve == false) 
        //{
        //    for (int id = num_test; id > 0; id--) 
        //    {
        //        int i = sort_tests[id];
        //        if (init_sol[0][i] > 0 && init_sol[0][0] - profit[i] >= threshold + best_obj)
        //        {
        //            bool flag = true;
        //            for (int j = 0; j < pre_sucs[i].size(); j++) 
        //            {
        //                if (init_sol[0][pre_sucs[i][j]] > 0) 
        //                {
        //                    flag = false;
        //                    break;
        //                }
        //            }
        //            if (flag) 
        //            {
        //                init_sol[0][i] = 0;
        //                set_fea[i] = 2;
        //                init_sol[0][0] = init_sol[0][0] - profit[i];
        //                improve = true;
        //            }
        //        }
        //    }
        //}

        loop = loop + 1;
        if (loop > iter_2)
        {
            break;
        }
    }
    return best_sol;
}

vector<vector<int>> pertubation(vector<vector<int>> sol, int rho, int rand_seed)
{
    srand(rand_seed);
    vector<int> set_fea(num_test + 1);
    for (int i = 1; i <= num_test; i++)
    {
        if (sol[0][i] > 0)
        {
            set_fea[i] = 1;
        }
        if (sol[0][i] == 0 && fea_set_map[i])
        {
            set_fea[i] = 2;
        }
    }
    int delete_test = 0;
    //swap configuration
    for (int exchange = 0; exchange <= max(1, 0.3 * (double)m_flight); exchange++)
    {
        int j_1 = rand() % (m_flight - 1) + 1;
        int j_2 = rand() % (m_flight - j_1) + 1 + j_1;
        if (sol[1][j_1] != sol[1][j_2])
        {
            bool flag = true;
            // check
            for (int i = 1; i <= num_test; i++)
            {
                if (!test_config[sol[1][j_1]][i] && sol[0][i] == j_2)
                {
                    //try to drop this test 
                    for (int i_id = 0; i_id < pre_sucs[i].size(); i_id++)
                    {
                        if (sol[0][pre_sucs[i][i_id]] > 0)
                        {
                            flag = false;
                            break;
                        }
                    }
                }
                if (!test_config[sol[1][j_2]][i] && sol[0][i] == j_1)
                {
                    //try to drop this test 
                    for (int i_id = 0; i_id < pre_sucs[i].size(); i_id++)
                    {
                        if (sol[0][pre_sucs[i][i_id]] > 0)
                        {
                            flag = false;
                            break;
                        }
                    }
                }
            }
            if (flag)
            {
                for (int i = 1; i <= num_test; i++)
                {
                    if (!test_config[sol[1][j_2]][i] && sol[0][i] == j_1)
                    {
                        sol[0][i] = 0;
                        set_fea[i] = 2;
                        sol[0][0] = sol[0][0] - profit[i];
                        delete_test = delete_test + 1;
                    }
                    if (!test_config[sol[1][j_1]][i] && sol[0][i] == j_2)
                    {
                        sol[0][i] = 0;
                        set_fea[i] = 2;
                        sol[0][0] = sol[0][0] - profit[i];
                        delete_test = delete_test + 1;
                    }
                }
                int tmp = sol[1][j_1];
                sol[1][j_1] = sol[1][j_2];
                sol[1][j_2] = tmp;
                exchange = exchange + 1;
            }
        }
    }

    //delete flight tests
    vector<double> profit_unit;
    profit_unit.resize(num_test + 1);
    int profit_sum = 0;
    for (int i = 1; i <= num_test; i++)
    {
        profit_unit[i] = (double)(profit[i]) / (t_test[i] + d_test[i]);
    }
    vector<int> sort_tests = bubble_sort(profit_unit);

    for (int id = 1; id <= num_test; id++)
    {
        int i = sort_tests[id];
        double r = (double)rand() / RAND_MAX;
        if (sol[0][i] > 0 && r <= 0.5)
        {
            sol[0][i] = 0;
            sol[0][0] = sol[0][0] - profit[i];
            delete_test = delete_test + 1;
            for (int j = 0; j < pre_sucs[i].size(); j++)
            {
                if (sol[0][pre_sucs[i][j]] > 0)
                {
                    sol[0][pre_sucs[i][j]] = 0;
                    sol[0][0] = sol[0][0] - profit[pre_sucs[i][j]];
                    delete_test = delete_test + 1;
                }
            }
        }
        if (delete_test >= rho)
        {
            break;
        }
    }
    return sol;
}

void LNS_1(vector<vector<int>>& best_sol)
{
    vector<vector<int>> init_sol = sol_pool[0];
    int loop = 0;
    int num_loop = 0;
    int loop_bound_old = 0;
    vector<vector<int>> best_solution = init_sol;
    while (loop < iter_1)
    {
        if (clock() - starting_time > 60 * CLOCKS_PER_SEC)
        {
            //break;
        }
        int num = 0;
        int loop_1 = 0;
        int loop_bound = init_sol[0][0];

        while (loop_1 < 100)
        {
            if (clock() - starting_time > 60 * CLOCKS_PER_SEC)
            {
                break;
            }
            loop_bound = init_sol[0][0];
            init_sol = sub_BPP(init_sol);
            init_sol = local_search(init_sol, p_min, sort_tests, seed);
            //cout << init_sol[0][0] << endl;
            if (loop_bound == init_sol[0][0])
            {
                num = num + 1;
                if (num >= 5)
                {
                    break;
                }
            }
            else
            {
                loop_bound = init_sol[0][0];
                if (best_solution[0][0] < init_sol[0][0])
                {
                    best_solution = init_sol;
                }
                num = 0;
            }
            loop_1 = loop_1 + 1;
        }
        if (loop_bound_old == init_sol[0][0])
        {
            num_loop = num_loop + 1;
        }
        else
        {
            num_loop = 0;
        }
        if (num_loop >= 3)
        {
            break;
        }
        loop_bound_old = init_sol[0][0];
        init_sol = pertubation(init_sol, delete_point_num, seed);
        loop = loop + 1;
    }
    best_sol = best_solution;
}
#pragma endregion

#pragma region Cplex solve
int cpxSolve(int T_solve)
{
    int z = BIG_INT;

    IloEnv env;
    try
    {
        IloModel model(env);
        // Variables
        NumVarMatrix x(env, num_test + 1);
        for (int i = 0; i <= num_test; i++)
        {
            x[i] = IloNumVarArray(env);
            for (int j = 0; j <= m_flight; j++)
            {
                string name = "x" + to_string(i) + to_string(j);
                x[i].add(IloNumVar(env, 0, 1, ILOINT, name.c_str()));
                if (i == 0 || j == 0)
                {
                    x[i][j].setBounds(0, 0);
                }
            }
        }
        NumVarMatrix y(env, m_flight + 1);
        for (int j = 0; j <= m_flight; j++)
        {
            y[j] = IloNumVarArray(env);
            for (int k = 0; k <= K_config; k++)
            {
                string name = "y" + to_string(j) + to_string(k);
                y[j].add(IloNumVar(env, 0, 1, ILOINT, name.c_str()));
                if (k == 0 || j == 0)
                {
                    y[j][k].setBounds(0, 0);
                }
            }
        }

        // Objective
        IloExpr obj(env);
        for (int i = 1; i <= num_test; i++)
        {
            for (int j = 1; j <= m_flight; j++)
            {
                obj += profit[i] * x[i][j];
            }
        }
        //cout << obj << endl;
        model.add(IloMaximize(env, obj));
        obj.end();

        // Constraints
        IloRangeArray constraints(env);
        // time & oil comsuption capacity constraint
        for (int j = 1; j <= m_flight; j++)
        {
            IloExpr sum_xt(env);
            IloExpr sum_xd(env);
            for (int i = 1; i <= num_test; i++)
            {
                sum_xt += t_test[i] * x[i][j];
                sum_xd += d_test[i] * x[i][j];
            }
            for (int k = 1; k <= K_config; k++)
            {
                sum_xt -= T_cap * y[j][k];
                sum_xd -= D_cap * y[j][k];
            }
            constraints.add(sum_xt <= 0);
            constraints.add(sum_xd <= 0);
            //cout << sum_xt << endl << sum_xd << endl;

            sum_xt.end();
            sum_xd.end();
        }

        // Precedence constraints
        for (int pre = 0; pre < pre_rela.size(); pre++)
        {
            x[pre_rela[pre].second][1].setBounds(0, 0);
            for (int j = 1; j <= m_flight - 1; j++)
            {
                IloExpr prece_expr(env);
                for (int k = 1; k <= j; k++)
                {
                    prece_expr += x[pre_rela[pre].first][k];
                    prece_expr -= x[pre_rela[pre].second][k];
                }
                prece_expr -= x[pre_rela[pre].second][j + 1];
                constraints.add(prece_expr >= 0);
                prece_expr.end();
            }
        }

        // Conflict constraints
        for (int clt = 0; clt < conflict_rela.size(); clt++)
        {
            for (int j = 1; j <= m_flight; j++)
            {
                IloExpr clt_cons(env);
                clt_cons += (x[conflict_rela[clt].first][j] + x[conflict_rela[clt].second][j] - 1);
                constraints.add(clt_cons <= 0);
                //cout << clt_cons << endl;
                clt_cons.end();
            }
        }

        // Configuration constraints 1
        for (int i = 1; i <= num_test; i++)
        {
            for (int j = 1; j <= m_flight; j++)
            {
                IloExpr config(env);
                config -= x[i][j];
                for (int k = 1; k <= K_config; k++)
                {
                    if (test_config[k][i] == 1)
                    {
                        config += y[j][k];
                    }
                }
                constraints.add(config >= 0);
                //cout << config << endl;
                config.end();
            }
        }
        // Configuration constraints 2
        for (int j = 1; j <= m_flight; j++)
        {
            for (int k = 1; k <= K_config; k++)
            {
                IloExpr config(env);
                config -= y[j][k];
                for (int i = 1; i <= num_test; i++)
                {
                    if (test_config[k][i] == 1)
                    {
                        config += x[i][j];
                    }
                }
                constraints.add(config >= 0);
                //cout << config << endl;
                config.end();
            }
        }

        // Other constraints
        for (int i = 1; i <= num_test; i++)
        {
            IloExpr ct1(env);
            for (int j = 1; j <= m_flight; j++)
            {
                ct1 += x[i][j];
            }
            constraints.add(ct1 - 1 <= 0);
            //cout << ct1 << endl;
            ct1.end();
        }

        for (int j = 1; j <= m_flight; j++)
        {
            IloExpr ct2(env);
            for (int k = 1; k <= K_config; k++)
            {
                ct2 += y[j][k];
            }
            constraints.add(ct2 - 1 <= 0);
            //cout << ct2 << endl;
            ct2.end();
        }

        for (int k = 1; k <= K_config; k++)
        {
            IloExpr ct3(env);
            for (int j = 1; j <= m_flight; j++)
            {
                ct3 += y[j][k];
            }
            constraints.add(ct3 - m_flight / K_config == 0);
            //cout << ct3 << endl;
            ct3.end();
        }

        for (int j = 1; j <= m_flight; j++)
        {
            IloExpr ct4(env);
            for (int k = 1; k <= K_config; k++)
            {
                ct4 += y[j][k];
            }
            for (int i = 1; i <= num_test; i++)
            {
                ct4 -= x[i][j];
            }
            constraints.add(ct4 <= 0);
            ct4.end();
        }

        // Cut
        for (int i = 1; i <= num_test; i++)
        {
            if (fea_set_map[i])
            {
                for (int j = 1; j <= m_flight; j++)
                {
                    if (j < test_early_flight[i])
                    {
                        //x[i][j].setBounds(0, 0);
                    }
                }
            }
            else
            {
                for (int j = 1; j <= m_flight; j++)
                {
                    //x[i][j].setBounds(0, 0);
                }
            }
        }

        model.add(constraints);
        constraints.end();

        // Solve & output
        IloCplex cplex(model);
        cplex.setParam(IloCplex::TiLim, T_solve);
        cplex.setOut(env.getNullStream());
        cplex.setParam(IloCplex::Param::Threads, 1);
        //cplex.setParam(IloCplex::Param::Preprocessing::Presolve, false);
        //cplex.setParam(IloCplex::Param::MIP::Strategy::RINSHeur, -1);
        cplex.solve();
        z = cplex.getObjValue();
        cpx_ub = int(cplex.getBestObjValue());
        //cout << "the solution solved by cplex is " << cplex.getStatus() << endl;
        //cout << "the global upperbound is " << cplex.getBestObjValue() << endl;
    }
    catch (IloException& e)
    {
        cerr << "Concert exception caught: " << e << endl;
    }
    catch (...)
    {
        cerr << "Unknown exception caught" << endl;
    }
    env.end();

    return z;
}

void check(vector<vector<int>> sol)
{
    // capacity constraints
    for (int j = 1; j <= m_flight; j++)
    {
        int t_c = 0;
        int d_c = 0;
        for (int i = 1; i <= num_test; i++)
        {
            if (sol[0][i] == j)
            {
                t_c = t_c + t_test[i];
                d_c = d_c + d_test[i];
            }
        }
        if (t_c > T_cap || d_c > D_cap)
        {
            cout << "capacity! " << j << endl;
            exit(11);
        }
    }
    // precedence constraints
    for (int i = 0; i < pre_rela.size(); i++)
    {
        int i_1 = pre_rela[i].first;
        int i_2 = pre_rela[i].second;
        if (sol[0][i_2] - sol[0][i_1] < 1 && sol[0][i_2] > 0 && sol[0][i_1] > 0)
        {
            cout << "precedence! " << i_1 << "\t" << i_2 << endl;
            exit(12);
        }
        if (sol[0][i_2] > 0 && sol[0][i_1] == 0)
        {
            cout << "precedence! " << i_1 << "\t" << i_2 << endl;
            exit(12);
        }
    }

    // conflict constraints
    for (int i = 0; i < conflict_rela.size(); i++)
    {
        int i_1 = conflict_rela[i].first;
        int i_2 = conflict_rela[i].second;
        if (sol[0][i_2] - sol[0][i_1] == 0 && sol[0][i_2] >= 1 && sol[0][i_1] >= 1)
        {
            cout << "conflict! " << i_1 << "\t" << i_2 << endl;
            exit(13);
        }
    }

    //configuration constraints
    for (int i = 1; i <= num_test; i++)
    {
        int j = sol[0][i];
        if (j > 0)
        {
            int k = sol[1][j];
            if (test_config[k][i] == false)
            {
                cout << "configuration! " << i << "\t" << sol[0][i] << endl;
                exit(14);
            }
        }
    }
}
#pragma endregion

#pragma region Free memory
void WriteLogFile(const string& szString, const string& filename)
{
    ofstream fout;
    fout.open(filename, std::ios_base::app);
    fout << szString << endl;
    fout.close();
}

void WriteLogFile_byEntry(const string& szString, const string& filename)
{
    ofstream fout;
    fout.open(filename, std::ios_base::app);
    fout << szString << '\t';
    fout.close();
}

void output()
{
    WriteLogFile("", output_file);
    WriteLogFile_byEntry(to_string(K_config), output_file);
    WriteLogFile_byEntry(to_string(density), output_file);
    //WriteLogFile_byEntry(to_string(init_ub), output_file);
    //WriteLogFile_byEntry(to_string(constructing_time), output_file);
    //WriteLogFile_byEntry(to_string(upperbound), output_file);
    //WriteLogFile_byEntry(to_string(heuristic_time), output_file);
    WriteLogFile_byEntry("  ", output_file);
    WriteLogFile_byEntry(to_string(z_cpx), output_file);
    WriteLogFile_byEntry(to_string(cpx_ub), output_file);
    //WriteLogFile_byEntry(to_string(sol_pool[0][0][0]), output_file);
    //WriteLogFile_byEntry(to_string(sol_no1[0][0]), output_file);
    //WriteLogFile_byEntry(to_string(sol_no2[0][0]), output_file);
}

void free_memory()
{
    vector<int>().swap(m_k_flight);
    vector<int>().swap(profit);
    vector<int>().swap(t_test);
    vector<int>().swap(d_test);
    vector<vector<bool>>().swap(adj_matr);
    vector<vector<bool>>().swap(trans_matr);
    vector<pair<int, int>>().swap(pre_rela);
    vector<vector<int>>().swap(pre_sucs);
    vector<vector<int>>().swap(pre_preds);
    vector<int>().swap(test_early_flight);
    vector<int>().swap(test_late_flight);
    vector<vector<bool>>().swap(conflict_matr);
    vector<pair<int, int>>().swap(conflict_rela);
    vector<vector<int>>().swap(conflict_sucs);
    vector<vector<int>>().swap(conflict_preds);
    vector<vector<int>>().swap(test_config);
    vector<bool>().swap(fea_set_map);
    vector<int>().swap(fea_point);
    vector<int>().swap(sol_test);
    vector<int>().swap(sol_flight);
    vector<vector<vector<int>>>().swap(sol_pool);
    vector<int>().swap(sort_tests);
    vector<int>().swap(random_list);
    vector<vector<int>>().swap(sol_no1);
    vector<vector<int>>().swap(sol_no2);
}
#pragma endregion