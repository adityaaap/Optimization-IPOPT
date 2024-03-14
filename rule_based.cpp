// g++ -std=c++11 -O2 -o rule_based rule_based.cpp -lcasadi

#include <casadi/casadi.hpp>
#include <iostream>
#include <vector>
#include <cmath>

using namespace casadi;
using namespace std;

// Function createNlConstraints(const vector<double>& init_pose_e,
//                              const vector<double>& ref_traj_x,
//                              const vector<double>& ref_traj_y,
//                              const vector<double>& ref_traj_theta,
//                              const vector<double>& pedestrian_position,
//                              double T, int N, const vector<vector<int>>& enforced_rules,
//                              double dist_from_ped) {
//     // Define control and slack variables
//     MX v = MX::sym("v", N);
//     MX omega = MX::sym("omega", N);
//     MX delta_e = MX::sym("delta_e", N);
//     MX slack_cbf_1 = MX::sym("slack_cbf_1", N);
//     MX slack_cbf_2 = MX::sym("slack_cbf_2", N);
//     MX slack_cbf_3 = MX::sym("slack_cbf_3", N);
    
//     // Vehicle dynamics placeholders (need defining with actual vehicle state update equations)
//     MX x = MX::sym("x", N + 1);
//     MX y = MX::sym("y", N + 1);
//     MX theta = MX::sym("theta", N + 1);
    
//     // Other constants and placeholders
//     double gamma1 = 0.1;
//     double gamma2 = 100;
//     double gamma3 = 1;
    
//     // Container for inequality constraints
//     std::vector<MX> ineq_constraint;
    
//     // Add dynamic constraints and other inequality constraints
//     for (int k = 0; k < N; ++k) {
//         // Dynamics
//         ineq_constraint.push_back(x[k + 1] - (x[k] + v[k] * cos(theta[k]) * T));
//         ineq_constraint.push_back(y[k + 1] - (y[k] + v[k] * sin(theta[k]) * T));
//         ineq_constraint.push_back(theta[k + 1] - (theta[k] + omega[k] * T));
        
//         // Add other inequality constraints from MATLAB code
//         // ...
//     }
    
//     // Container for equality constraints
//     std::vector<MX> eq_constraint;
    
//     // Enforce end pose (add other constraints from the for loop)
//     eq_constraint.push_back(v[N - 1]);
//     eq_constraint.push_back(omega[N - 1]);
  
//     // Apply enforcement rules, if any
//     for (auto rule : enforced_rules) {
//         if (rule[0] == 1) {
//             for (int i = 0; i < N; ++i) {
//                 eq_constraint.push_back(slack_cbf_1[i]);
//             }
//         }
//         // Add checks for other rules
//     }

//     // Concatenate all constraints
//     MX ineq = vertcat(ineq_constraint);
//     MX eq = vertcat(eq_constraint);
    
//     // Return as a CasADi function
//     return Function("nl_constraints", {v, omega, delta_e, slack_cbf_1, slack_cbf_2, slack_cbf_3}, {ineq, eq});
// }

// Helper function to convert degrees to radians
double deg2rad(double degrees) {
    return degrees * M_PI / 180.0;
}

MX objective_function(MX &v, MX &omega, MX &delta_e, MX &slack_cbf_1, MX &slack_cbf_2, MX &slack_cbf_3, int N){
    // MX v = u(Slice(0 * N, 1 * N));
    // MX omega = u(Slice(1 * N, 2 * N));
    // MX delta_e = u(Slice(2 * N, 3 * N));
    // MX slack_cbf_1 = u(Slice(3 * N, 4 * N));
    // MX slack_cbf_2 = u(Slice(4 * N, 5 * N));
    // MX slack_cbf_3 = u(Slice(5 * N, 6 * N));

    /*************************************************/
    // MX v = v; //u(Slice(0, N));
    // MX omega = omega; //u(Slice(N, 2*N));
    // MX delta_e = delta_e; //u(Slice(2*N, 3*N));
    // MX slack_cbf_1 = slack_cbf_1; //u(Slice(3*N, 4*N));
    // MX slack_cbf_2 = slack_cbf_2; //u(Slice(4*N, 5*N));
    // MX slack_cbf_3 = slack_cbf_3; //u(Slice(5*N, 6*N));

    // Tuning parameters
    double cost_coeff_v = 0.1;
    double cost_coeff_omega = 0.1;
    double cost_coeff_clf_slack = 0.01;
    double cost_coeff_cbf_slack_1 = 100;
    double cost_coeff_cbf_slack_2 = 1.5;
    double cost_coeff_cbf_slack_3 = 0.1;

    MX cost = cost_coeff_v * sumsqr(v) +
              cost_coeff_omega * sumsqr(omega) +
              cost_coeff_clf_slack * sumsqr(delta_e) +
              cost_coeff_cbf_slack_1 * sumsqr(slack_cbf_1) +
              cost_coeff_cbf_slack_2 * sumsqr(slack_cbf_2) +
              cost_coeff_cbf_slack_3 * sumsqr(slack_cbf_3);
    cout<<"obj end"<<endl;
    return cost;
}

void nl_constraints(Opti &opti, MX &v, MX &omega, MX &delta_e, MX &slack_cbf_1, MX &slack_cbf_2, MX &slack_cbf_3, vector<double> init_pose_e, vector<double> ref_traj_x, vector<double> ref_traj_y, vector<double> ref_traj_theta, vector<double> pedestrian_position, double N, double T, vector<vector<int>> enforced_rules, double dist_from_ped){
    // Extract control variables from `u`
    //MX vari = u(0);
    // MX v_store = u(0); //u(Slice(0, N));
    // MX omega_store = u(1); //u(Slice(N, 2*N));
    // MX delta_e_store = u(2); //u(Slice(2*N, 3*N));
    // MX slack_cbf_1_store = u(3); //u(Slice(3*N, 4*N));
    // MX slack_cbf_2_store = u(4); //u(Slice(4*N, 5*N));
    // MX slack_cbf_3_store = u(5); //u(Slice(5*N, 6*N));
    // MX v_store = u(Slice(0 * N, 1 * N));
    // MX omega = u(Slice(1 * N, 2 * N));
    // MX delta_e = u(Slice(2 * N, 3 * N));
    // MX slack_cbf_1 = u(Slice(3 * N, 4 * N));
    // MX slack_cbf_2 = u(Slice(4 * N, 5 * N));
    // MX slack_cbf_3 = u(Slice(5 * N, 6 * N));
    MX v_store = v;
    MX omega_store = omega;
    MX delta_e_store = delta_e;
    MX slack_cbf_1_store = slack_cbf_1;
    MX slack_cbf_2_store = slack_cbf_2;
    MX slack_cbf_3_store = slack_cbf_3;

    // Constants
    double eps = 1.5; // tuning constant for CLF
    double gamma1 = 0.1; // tuning constant for CBF 1
    double gamma2 = 100; // tuning constant for CBF 2
    double gamma3 = 1; // tuning constant for CBF 3

    // Define state variables
    MX x = MX::zeros(N,1);
    MX y = MX::zeros(N,1);
    MX theta = MX::zeros(N,1);

    // Set initial state
    x(0) = init_pose_e[0];
    y(0) = init_pose_e[1];
    theta(0) = init_pose_e[2];

    // System dynamics and state update, use map to simulate for-loop behavior
    // MXDict f = {{"x", x}, {"y", y}, {"theta", theta}, {"v", v_store}, {"omega", omega_store}};
    // auto sys_dynamics = [&](const MXDict& arg) {
    //     return MXDict{{"x_next", arg.at("x") + arg.at("v")*cos(arg.at("theta"))*T},
    //                   {"y_next", arg.at("y") + arg.at("v")*sin(arg.at("theta"))*T},
    //                   {"theta_next", arg.at("theta") + arg.at("omega")*T}};
    // };
    // Function dynamics_func("dynamics_func", f, sys_dynamics(f));
    // Map dynamics_map("dynamics_map", dynamics_func, N);
    // cout<<N<<endl;
    // cout<<x(0) << y(0)<<theta(0)<<endl;
    // cout<<v_store<<endl;

    for (int i = 0; i < N-1; i++) {
        x(i+1) = x(i) + v_store(i) * cos(theta(i)) * T;
        y(i+1) = y(i) + v_store(i) * sin(theta(i)) * T;
        theta(i+1) = theta(i) + omega_store(i) * T;
    }
    //cout<<"states initiated"<<endl;

    // Apply the map function for state update
    // MXDict result = dynamics_map(f);
    // x = vertcat(x(0), result.at("x_next"));
    // y = vertcat(y(0), result.at("y_next"));
    // theta = vertcat(theta(0), result.at("theta_next"));

    // Error calculation for CLF constraint
    MX error = pow(x - ref_traj_x, 2) + pow(y - ref_traj_y, 2) + pow(theta - ref_traj_theta, 2);

    // Inequality constraints
    opti.subject_to(v_store * (2*(x - ref_traj_x) * cos(theta) + 2*(y - ref_traj_y) * sin(theta) +
                     2*(theta - ref_traj_theta) * omega_store) + eps*error - delta_e_store <= 0);
    
    
    opti.subject_to(delta_e_store >= 0);
    opti.subject_to(v_store >= 0);
    opti.subject_to(v_store <= 20);
    opti.subject_to(omega_store >= -deg2rad(50));
    opti.subject_to(omega_store <= deg2rad(50));
    
    // opti.subject_to(pow(pedestrian_position[0] - x, 2) + pow(pedestrian_position[1] - y, 2) -
    //                  gamma1*(2*(x - pedestrian_position[0])*v_store*cos(theta) +
    //                  2*(y - pedestrian_position[1])*v_store*sin(theta)) - pow(dist_from_ped, 2) +
    //                  slack_cbf_1_store >= 0);
    opti.subject_to(-v_store * cos(theta) + gamma2*(-y + 5) - slack_cbf_2_store <= 0);
    opti.subject_to(gamma3*(v_store - 5) - slack_cbf_3_store >= 0);
    opti.subject_to(y <= 0.75);
    cout<<"before equality const"<<endl;
    // Equality constraints
    //cout<<x(N-1)<<endl;
    //cout<<"first ineq error"<<endl;
    double x_end = ref_traj_x[N-1];
    opti.subject_to(x(N-1) == x_end);
    cout<<"first ineq error"<<endl;
    opti.subject_to(y(N-1) == ref_traj_y[N-1]);
    opti.subject_to(theta(N-1) == ref_traj_theta[N-1]);
    opti.subject_to(v_store(N-1) == 0);
    opti.subject_to(omega_store(N-1) == 0);
    
    // Enforce CBF slack variables if rules are applied
    if (enforced_rules[0][0] == 1) {
        opti.subject_to(slack_cbf_1_store == 0);
    }
    if (enforced_rules[0][1] == 1) {
        opti.subject_to(slack_cbf_2_store == 0);
    }
    if (enforced_rules[0][2] == 1) {
        opti.subject_to(slack_cbf_3_store == 0);
    }
}


int main() {
    // Parameters
    //cout<< '0' - '0' <<endl;
    double a = 1;
    double L = 2;
    double tf = 2;
    double T = 0.05;
    int N = static_cast<int>(tf / T);
    vector<double> init_pose_e = {0, 2, 0};
    vector<double> pedestrian_position = {10, 1.75};
    double dist_from_ped = 1;
    vector<vector<int>> enforced_rules = {{1,1,0}}; // Just one for the purpose of this code

    // Reference trajectory
    vector<double> ref_traj_x(N);
    vector<double> ref_traj_y(N, 1.75); // same value for all N
    vector<double> ref_traj_theta(N, 0); // same value for all N
    double ref_traj_start = 1;
    double ref_traj_end = 20;
    double ref_traj_step = (ref_traj_end - ref_traj_start) / (N - 1);
    for (int i = 0; i < N; ++i) {
        ref_traj_x[i] = ref_traj_start + i * ref_traj_step;
    }

    // CasADi Opti environment
    Opti opti;

    // Define decision variables for control inputs and slack variables
    MX v = opti.variable(N);
    MX omega = opti.variable(N);
    MX delta_e = opti.variable(N);
    MX slack_cbf_1 = opti.variable(N);
    MX slack_cbf_2 = opti.variable(N);
    MX slack_cbf_3 = opti.variable(N);

    // Pack all controls and slack variables into a single decision variable
    MX u = vertcat(v, omega, delta_e, slack_cbf_1, slack_cbf_2, slack_cbf_3);

    // Initial guess
    // vector<double> init_u(N * 6);
    // fill(begin(init_u), end(init_u), 0.0);
    // fill(begin(init_u), begin(init_u) + N, 10.0); // v init
    // fill(begin(init_u) + N, begin(init_u) + 2 * N, 0.3); // omega init
    DM v_init = DM::ones(N, 1) * 10;    // 10*ones(N, 1)
    DM omega_init = DM::ones(N, 1) * 0.3;// ones(N,1)*0.3
    DM delta_e_init = DM::ones(N, 1) * 10;// ones(N, 1)*10
    DM slack_cbf_1_init = DM::zeros(N, 1);// ones(N, 1)*0
    DM slack_cbf_2_init = DM::ones(N, 1); // ones(N, 1)*1
    DM slack_cbf_3_init = DM::zeros(N, 1);// ones(N, 1)*0
    //cout<<"here"<<endl;
    // Concatenate all initial guesses into a single vector
    DM u0 = vertcat(v_init, omega_init, delta_e_init, slack_cbf_1_init, slack_cbf_2_init, slack_cbf_3_init);

    // Objective function
    //Function objective_fun = Function("objective_fun", {u}, { /* Define the objective here */ });
    MX cost = objective_function(v, omega, delta_e, slack_cbf_1, slack_cbf_2, slack_cbf_3, N);
    // Nonlinear constraints function
    //Function nl_constraints = Function("nl_constraints", {u}, { /* Define the constraints here */ });
    nl_constraints(opti, v, omega, delta_e, slack_cbf_1, slack_cbf_2, slack_cbf_3, init_pose_e, ref_traj_x, ref_traj_y, ref_traj_theta, pedestrian_position, N, T, enforced_rules, dist_from_ped);
    // Call the solver
    opti.minimize(cost);

    // Set initial condition
    opti.set_initial(u, u0);

    // Bounds and other constraints
    opti.subject_to(0 <= v);
    opti.subject_to(v <= 20);
    opti.subject_to(-deg2rad(50) <= omega);
    opti.subject_to(omega <= deg2rad(50));
    // Add other constraints defined in the nl_constraints
    //opti.subject_to(nl_constraints(u));

    // Solver options
    opti.solver("ipopt");

    // Solve the problem
    try {
        OptiSol sol = opti.solve();

        // Get the optimal solution
        DM u_opt = sol.value(u);

        // TODO: Post-process the solution to get state trajectories and other desired outputs
        // TODO: Include the plotting part if necessary or handle it differently in C++
        // Note that plotting likely needs a separate library such as matplotlib-cpp or an equivalent

        // Prints the optimized values (might need formatting according to desire)
        cout << "Optimal control trajectory: " << u_opt << endl;

    } catch (std::exception& e) {
        cerr << "Exception during optimization: " << e.what() << endl;
    }

    return 0;
}