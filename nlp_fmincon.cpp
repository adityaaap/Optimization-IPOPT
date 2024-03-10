/***************************************************************************/
// - This file contains C++ equivalent implementation of a matlab code which uses fmincon.
// - This file can be used as an example reference to formulate optimization problems in C++ 
//   as one would have done in Matlab using fmincon.
// - The matlab code is give at the end in commented format
/***************************************************************************/

#include <iostream>
#include <casadi/casadi.hpp>

using namespace casadi;
using namespace std;

void nonlcon(Opti &opti, MX &x){
    opti.subject_to(2*x(0)*x(0) + x(1) == 4);
}

int main() {
    // Optimization variable
    Opti opti;
    MX x = opti.variable(2);

    // Define the objective function
    MX f = pow(x(0) - 2, 2) + pow(x(1) - 3, 2);
    opti.minimize(f);
    cout<<"Obj set"<<endl;

    // Linear inequality constraints: A * x <= b
    // Matrix<double> A = {{1, 1}, {-1, -1}};
    // std::vector<double> b = {2, -1};
    // opti.subject_to(mtimes(A, x) <= b);
    // cout<<"A * x <= b"<<endl;
    DM A = DM(2, 2);
    A(0, 0) = 1; A(0, 1) = 1;
    A(1, 0) = -1; A(1, 1) = -1;
    DM b = DM(2, 1);
    b(0, 0) = 2;
    b(1, 0) = -1;
    opti.subject_to(mtimes(A, x) - b <= vector<double>({0.0, 0.0}));
    cout<<"A * x <= b"<<endl;


    // Linear equality constraint: Aeq * x == beq
    // Matrix<double> Aeq = {{2, 2}};
    // double beq = 4;
    // opti.subject_to(mtimes(Aeq, x) == beq);
    DM Aeq = DM(1, 2);  // Create a 1x2 dense matrix
    Aeq(0, 0) = 2; Aeq(0, 1) = 2;
    double beq = 4;
    opti.subject_to(mtimes(Aeq, x) == beq);
    cout<<"Aeq * x == beq"<<endl;

    // Lower and upper bounds on x
    opti.subject_to(opti.bounded(0, x, 5)); // x vector bounded between 0 and 5

    // Initial guess for the decision variables
    std::vector<double> x0 = {1, 2};
    opti.set_initial(x, x0);

    // Nonlinear constraint
    // opti.subject_to(opti.callback([=](const MX& x){
    //     MX ceq = 2 * pow(x(0), 2) + x(1) - 4;
    //     return std::vector<MX>{ceq}; // No inequality constraints, just the equality constraint
    // }, "nlcon", x));
    cout<<"before nonlinear"<<endl;
    nonlcon(opti, x);


    // Solver options
    opti.solver("ipopt", Dict{{"verbose", true}, {"ipopt.print_level", 5}});

    try {
        // Solve the problem
        OptiSol sol = opti.solve();

        // Extract the optimized values
        DM x_opt = sol.value(x);
        DM f_opt = sol.value(f);

        std::cout << "Optimal x: " << x_opt << std::endl;
        std::cout << "Optimal objective: " << f_opt << std::endl;

    } catch (std::exception &e) {
        std::cerr << "Exception during optimization: " << e.what() << std::endl;
    }

    return 0;
}

/*********************************************************************************/
// Corresponding Matlab Code using fmincon

// %% Define the objective function
// fun = @(x) (x(1) - 2)^2 + (x(2) - 3)^2; % define f(x)

// %% Define the Linear inequality constraints: A*x <= b
// A = [1, 1; -1, -1]; % x1 + x2 <= 2 and -x1 -x2 <= -1
// b = [2; -1];

// %% Define the Linear equality constraint: Aeq*x = beq
// Aeq = [2, 2]; % 2x1 + 2x2 = 4
// beq = 4;

// %% Define lower and upper bounds for the decision variables
// lb = [0; 0]; % 0 <= x1 <= 5
// ub = [5; 5]; % 0 <= x2 <= 5

// %% Initial guess
// x0 = [1,2];

// %% Customize optimization process
// options = optimoptions('fmincon', 'Display', 'iter','Algorithm','interior-point','MaxIterations',500, 'MaxFunctionEvaluations',3000,'PlotFcn', {@optimplotfval, @optimplotx});

// %% Call fmincon to minimize the function with constraints
// [x_opt, f_opt] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, @nonlcon, options);

// %%
// function [c, ceq] = nonlcon(x)
//     % Nonlinear equality constraint: 2*(x1)^2 + x2 = 4
//     ceq = 2 * x(1)^2 + x(2) - 4;
    
//     % There are no inequality constraints, so c is empty
//     c = [];
// end