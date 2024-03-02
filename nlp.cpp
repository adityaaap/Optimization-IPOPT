#include <iostream>
#include <fstream>
#include <ctime>
#include <casadi/casadi.hpp>

using namespace casadi;

int main(){
    // Parameters similar to the MATLAB problem
    std::vector<double> E_dem = {100, 50, 70, 90, 30, 150};
    std::vector<double> E_wind = {80, 130, 30, 90, 20, 20};

    // Create an optimization environment
    Opti opti;
    MX x = opti.variable(13);

    // Objective: sum of the first 6 elements of x
    //opti.minimize(sum1(x(Slice(0, 6))));
    opti.minimize(dot(x(Slice(0, 6)), x(Slice(0, 6))));
    // Constraints
    opti.subject_to(x(12) * x(12) == 10);
    opti.subject_to(x(6) == 0);

    for (int i = 0; i < 6; ++i) {
        // Corresponds to the inequality constraint expression in MATLAB
        opti.subject_to(-x(i + 6) + x(i + 7) - E_wind[i] + E_dem[i] - x(i) <= 0);
    }

    // Upper and Lower bounds for x
    for (int i = 0; i < 6; ++i) {
        opti.subject_to(x(i) >= 0); 
        opti.subject_to(x(i) <= 50); 
    }
    for (int i = 6; i < 13; ++i) {
        opti.subject_to(x(i) >= 0);
        opti.subject_to(x(i) <= 100); 
    }

    // Define the solver options and solve the problem
    opti.solver("ipopt"); // Using the IPOPT solver

    try {
        OptiSol sol = opti.solve();
        DM x_opt = sol.value(x);

        std::cout << "The optimal value of x is: " << x_opt << std::endl;

    } catch (std::exception &e) {
        std::cerr << "Exception during optimization: " << e.what() << std::endl;
    }

    return 0;
}