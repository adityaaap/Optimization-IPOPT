#include <casadi/casadi.hpp>
#include <iostream>

using namespace casadi;

int main() {
    // System dynamics
    // DM A = {{.9988, .0494},
    //         {-.0494, .9741}};

    DM A = DM::zeros(2,2);
    A(0,0) = 0.99;
    A(0,1) = 0.49;
    A(1,0) = -0.05;
    A(1,1) = 0.97;
    DM B = {{0.0012}, {0.05}};

    // Initial condition and horizon length
    DM x0 = {{10}, {0}};
    int N = 300; // Horizon length

    // LQR weighting matrices
    DM Q_bar = {{1, 0},
                {0, 0.1}};
    double r_bar = 10;

    // Bounds
    double x1Max = 10, x1Min = -10;
    double x2Max = 10, x2Min = -10;
    double uMax = 5, uMin = -5;

    // Combined state and control vector size
    int n = Q_bar.size1(), m = 1; // A is n-by-n, and B is n-by-m

    // Form the new Q matrix
    DM Q = DM::zeros(N * (n + m), N * (n + m));
    // for (int i = 0; i < N; ++i) {
    //     Q(Slice(n * i, n * (i + 1)), Slice(n * i, n * (i + 1))) = Q_bar;
    //     Q(Slice(N * n + m * i, N * n + m * (i + 1)),
    //       Slice(N * n + m * i, N * n + m * (i + 1))) = r_bar;
    // }

    for(int i=0; i<N; i++){
        for(int j =0; j<N ; j++){
            if(i == j){
                Q(Slice(N * n + m * (i) + 1, N * n + (i+1) * m), Slice(N * n + m * (i) + 1, N * n + (i+1) * m)) = r_bar;

                //Q(Slice(N * n + m * (i-1) + 1, N * n + i * m), Slice(N * n + m * (i-1) + 1, N * n + i * m)) = r_bar;
            }
        }
    }
    for(int i=0; i<N; i++){
        for(int j =0; j<N ; j++){
            if(i == j){
                //Q(Slice(n * (i-1) + 1, i * n), Slice(n * (i-1) + 1, i * n)) = Q_bar;
                Q(Slice(n * (i) + 1, (i+1) * n), Slice(n * (i) + 1, (i+1) * n)) = Q_bar;

            }
        }
    }


    // Form the inequality constraints (bounds on state and input)
    //DM A1 = DM::zeros(6 * N, (n + m) * N);
    //DM b1 = DM::repmat(DM::vertcat({x1Max, -x1Min, x2Max, -x2Min, uMax, -uMin}), N, 1);
    // Populate A1 with identity matrices for inequality constraints
    // for (int i = 0; i < N; ++i) {
    //     A1(i, i) = 1; // x1 max constraint
    //     A1(N + i, i) = -1; // x1 min constraint
    //     A1(2 * N + i, n * N + i) = 1; // x2 max constraint
    //     A1(3 * N + i, n * N + i) = -1; // x2 min constraint
    //     A1(4 * N + i, (n + m) * N - N + i) = 1; // u max constraint
    //     A1(5 * N + i, (n + m) * N - N + i) = -1; // u min constraint
    // }

    DM eye_N = DM::eye(N);
    DM neg_eye_N = -1 * eye_N;

    // Create zero matrices with appropriate dimensions
    DM zeros_N = DM::zeros(N, N);

    // Construct the blocks
    DM block1 = DM::horzcat({eye_N, zeros_N, zeros_N});
    DM block2 = DM::horzcat({neg_eye_N, zeros_N, zeros_N});
    DM block3 = DM::horzcat({zeros_N, eye_N, zeros_N});
    DM block4 = DM::horzcat({zeros_N, neg_eye_N, zeros_N});
    DM block5 = DM::horzcat({zeros_N, zeros_N, eye_N});
    DM block6 = DM::horzcat({zeros_N, zeros_N, neg_eye_N});

    // Concatenate all blocks to form A1
    DM A1 = DM::vertcat({block1, block2, block3, block4, block5, block6});

    DM ones_N = DM::ones(N,1);
    DM row1 = x1Max * ones_N;
    DM row2 = -x1Min * ones_N;
    DM row3 = x2Max * ones_N;
    DM row4 = -x2Min * ones_N;
    DM row5 = uMax * ones_N;
    DM row6 = -uMin * ones_N;

    DM b1 = DM::vertcat({row1, row2, row3, row4, row5, row6});

    // Formulate the equality constraint (dynamics)
    DM A2 = DM::zeros(n * N, (n + m) * N);
    // DM B2 = DM::zeros(n * N, 1);
    // B2(Slice(0, n)) = A * x0;

    DM row1_b2 = mtimes(A, x0);
    DM row2_b2 = DM::zeros(n * (N-1), 1);
    DM B2 = DM::vertcat({row1_b2, row2_b2});

    for(int i =0; i< N; i++){
        if(i ==1){
            //A2(Slice(n * (i-1) + 1, n * i), Slice(n * (i-1) + 1, n * i)) = eye_N;
            A2(Slice(n * (i) + 1, n * (i+1)), Slice(n * (i) + 1, n * (i+1))) = eye_N;

            //A2(Slice(n*(i-1)+1,n*i), Slice(N*n+m*(i-1)+1, N*n+m*i)) = -B;
            A2(Slice(n*(i)+1,n*(i+1)), Slice(N*n+m*(i)+1, N*n+m*(i+1))) = -B;

        }
        else{
            // A2(Slice(n*(i-1)+1, n*i), Slice(n*(i-1)+1, n*i)) = eye_N;
            // A2(Slice(n*(i-1)+1, n*i), Slice(N*n+m*(i-1)+1, N*n+m*i)) = -B;
            // A2(Slice(n*(i-1)+1, n*i), Slice(n*(i-2)+1, n*(i-1))) = -A;

            A2(Slice(n*(i)+1, n*(i+1)), Slice(n*(i)+1, n*(i+1))) = eye_N;
            A2(Slice(n*(i)+1, n*(i+1)), Slice(N*n+m*(i)+1, N*n+m*(i+1))) = -B;
            A2(Slice(n*(i)+1, n*(i+1)), Slice(n*(i-1)+1, n*(i))) = -A;
        }
    }


    // for (int i = 1; i < N; ++i) {
    //     // Use DM::eye instead of identity
    //     A2(Slice(n * i, n * (i + 1)), Slice(n * i, n * (i + 1))) = DM::eye(n);
    //     A2(Slice(n * i, n * (i + 1)), Slice((n + m) * i - m, (n + m) * i)) = -B;
    //     if (i != 1)
    //         A2(Slice(n * i, n * (i + 1)), Slice(n * (i - 1), n * i)) = -A;
    // }

    // Set up and solve the QP
    Opti opti;
    MX X = opti.variable((n + m) * N, 1);
    opti.minimize(0.5 * mtimes(X.T(), Q * X));
    opti.subject_to(mtimes(A1, X) <= b1);
    opti.subject_to(mtimes(A2, X) == B2);

    // Solving the QP
    opti.solver("osqp", {{"verbose", true}}); // Using OSQP solver
    OptiSol sol = opti.solve();

    // Access the optimal value of the decision variable
    DM u_opt = sol.value(X);

    // Print result
    std::cout << "Optimal control sequence: " << u_opt << std::endl;
    return 0;
}