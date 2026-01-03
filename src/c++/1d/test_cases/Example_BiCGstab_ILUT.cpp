#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>

int main() {
    using namespace Eigen;

    // Define a sparse matrix A (e.g., 5x5 diagonal + small off-diagonals)
    typedef SparseMatrix<double> SpMat;
    SpMat A(5, 5);
    A.insert(0, 0) = 4; A.insert(0, 1) = -1;
    A.insert(1, 0) = -1; A.insert(1, 1) = 4; A.insert(1, 2) = -1;
    A.insert(2, 1) = -1; A.insert(2, 2) = 4; A.insert(2, 3) = -1;
    A.insert(3, 2) = -1; A.insert(3, 3) = 4; A.insert(3, 4) = -1;
    A.insert(4, 3) = -1; A.insert(4, 4) = 4;

    A.makeCompressed();

    // Preconditioner type: ILUT
    typedef IncompleteLUT<double> ILUT;

    // Create and compute the preconditioner ONCE
    ILUT ilu;
    ilu.setDroptol(1e-3);      // Drop tolerance
    ilu.setFillfactor(10);     // Fill factor
    ilu.compute(A);            // Factorize

    if (ilu.info() != Success) {
        std::cerr << "Preconditioner computation failed!" << std::endl;
        return -1;
    }

    // Create the solver and assign A and the preconditioner
    BiCGSTAB<SpMat, ILUT> solver;
    solver.setPreconditioner(ilu); // reuse this ilu for all solves
    solver.compute(A);             // Analyze A structure (no refactorization of ilu)

    if (solver.info() != Success) {
        std::cerr << "Solver initialization failed!" << std::endl;
        return -1;
    }

    // Solve for multiple RHS vectors
    for (int i = 0; i < 3; ++i) {
        VectorXd b = VectorXd::Random(5);
        VectorXd x = solver.solve(b);

        std::cout << "Solve " << i << " | Iterations: " << solver.iterations()
                  << " | Error: " << solver.error() << std::endl;
    }

    return 0;
}
