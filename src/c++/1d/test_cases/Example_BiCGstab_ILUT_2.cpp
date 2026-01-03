#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>

int main() {
    using namespace Eigen;

    // Define a sparse matrix
    typedef SparseMatrix<double> SpMat;
    typedef Triplet<double> T;
    std::vector<T> tripletList;
    const int size = 5;

    // Fill matrix with simple values
    for (int i = 0; i < size; ++i) {
        tripletList.emplace_back(i, i, 4.0);
        if (i > 0)
            tripletList.emplace_back(i, i - 1, -1.0);
        if (i < size - 1)
            tripletList.emplace_back(i, i + 1, -1.0);
    }

    SpMat A(size, size);
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    // Right-hand side vector
    VectorXd b(size);
    b << 1, 2, 3, 4, 5;

    // Initial guess
    VectorXd x0 = VectorXd::Zero(size);

    // Create and initialize fixed preconditioner (IncompleteLUT)
    IncompleteLUT<double> ilu;
    ilu.setDroptol(1e-4);  // Optional tuning
    ilu.compute(A);        // Factorize the matrix once

    if (ilu.info() != Success) {
        std::cerr << "Preconditioner factorization failed." << std::endl;
        return -1;
    }

    // Create solver and assign fixed preconditioner
    BiCGSTAB<SpMat, IncompleteLUT<double>> solver;
    solver.setMaxIterations(1000);
    solver.setTolerance(1e-10);
    solver.preconditioner() = ilu;  // Use fixed preconditioner
    solver.compute(A);              // Only symbolic analysis; numerical factorization was done

    VectorXd x = solver.solveWithGuess(b, x0);

    std::cout << "Solution x:\n" << x << std::endl;
    std::cout << "Iterations: " << solver.iterations() << std::endl;
    std::cout << "Estimated error: " << solver.error() << std::endl;

    return 0;
}
