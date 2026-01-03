SparseMatrix<double> A(N, N);
// Only set structure once
A.reserve(VectorXi::Constant(N, 5)); // 5 nonzeros per row
// ... triplets setup ...
A.setFromTriplets(triplets.begin(), triplets.end()); // only once

// In time loop:
for (...) {
    // Update A’s values
    for (int k = 0; k < A.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
            it.valueRef() = compute_new_value(it.row(), it.col(), t);  // <-- update only
}
