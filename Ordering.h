#pragma once

#include <Eigen/Sparse>
#include "SparseMatrix.h"

namespace CholUp {

template<class T>
void
createMatrixPermuted(CholUp::SparseMatrix<T>& A,
                     std::vector<Eigen::Triplet<T>>& triplets)
{
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> p;

    const int N =  std::max(std::max_element(triplets.begin(), triplets.end(),
                            [](const Eigen::Triplet<T>& t0, const Eigen::Triplet<T>& t1){return t0.row() < t1.row();})->row(),

                            std::max_element(triplets.begin(), triplets.end(),
                            [](const Eigen::Triplet<T>& t0, const Eigen::Triplet<T>& t1){return t0.col() < t1.col();})->col()) + 1;

    Eigen::SparseMatrix<double> B;
    B.resize(N, N);
    B.setFromTriplets(triplets.begin(), triplets.end());
    assert(B.isCompressed());
  //  assert((B - B.transpose().eval()).squared_norm() < 1e-8 && "matrix should be symmetric!");

    Eigen::internal::minimum_degree_ordering(B, p);

    if(A.perm) delete[] A.perm;
    A.perm = new int[N];

    for(int i = 0; i < N; ++i)
    {
        A.perm[p.indices()[i]] = i;
    }

    // reorder triplets
    for(auto& t : triplets)
    {
        t = Eigen::Triplet<T>(A.perm[t.row()], A.perm[t.col()], t.value());
    }

    B.setFromTriplets(triplets.begin(), triplets.end());

    A = CholUp::fromEigen(B);
}

} /* namespace CholUp */