#pragma once

#include <Eigen/Sparse>
#include "SparseMatrix.h"

namespace CholUp {

template<class T>
CholUp::SparseMatrix<T>
permuteMatrix(const Eigen::SparseMatrix<T>& A, std::vector<int>& perm)
{
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> p;
    Eigen::SparseMatrix<T> B = A;
    
    Eigen::internal::minimum_degree_ordering(B, p);
    
    // TODO: optimize
    B = p.transpose() * A * p;
    auto ret = CholUp::fromEigen(B);
    
    perm.resize(A.rows());
    
    for(int i = 0; i < A.rows(); ++i)
    {
        perm[p.indices()[i]] = i;
    }
    
    return std::move(ret);
}

} /* namespace CholUp */
