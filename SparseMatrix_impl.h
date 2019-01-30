#ifndef SparseMatrix_impl_h
#define SparseMatrix_impl_h

#include <fstream>
#include <complex>
#include <assert.h>
#include <tuple>
#include <iostream>
#include <iomanip>
#include <iostream>

#include "Timer.hpp"

template<class T>
bool SparseMatrix<T>::compareTriplet(const Triplet& t0, const Triplet& t1)
{
    using namespace std;
    
    if(get<1>(t0) < get<1>(t1)) return true;
    else if(get<1>(t0) > get<1>(t1)) return false;
    else if(get<0>(t0) < get<0>(t1)) return true;
    else return false;
}

template<class T>
SparseMatrix<T>::SparseMatrix()
{
}

template<class T>
SparseMatrix<T>::SparseMatrix(const int nrows_, const int ncols_)
: nrows(nrows_), ncols(ncols_)
{
    col = new int[ncols + 1];
    fill_n(col, ncols + 1, 0);
}

template<class T>
SparseMatrix<T>::SparseMatrix(SparseMatrix&& A)
{
    *this = A;
}

template<class T>
SparseMatrix<T>::SparseMatrix(const SparseMatrix& A)
{
    *this = A;
}

template<class T>
SparseMatrix<T>::SparseMatrix(Eigen::SparseMatrix<T, Eigen::ColMajor, int>& eigenMatrix)
:   nrows((int)eigenMatrix.rows()),
    ncols((int)eigenMatrix.cols()),
    nnz((int)eigenMatrix.nonZeros()),
    col(eigenMatrix.outerIndexPtr()),
    row(eigenMatrix.innerIndexPtr()),
    vals(eigenMatrix.valuePtr()),
    dataBorrowed(true)
{
    setDiagonalIndizes();
}

template<class T>
SparseMatrix<T>::~SparseMatrix()
{
    if(!dataBorrowed)
    {
        if(col) delete[] col;
        if(row) delete[] row;
        if(vals) delete[] vals;
    }
    
    if(diag) delete[] diag;
}


template<class T>
SparseMatrix<T>& SparseMatrix<T>::operator=(const SparseMatrix& A)
{
    using namespace std;
    
    nnz = A.nnz;
    ncols = A.ncols;
    nrows = A.nrows;
    
    row = new int[nnz];
    copy_n(A.row, nnz, row);
    
    col = new int[ncols + 1];
    copy_n(A.col, ncols + 1, col);
    
    vals = new T[nnz];
    copy_n(A.vals, nnz, vals);
    
    if(A.diag)
    {
        diag = new int[min(ncols, nrows)];
        copy_n(A.diag, min(ncols, nrows), diag);
    }
    
    return *this;
}


template<class T>
SparseMatrix<T>& SparseMatrix<T>::operator=(SparseMatrix&& A)
{
    dataBorrowed = A.dataBorrowed;
    
    nnz = A.nnz;
    ncols = A.ncols;
    nrows = A.nrows;
    
    
    
    row = A.row;
    col = A.col;
    vals = A.vals;
    nrows = A.nrows;
    diag = A.diag;
    
    A.row = A.col = A.diag = nullptr;
    A.vals = nullptr;
    
    return *this;
}


template<class T>
void
SparseMatrix<T>::setDiagonalIndizes()
{
    if(diag) delete [] diag;
    diag = new int[ncols];
    std::fill_n(diag, ncols, -1);
    
    for(int i = 0; i < ncols; ++i)
    {
        for(int j = col[i]; j < col[i+1]; ++j)
        {
            if(row[j] == i)
            {
                diag[i] = j;
                break;
            }
        }
    }
}

template<class T>
void SparseMatrix<T>::addToDiagonal(const T val)
{
    for(int i = 0; i < ncols; ++i)
    {
        for(int j = col[i]; j < col[i+1]; ++j)
        {
            if(row[j] == i)
            {
                vals[j] += val;
                break;
            }
        }
        
    }
}

template<class T>
void
SparseMatrix<T>::setTriplets(std::vector<Triplet>& triplets, int ncols)
{
    using namespace std;
    
    assert(all_of(triplets.begin(), triplets.end(), [=](const Triplet& t){return ncols > get<2>(t);}));
    sort(triplets.begin(), triplets.end(), compareTriplet);
    
    if(col) delete[] col;
    col = new int [ncols + 1];
    
    int nt = triplets.size();
    
    if(row) delete[] row;
    row = new int[nt];
    
    if(vals) delete[] vals;
    vals = new T[nt];
    
    int cnt = 0;
    
    while(cnt <= get<1>(triplets.front()))
        col[cnt++] = 0;
    
    int i = get<0>(triplets.front());
    int j = get<1>(triplets.front());
    T sum = T();
    
    cnt = 0;
    int colCnt = 0;
    
    for(auto& t : triplets)
    {
        if(get<0>(t) == i && get<1>(t) == j)
        {
            sum += get<2>(t);
        } else
        {
            row[cnt] = i;
            vals[cnt] = sum;
            
            ++cnt;
            
            while(get<1>(t) > j++)
                col[colCnt++] = cnt;
            
            tie(i, j, sum) = t;
        }
    }
    
    row[cnt] = i;
    vals[cnt] = sum;
    
    ++cnt;
    
    while(j >= ncols)
    {
        col[colCnt++] = cnt;
    }

    col[colCnt++] = cnt;
}

template<class T>
void
SparseMatrix<T>::writeMatrixMarket(const std::string& filename, const bool symmetric) const
{
    using namespace std;
    
    ofstream file(filename);
    file << "%%MatrixMarket matrix coordinate real ";
    
    if(symmetric) file << "symmetric";
    else file << "general";
    
    file << endl;
    file << setiosflags(ios::fixed);
    file.precision(20);
    
    if(!nnz || !ncols)
    {
        file.close();
        return;
    }
    
    const int nrows2 = nrows == -1 ? ncols : nrows;
    
    file << nrows2 << " " << ncols << " " << nnz << "\n";
    for(int i = 0; i < ncols; ++i)
        for(int j = col[i]; j < col[i+1]; ++j)
            file << row[j] + 1 << " " << i + 1 << " " << vals[j] << "\n";
    
    file.close();
}

#endif /* SparseMatrix_impl_h */
