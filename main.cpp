#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>

#include "SparseMatrix.h"
#include "SupernodalCholesky.h"

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "Timer.hpp"
#include "Ordering.h"
#include <Eigen/OrderingMethods>

using namespace std;

// just loads a bunch of binary integers
vector<int>
loadIds(std::string fname)
{
    ifstream file(fname, std::fstream::binary);

    if(!file.good())
    {
        std::cout << "file not found" << std::endl;
        return {};
    }

    std::vector<char> ids((std::istreambuf_iterator<char>(file)),
                          (std::istreambuf_iterator<char>()));

    int* idPtr = (int*)ids.data();
    int nids = (int)ids.size() / sizeof(int);

    return std::vector<int>(idPtr, idPtr + nids);
}




int main(int argc, const char * argv[])
{
    Eigen::SparseMatrix<double> A;
    Eigen::loadMarket(A, "../data/LTL.mtx");
    assert(A.isCompressed());

    auto roiIds = loadIds("../data/ids");

    // factorize & update
    CholUp::Timer t0("Factor");
    CholUp::SupernodalCholesky<CholUp::SparseMatrix<double>> chol(A);
    t0.printTime("full");
    t0.reset();
    
    auto cholPart0 = chol.dirichletPartialFactor(roiIds);
    t0.printTime("partial");

    // solve linear system involving cholPart0
    // setup rhs
    CholUp::Matrix<double> rhs(A.cols(), 1);
    rhs.fill();
   
    rhs(roiIds[0], 0) = 1;
    auto rhs0 = rhs;

    // solve. Permutation is automatically taken care of. Rhs defines values for boundary conditions, values in rhs[roiIds] are replaced by the solve.
    cholPart0.solveL_RowMajor<1>(rhs.data);
    cholPart0.solveLT_RowMajor<1>(rhs.data);
    


    return 0;
}
