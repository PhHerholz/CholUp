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

vector<Eigen::Triplet<double>>
loadTriplets(std::string fname)
{
    vector<Eigen::Triplet<double>> triplets;
    
    ifstream file(fname);
    int n;
    file >> n;
    
    for(int k = 0; k < n; ++k)
    {
        int i, j;
        double v;
        
        file >> i;
        file >> j;
        file >> v;
    
        triplets.emplace_back(i, j, v);
    }
    
    file.close();
    
    return triplets;
}

int
tripletDimensions(std::vector<Eigen::Triplet<double>>& triplets)
{
    const int N =  1 + std::max(std::max_element(triplets.begin(), triplets.end(),
                                             [](const Eigen::Triplet<double>& t0, const Eigen::Triplet<double>& t1){return t0.row() < t1.row();})->row(),

                                std::max_element(triplets.begin(), triplets.end(),
                                             [](const Eigen::Triplet<double>& t0, const Eigen::Triplet<double>& t1){return t0.col() < t1.col();})->col());

    return N;
}



int main(int argc, const char * argv[])
{
    // load data for update
    CholUp::SparseMatrix<double> A;
    auto triplets = loadTriplets("./data/tripletsLTL");
    auto roiIds = loadIds("./data/ids");


    // create matrix from triplets using amd reordering
    createMatrixPermuted(A, triplets);

    {
        // test multiply
        CholUp::Matrix<double> rhs(A.ncols, 3);
        rhs.fill();

        for(int i = 0; i < A.ncols * 3; ++i)
            rhs.data[i] = 0.01 * (rand() % 100);

        auto x = A * rhs;


        x.write("./data/y1");
        rhs.write("./data/y0");


        A.writeMatrixMarket("./data/AA.mtx");
    }



    // factorize & update
    Timer t0("Factor");
    CholUp::SupernodalCholesky<CholUp::SparseMatrix<double>> chol(A);
    t0.printTime("full");
    t0.reset();
    
    auto cholPart0 = chol.dirichletPartialFactor(A, roiIds);
    t0.printTime("partial");

    // solve linear system involving cholPart0
    // setup rhs
    CholUp::Matrix<double> rhs(A.ncols, 3);
    rhs.fill();
    rhs(roiIds[0], 0) = rhs(roiIds[0], 1) = rhs(roiIds[0], 2) = 1;
    auto rhs0 = rhs;

    // solve. Permutation is automatically taken care of. Rhs defines values for boundary conditions, values in rhs[roiIds] are replaced by the solve.
    cholPart0.solve(rhs);

    // write out solution
    rhs0.write("./data/b");
    rhs.write("./data/x");

    // refactor
    Eigen::SparseMatrix<double> eigenM2;
    Eigen::loadMarket(eigenM2, "./data/subLTL.mtx");
    CholUp::SparseMatrix<double> A2(eigenM2);

    Timer t3("Full refactor");
    CholUp::SupernodalCholesky<CholUp::SparseMatrix<double>> chol2(A2);
    t3.printTime();

    return 0;
}
