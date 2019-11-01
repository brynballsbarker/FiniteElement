//---------------------------------------------------------------------------//
/*!
 * \file Main.cc
 */
//---------------------------------------------------------------------------//

#include "Solver.hh"

#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <fstream>

int main()
{

    // Create csv file for storing errors.
    std::string filename = "results/errors.csv";
    std::ofstream file( filename );
    file << "n,forder,error" << std::endl;

    // Initialize n, f_order, and the constant c
    std::array<int,4> n_vals = {10,100,1000,10000};
    int f_order;
    double c;

    // F is order 0
    f_order = 0;
    c = 3;

    std::cout << "f(x) = c\n";

    for ( int n : n_vals )
    {
        std::cout << "\tNumber of elements: " << n << std::endl;
        FEA::Solver solver( n, f_order, c );
        solver.initialize();
        solver.solve();
    }
    
    // F is order 1
    f_order = 1;
    c = 1;

    std::cout << "f(x) = x\n";

    for ( int n : n_vals )
    {
        std::cout << "\tNumber of elements: " << n << std::endl;
        FEA::Solver solver( n, f_order, c );
        solver.initialize();
        solver.solve();
    }

    // F is order 2
    f_order = 2;
    c = 1;

    std::cout << "f(x) = x^2\n";

    for ( int n : n_vals )
    {
        std::cout << "\tNumber of elements: " << n << std::endl;
        FEA::Solver solver( n, f_order, c );
        solver.initialize();
        solver.solve();
    }

    return 0;
}
