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

int main(int argc, char** argv)
{
    // Get command line arguements: p val and beam bool.
    int p = ( argc > 1 ) ? std::atoi( argv[1] ) : 2;
    bool beam = ( argc ==3 && ( std::string( argv[2]) == "yes" ) ) ? true : false;

    // Send back command line arguments.
    std::cout << "your p value is " << p << std::endl;
    std::cout << "your beam values is " << beam << std::endl;

    // Define beam constants [M,Q,f,E,I].
    std::array<double,5> consts = {10.0,10.0,10.*.005*.005*.005,1000000.,10.};

    // Create csv file for storing errors.
    std::string filename = "results/errors.csv";
    std::ofstream file( filename );
    file << "n,forder,error" << std::endl;

    // Initialize n, f_order, and the constant c
    std::array<int,4> n_vals = {10,100,1000,10000};
    int f_order;
    double c;

    if ( beam )
    {
        // Solving the beam
        f_order = 0;
        c = 0;

        std::cout << "Beam\n";

        for ( int n : n_vals )
        {
            std::cout << "\tNumber of elements: " << n << std::endl;
            FEA::Solver solver( n, f_order, c, p, beam, consts );
            solver.initialize();
            solver.solve();
        }

    }


    else
    {

        // F is order 0
        f_order = 0;
        c = 3;

        std::cout << "f(x) = c\n";

        for ( int n : n_vals )
        {
            std::cout << "\tNumber of elements: " << n << std::endl;
            FEA::Solver solver( n, f_order, c, p, beam, consts );
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
            FEA::Solver solver( n, f_order, c, p, beam, consts );
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
            FEA::Solver solver( n, f_order, c, p, beam, consts );
            solver.initialize();
            solver.solve();
        }

    }

    return 0;
}
