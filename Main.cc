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
    if ( beam )
        std::cout << "\nSolving the beam problem with p = " << p << "\n\n";
    else
        std::cout << "\nSolving u_xx = f with p = " << p << "\n\n";

    // Define beam constants [M,Q,f,E,I].
    std::array<double,5> consts = {0.0,0.0,10.*pow(.005,3),pow(10.,6),pow(.005,4)/12.};

    // Create csv file for storing errors.
    std::string filename = "results/errors.csv";
    std::ofstream file( filename );
    file << "n,forder,hval,error" << std::endl;

    // Initialize n, f_order, and the constant c
    std::array<int,3> n_vals = {1,10,100};
    std::array<double,5> h_vals = {0.1, 0.01, 0.005, 0.002, 0.001};
    std::array<std::string, 3> f_labs = {"f(x) = c\n","f(x) = x\n","f(x) = x^2\n"};
    int f_order = 0;
    double c = 3;

    if ( beam )
    {
        // Solving the beam

        for ( double h : h_vals )
        {
            std::cout << "    h = " << h << "\n";
            consts[2] = 10.*pow(h,3);
            consts[4] = .005*pow(h,3)/12.;

            for ( int n : n_vals )
            {
                std::cout << "\tNumber of elements: " << n << std::endl;
                FEA::Solver solver( n, f_order, c, p, beam, consts );
                solver.initialize();
                solver.solve();
            }
        }

    }


    else
    {
        for ( int f = 0; f < 3; ++f )
        {
            std::cout << f_labs[f];
            for ( int n : n_vals )
            {
                std::cout << "\tNumber of elements: " << n << std::endl;
                FEA::Solver solver( n, f, c, p, beam, consts );
                solver.initialize();
                solver.solve();
            }
        }

    }

    return 0;
}
