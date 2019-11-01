//---------------------------------------------------------------------------//
/*!
 * \file Solver.cc
 */
//---------------------------------------------------------------------------//

#include "Solver.hh"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <fstream>

#include <chrono>

namespace FEA
{
//---------------------------------------------------------------------------//
// Constructor.
Solver::Solver( const int n,
                const int f_order,
                const double c )
    : d_n( n )
    , d_f_order( f_order )
{
    // Define f.
    d_f = std::make_shared<Force>( f_order, c );

    // Define the solution.
    d_u = std::make_shared<TrueSol>( f_order, c );

    // Create the mesh.
    d_mesh = std::make_shared<Mesh>( n );
}

//---------------------------------------------------------------------------//
// Initialize the problem by initializing the elements.
void Solver::initialize()
{
    // Create a vector of elements.
    int num_elements = d_mesh->totalNumElements();

    // Loop through the number of elements.
    for ( int c = 0; c < num_elements; ++c )
    {
        // Create new element with given id.
        FEA::Element element;
        d_mesh->initializeElement( c, element );

        // Add element to storage vector. 
        d_elements.push_back( element );
    }
}

//---------------------------------------------------------------------------//
// Solve for the FEA Approximate Solution.
void Solver::solve() 
{
    // Allocate mesh fields.
    int num_nodes = d_mesh->totalNumNodes();

    // Initialize K.
    std::vector<std::vector<double> > K(
            num_nodes,std::vector<double>(num_nodes, 0.0));
    K[num_nodes-1][num_nodes-1] = 1.0;

    // Initialize F.
    std::vector<double> F(num_nodes, 0.0);

    // Loop over elements.
    for ( auto& e : d_elements )
    {
        // Add local stiffness matrix to global
        updateK( K, e );

        // Add local force vector to global
        updateF( F, e );

    }

    // Solve Kd = F using Gauss Elimination.
    std::vector<double> d( F.begin(), F.end() );;
    gauss( d, K );

    // Initialize domain and solution vectors.
    std::vector<double> domain(2*num_nodes - 1,0.0);
    std::vector<double> uh(2*num_nodes - 1, 0.0);
    std::vector<double> u(2*num_nodes - 1, 0.0);

    // Evaluate the solution across domain.
    createDomain( domain );
    evaluateApprox( domain, d, uh );
    evaluateTrue( domain, u );

    // Compute the error.
    double error, nodes, midpoints;
    computeError( error, d );
    checkBehavior( nodes, midpoints, u, uh );

    std::cout << "\t\tError on nodes: \t\t" << nodes << std::endl;
    std::cout << "\t\tError on midpoints: \t" << midpoints << std::endl;

    // Store the solutions for plotting.
    plot( domain, uh, u );
    std::cout << "\t\tL2 Error: \t\t\t\t" << error << std::endl;

}

//---------------------------------------------------------------------------//
// Update the global stiffness matrix.
void Solver::updateK( std::vector<std::vector<double> >& K, 
                      FEA::Element& e) 
{
        // Get global element position
        int i = e.ind1;
        int j = e.ind2;

        // Update stiffness matrix.
        K[i][i] += e.k[0][0];
        K[i][j] += e.k[0][1];
        K[j][i] += e.k[1][0];
        K[j][j] += e.k[1][1];

}

//---------------------------------------------------------------------------//
// Update the global force vector.
void Solver::updateF( std::vector<double>& F,
                      FEA::Element& e) 
{
        // Get global element position
        int i = e.ind1;
        int j = e.ind2;

        // Calculate integral.
        std::array<double,2> f;
        integrateRHS( e, f);

        // Update force vector. 
        F[i] += f[0];
        F[j] += f[1];

}

//---------------------------------------------------------------------------//
// Compute the inner products (N_1,f) and (N_2,f).
void Solver::integrateRHS( const FEA::Element& e,
                           std::array<double,2>& f) 
{
    // Get global element endpoints.
    double left = e.n1;
    double right = e.n2;

    // Compute the force value at the endpoints.
    double f_1;
    double f_2;
    d_f->forceValue( left, f_1 );
    d_f->forceValue( right, f_2 );

    // Approximate the inner products. 
    double h = right - left;
    f[0] = h * (2.*f_1 + f_2)/6.;
    f[1] = h * (f_1 + 2.*f_2)/6.;

    // Check if this is the last element.
    if ( right == 1. )
        f[1] = 0.0;
}

//---------------------------------------------------------------------------//
// Use Gauss Elimination to solve the system Kd = F.
void Solver::gauss( std::vector<double>& d,
                    std::vector<std::vector<double> > K ) 
{
    double scale;
    int n = d_mesh->totalNumNodes();
    
    // Forward reduction.
    for ( int i = 0; i < n-2; ++i )
    {
        // Scale row so first entry is 1.
        scale = K[i][i];
        for ( int j = i; j < i+2; ++j )
            K[i][j] /= scale;
        d[i] /= scale;
        
        // Eliminate other leading entries.
        scale = K[i+1][i];
        for ( int j = i; j < i+3; ++j )
            K[i+1][j] -= scale * K[i][j];
        d[i+1] -= scale * d[i];
    }

    // Handle the second to last row.
    scale = K[n-2][n-2];
    K[n-2][n-2] /= scale;
    d[n-2] /= scale;

    // Backward reduction. 
    for ( int i = n-3; i >= 0; --i )
    {
        scale = K[i][i+1];
        K[i][i+1] -= scale * K[i+1][i+1];
        d[i] -= scale * d[i+1];
    }
}

//---------------------------------------------------------------------------//
// Create the domain which consists of each node and each element midpoint.
void Solver::createDomain( std::vector<double>& domain )
{
    int dom_len = domain.size();
    double h_step = 1. / ( dom_len - 1. );

    for ( int i = 0; i < dom_len; ++i )
        domain[i] = i * h_step;
}

//---------------------------------------------------------------------------//
// Evaluate the FEA solution over the domain. 
void Solver::evaluateApprox( const std::vector<double>& domain,
                             const std::vector<double>& d,
                             std::vector<double>& uh ) 
{
    int dom_len = domain.size();

    double x;
    int element_id;
    double ref_coord;
    std::array<double,2> values;

    // Loop through each value in domain
    for ( int i = 0; i < dom_len; ++i )
    {
        x = domain[i];

        // Determine which element x is in.
        d_mesh->locateX( x, element_id );

        // Find local position of x in the given element.
        d_mesh->mapGlobalToLocalFrame( x, element_id, ref_coord );

        // Get the value of the shape functions at x.
        d_mesh->shapeFunctionValues( ref_coord, values );

        // Compute the value of uh at x.
        auto& e = d_elements[ element_id ];
        uh[i] = d[e.ind1] * values[0] + d[e.ind2] * values[1];
    }
}

//---------------------------------------------------------------------------//
// Evalute the true solution over the domain.
void Solver::evaluateTrue( const std::vector<double>& domain, 
                           std::vector<double>& u )
{
    int dom_len = domain.size();

    for ( int i = 0; i < dom_len; ++i )
    {
        u[i] = d_u->trueValue( domain[i] );
    }
}

//---------------------------------------------------------------------------//
// Compute the L2 norm of |u - uh|.
void Solver::computeError( double& error,
                    const std::vector<double>& d )
{
    double u_ksi;
    double uh_ksi;

    double x;
    std::array<double,2> values;

    // Define gaussian quadrature interpolation points. 
    std::array<double,3> ksi = {-pow(3./5.,.5), 0.0, pow(3./5.,.5)};
    std::array<double,3> w = {5./9.,8./9.,5./9.};

    double h = d_mesh->getElementWidth();
    double dx = h / 2.;
    double sum = 0.0;

    // Loop over elements.
    for ( auto& e : d_elements )
    {
        for ( int i = 0; i < 3; ++i )
        {
            // Map local to global and compute u(x).
            d_mesh->mapLocalToGlobalFrame( ksi[i], e.n1, e.n2, x );
            u_ksi = d_u->trueValue( x );

            // Compute uh(ksi)
            d_mesh->shapeFunctionValues( ksi[i], values );
            uh_ksi = d[e.ind1]*values[0] + d[e.ind2]*values[1];

            // Add the approximation to the sum
            sum += pow( u_ksi-uh_ksi, 2 ) * dx * w[i];
        }
    }

    error = pow(sum,.5);

    // Write the error to output file. 
    std::string filename = "results/errors.csv";
    std::ofstream file( filename, std::fstream::in | std::fstream::out | std::fstream::app );
    file << d_n << ',' << d_f_order << ',' << error << std::endl;
}

//---------------------------------------------------------------------------//
// Compute the L2 norm of |u - uh|.
void Solver::checkBehavior( double& nodes,
                            double& midpoints,
                            const std::vector<double>& u,
                            const std::vector<double>& uh )
{
    for ( int i = 0; i < u.size(); i+=2 )
        nodes += pow( u[i] - uh[i] , 2 );
    nodes = sqrt(nodes);

    for ( int i = 1; i < u.size(); i+=2 )
        midpoints += pow( u[i] - uh[i] , 2 );
    midpoints = sqrt(midpoints);

}

//---------------------------------------------------------------------------//
// Write the solutions to output file for plotting. 
void Solver::plot( const std::vector<double>& domain,
                   const std::vector<double>& uh,
                   const std::vector<double>& u )
{

    int n = uh.size();

    // Open new file for this combination of n and f.
    std::string filename = "results/results_" + std::to_string(d_n) +
                           "_" + std::to_string(d_f_order) + ".csv";
    std::ofstream file( filename );

    // Add column names.
    file << "d,uh,u" << std::endl;

    // Write the solutions to the file.
    for ( int i = 0; i < n; ++i )
    {
        file << domain[i] << ","
             << uh[i] << "," << u[i] << std::endl;
    }

}

//---------------------------------------------------------------------------//

} 

//---------------------------------------------------------------------------//
// end Solver.cc
//---------------------------------------------------------------------------//
