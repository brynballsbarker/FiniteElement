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
                const double c,
                const int p,
                const bool beam,
                const std::array<double,5> consts )
    : d_n( n )
    , d_f_order( f_order )
    , d_p( p )
    , d_beam( beam )
    , d_consts( consts )
{
    // Define f.
    d_f = std::make_shared<Function>( f_order, c, false, true, beam, consts );

    // Define the solution.
    d_u = std::make_shared<Function>( f_order, c, true, false, beam, consts );

    // Create the mesh.
    d_mesh = std::make_shared<Mesh>( n, p );
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
    int num_elements = d_mesh->totalNumElements();
    int num_nodes = d_mesh->totalNumNodes();

    int N_size = num_elements + d_p;


    // Storing for plotting
    //
    // domain
    //

    std::string filename = "toplot.csv";
    std::ofstream file( filename );

    // Add column names.
    file << "x,y" << std::endl;

    std::array<double,100> dom;

    int dom_len = dom.size();
    double h_step = 2. / ( dom_len - 1. );

    for ( int i = 0; i < dom_len; ++i )
        dom[i] = -1. + i * h_step;

    for ( auto& e : d_elements )
    {
        std::array<double,100> x;
        int x_len = x.size();
        h_step = (e.n2 - e.n1)/(x_len-1.);

        for ( int i = 0; i < x_len; ++i )
            x[i] = e.n1 + i * h_step;

        for ( auto& s : e.shape_functions )
            for ( int i = 0; i < x_len; ++i )
                file << x[i] << "," << s->evaluate(dom[i]) << std::endl;
        
    }


    // Initialize K.
    std::vector<std::vector<double> > K(
            N_size,std::vector<double>(N_size, 0.0));

    // Initialize F.
    std::vector<double> F(N_size, 0.0);

    // Loop over elements.
    for ( auto& e : d_elements )
    {
        // Add local stiffness matrix to global
        updateK( K, e );

        // Add local force vector to global
        updateF( F, e );

    }

    // Adjust F to match boundary data.
    F[N_size-1] = 0.;

    // Adjust K to match boundary data.
    for ( int i = 0; i < N_size; ++i )
    {
        K[N_size-1][i] = 0.;
        K[i][N_size-1] = 0.;
    }
    K[N_size-1][N_size-1] = 1.0;

    if ( d_beam )
    {
        // Adjust F to match beam boundary data.
        double end = -1.0;
        F[0] -= d_elements[0].shape_functions[0]->evaluateDeriv(end)*d_consts[0];
        F[0] += d_elements[0].shape_functions[0]->evaluate(end)*d_consts[1];
        F[1] -= d_elements[0].shape_functions[1]->evaluateDeriv(end)*d_consts[0];
        F[N_size-2] = 0.;

        // Adjust K to match beam boundary data.
        for ( int i = 0; i < N_size; ++i )
        {
            K[N_size-2][i] = 0.;
            K[i][N_size-2] = 0.;
        }
        K[N_size-2][N_size-2] = 1.0;
    }

    // Solve Kd = F using Gauss Elimination.
    std::vector<double> d( F.begin(), F.end() );;
    gauss( d, K );

    // Initialize domain and solutions.
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

        // Update stiffness matrix.
        for ( int j = 0; j < d_p+1; ++j )
        {
            for ( int k = j; k < d_p+1; ++k )
            {
                // Compute energy inner product.
                double val = energyInnerProduct(e.shape_functions[j],
                                                e.shape_functions[k]);

                // Store in correct K location.
                K[i+j][i+k] += val;
                if ( j != k )
                    K[i+k][i+j] += val;
            }
        }

}

//---------------------------------------------------------------------------//
// Update the global force vector.
void Solver::updateF( std::vector<double>& F,
                      FEA::Element& e) 
{
        // Get global element position
        int i = e.ind1;

        // Update stiffness matrix.
        for ( int j = 0; j < d_p+1; ++j )
            F[i+j] +=  innerProduct(e.shape_functions[j],d_f);

}

//---------------------------------------------------------------------------//
// Use Gauss Elimination to solve the system Kd = F.
void Solver::gauss( std::vector<double>& d,
                    std::vector<std::vector<double> > K ) 
{
    double scale;
    int N_size = d_n + d_p;
    
    // Forward reduction.
    for ( int i = 0; i < N_size-2; ++i )
    {
        // Scale row so first entry is 1.
        scale = K[i][i];
        int end1 = ( i+d_p+1 > N_size ) ? N_size : i+d_p+1;
        for ( int j = i; j < end1; ++j )
            K[i][j] /= scale;
        d[i] /= scale;
        
        // Eliminate other leading entries.
        int end2 = ( d_p+1 > N_size-i ) ? N_size-i : d_p+1;
        for ( int k = 1; k < end2; ++k )
        {
            scale = K[i+k][i];
            int end3 = ( i+d_p+1 > N_size ) ? N_size : i+d_p+1;
            for ( int j = i; j < end3; ++j )
                K[i+k][j] -= scale * K[i][j];
            d[i+k] -= scale * d[i];
        }
    }

    // Handle the second to last row.
    scale = K[N_size-2][N_size-2];
    K[N_size-2][N_size-2] /= scale;
    d[N_size-2] /= scale;

    // Backward reduction. 
    for ( int i = N_size-2; i >= 1; --i )
    {
        int end4 = ( 0 < i-d_p ) ? i-d_p : 0;
        for ( int j = i-1; j >= end4; --j )
        {
            scale = K[j][i];
            K[j][i] -= scale * K[i][i];
            d[j] -= scale * d[i];
        }
    }

}

//---------------------------------------------------------------------------//
// Create the domain which consists of each node and each element midpoint.
void Solver::createDomain( std::vector<double>& domain )
{
    // Determine h.
    int dom_len = domain.size();
    double h_step = 1. / ( dom_len - 1. );

    // Initialize the domain.
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

    // Initialize necessary variables.
    double x;
    int element_id;
    double ref_coord;
    std::array<double,2> values;

    // Loop through each value in domain
    for ( int i = 0; i < dom_len; ++i )
    {
        // Get relevant domain value.
        x = domain[i];

        // Determine which element x is in.
        d_mesh->locateX( x, element_id );

        // Find local position of x in the given element.
        d_mesh->mapGlobalToLocalFrame( x, element_id, ref_coord );

        // Compute the value of uh at x.
        // Get the element x is contained in.
        auto& e = d_elements[ element_id ];

        // Loop through the relevant shape functions to compute uh(x).
        double val = 0.;
        for ( int j = 0; j < d_p + 1; ++j )
            val += d[element_id+j] * e.shape_functions[j]->evaluate(ref_coord);
        uh[i] = val;    
    }
}

//---------------------------------------------------------------------------//
// Evalute the true solution over the domain.
void Solver::evaluateTrue( const std::vector<double>& domain, 
                           std::vector<double>& u )
{
    int dom_len = domain.size();

    // Evaluate u at each point in the domain.
    for ( int i = 0; i < dom_len; ++i )
        u[i] = d_u->evaluate( domain[i] );
}

//---------------------------------------------------------------------------//
// Compute the L2 norm of |u - uh|.
void Solver::computeError( double& error,
                    const std::vector<double>& d )
{
    double u_xi;
    double uh_xi;

    double x;
    std::array<double,2> values;

    // Define gaussian quadrature interpolation points. 
    std::array<double,3> xi = {-pow(3./5.,.5), 0.0, pow(3./5.,.5)};
    std::array<double,3> w = {5./9.,8./9.,5./9.};

    // Compute jacobian J.
    double J = d_mesh->getJacobian();

    // Initialize sum. 
    double sum = 0.0;

    // Loop over elements.
    for ( auto& e : d_elements )
    {
        for ( int i = 0; i < 3; ++i )
        {
            // Map local to global and compute u(x).
            d_mesh->mapLocalToGlobalFrame( xi[i], e.n1, e.n2, x );
            u_xi = d_u->evaluate( x );

            // Compute uh(xi)
            double val = 0.;
            for ( int j = 0; j < d_p + 1; ++j )
                val += d[e.ind1+j] * e.shape_functions[j]->evaluate(xi[i]);
            uh_xi = val;    

            // Add the approximation to the sum
            sum += pow( u_xi-uh_xi, 2 ) * J * w[i];
        }
    }

    error = pow(sum,.5);

    // Write the error to output file. 
    std::string filename = "results/errors.csv";
    std::ofstream file( filename, std::fstream::in | std::fstream::out | std::fstream::app );
    file << d_n << ',' << d_f_order << ',' << error << std::endl;
}

//---------------------------------------------------------------------------//
// Compute the inner product (f,N_A).
double Solver::innerProduct( std::shared_ptr<FEA::Shape> func1, 
                             std::shared_ptr<FEA::Function> func2 )
{
    std::array<double,2> endpoints = func1->getEndpoints();

    double x;
    double val1, val2;

    // Define gaussian quadrature interpolation points. 
    std::array<double,3> xi = {-pow(3./5.,.5), 0.0, pow(3./5.,.5)};
    std::array<double,3> w = {5./9.,8./9.,5./9.};

    // Compute jacobian J.
    double J = d_mesh->getJacobian();

    // Initialize sum. 
    double sum = 0.0;

    // Loop over interpolation points.
    for ( int i = 0; i < 3; ++i )
    {
        // Evaluate the shape function.
        val1 = func1->evaluate(xi[i]);

        // Evaluate the force function.
        d_mesh->mapLocalToGlobalFrame( xi[i], endpoints[0], endpoints[1], x );
        val2 = func2->evaluate(x);

        sum += val1*val2*J*w[i];
    }

    return sum;
    
}
//---------------------------------------------------------------------------//
// Compute the energy inner product (N_A, N_B).
double Solver::energyInnerProduct( std::shared_ptr<FEA::Shape> func1, 
                                   std::shared_ptr<FEA::Shape> func2)
{
    double val1, val2;

    // Define gaussian quadrature interpolation points. 
    std::array<double,3> xi = {-pow(3./5.,.5), 0.0, pow(3./5.,.5)};
    std::array<double,3> w = {5./9.,8./9.,5./9.};

    // Compute jacobian J.
    double J = d_mesh->getJacobian();

    // Initialize sum. 
    double sum = 0.0;

    // Loop over interpolation points.
    for ( int i = 0; i < 3; ++i )
    {
        if ( d_beam )
        {
            // Store E*I.
            double EI = d_consts[3]*d_consts[4];

            // Evaluate both shape functions.
            val1 = func1->evaluateSecondDeriv(xi[i]);
            val2 = func2->evaluateSecondDeriv(xi[i]);
            sum += val1*val2*J*w[i]*EI;
        }
        else
        {
            // Evaluate both shape functions.
            val1 = func1->evaluateDeriv(xi[i]);
            val2 = func2->evaluateDeriv(xi[i]);
            sum += val1*val2*J*w[i];
        }
    }
    return sum;
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
        file << domain[i] << ","
             << uh[i] << "," << u[i] << std::endl;

}

//---------------------------------------------------------------------------//

} 

//---------------------------------------------------------------------------//
// end Solver.cc
//---------------------------------------------------------------------------//
