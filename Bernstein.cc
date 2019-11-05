//---------------------------------------------------------------------------//
/*!
 * \file Bernstein.cc
 */
//---------------------------------------------------------------------------//

#include "Bernstein.hh"
#include "Functions.hh"

#include <array>
#include <vector>
#include <cmath>

namespace FEA
{

//---------------------------------------------------------------------------//
// Constructor.
Bernstein::Bernstein( const int p_val,
                  const int a_val )
    : d_p( p_val )
    , d_a( a_val )
{ 
    // Compute p choose a.
    d_choose = 1.0 * factorial(d_p) / (1.0 * factorial( d_a - 1) * factorial( d_p + 1 - d_a ));

    d_local = true;
}

//---------------------------------------------------------------------------//
// Get the value of u at a point x.
double Bernstein::evaluate( const double& coord )
{   
    return 1./pow(2,d_p) * d_choose * pow(1.-coord,d_p+1-d_a) * pow(1.+coord,d_a-1)
}

//---------------------------------------------------------------------------//
// Get the value of u at a point x.
double Bernstein::evaluateDeriv( const double& coord )
{   
    double coef = 1./pow(2,d_p) * d_choose;
    double part1 = pow(1.-coord,d_p+1-d_a);
    double part2 = pow(1.+coord,d_a-1);
    double part1_prime, part2_prime;

    if ( d_p + 1 == d_a )
        part1_prime = 0.0;
    else if ( d_p == d_a )
        part1_prime = -1.0;
    else 
        part1_prime = -1.0*(d_p-d_a+1.)*pow(1.-coord,d_p-d_a);

    if ( d_a == 1 )
        part2_prime = 0.0;
    else if ( d_a == 2 )
        part2_prime = 1.0;
    else
        part2_prime = (d_a-1.)* pow(1.+coord,d_a-2);

    return coef * (part1_prime*part2 + part1*part2_prime);
}

//---------------------------------------------------------------------------//

} 

//---------------------------------------------------------------------------//


//---------------------------------------------------------------------------//
// end Bernstein.cc
//---------------------------------------------------------------------------//
