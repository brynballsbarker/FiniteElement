//---------------------------------------------------------------------------//
/*!
 * \file Function.cc
 */
//---------------------------------------------------------------------------//

#include "Function.hh"

#include <cassert>
#include <cmath>

namespace FEA
{
//---------------------------------------------------------------------------//
// Constructor.
Function::Function(const int f_order,
                   const double c,
                   const bool u,
                   const bool f,
                   const bool beam,
                   std::array<double,5> consts )
    : d_f_order( f_order )
    , d_c( c )
    , d_u( u )
    , d_f( f )
    , d_beam( beam )
    , d_M( consts[0] )
    , d_Q( consts[1] )
    , d_r( consts[2] )
    , d_E( consts[3] )
    , d_I( consts[4] )
{
    /*.......*/
}
//---------------------------------------------------------------------------//
// Get the value of f at a point x.
double Function::evaluate( const double& coord )
{
    // Evaluate Force function.
    if ( d_f )
    {
        // beam problem.
        if ( d_beam )
            return d_r;

        // f(x) = c.
        else if ( d_f_order == 0 )
            return  d_c;

        // f(x) = x.
        else if (d_f_order == 1 )
            return coord;
         
        // f(x) = x^2;
        else
           return coord * coord;
    }

    // Evaluate true solution.
    else
    {
        // beam problem.
        if ( d_beam )
            return ( 1/24.*d_r*(pow(coord,4)-4*coord+3)
                    +1./6.*d_Q*(pow(coord,3)-3*coord+2)
                    +1./2.*d_M*(pow(coord,2)-2*coord+1) )/(d_E*d_I);

        // u_xx = c.
        else if ( d_f_order == 0 )
            return  d_c * (1. - pow(coord,2) )/ 2.;

        // u_xx = x.
        else if (d_f_order == 1 )
            return (1 - pow(coord,3))/6.;

        // u_xx = x^2.
        else
           return (1 - pow(coord,4))/12.;
    }
}

//---------------------------------------------------------------------------//

} // end namespace FEA

//---------------------------------------------------------------------------//
// end Function.cc
//---------------------------------------------------------------------------//
