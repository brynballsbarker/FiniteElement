#include <cassert>
#include <cmath>

namespace FEA
{
//---------------------------------------------------------------------------//
// Constructor.
Function::Function(const int f_order,
                   const double c,
                   const bool u,
                   const bool f )
    : d_f_order( f_order )
    , d_c( c )
    , d_u( u )
    , d_f( f )
{
    /*.......*/
}
//---------------------------------------------------------------------------//
// Get the value of f at a point x.
double Function::evaluate( const double& coord )
{
    if ( d_f )
    {
        if ( d_f_order == 0 )
            return  d_c;
        else if (d_f_order == 1 )
            return coord;
        else
           return coord * coord;
    }
    else
    {
        if ( d_f_order == 0 )
            return  d_c * (1. - pow(coord,2) )/ 2.;
        else if (d_f_order == 1 )
            return (1 - pow(coord,3))/6.;
        else
           return (1 - pow(coord,4))/12.;
    }
}

//---------------------------------------------------------------------------//

} // end namespace FEA

//---------------------------------------------------------------------------//
// end Function.cc
//---------------------------------------------------------------------------//


