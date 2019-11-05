//---------------------------------------------------------------------------//
/*!
 * \file Shape.cc
 */
//---------------------------------------------------------------------------//

#include "Shape.hh"

#include <array>
#include <vector>
#include <cmath>

namespace FEA
{

//---------------------------------------------------------------------------//
// Constructor.
Shape::Shape( const std::vector<double> coef,
              std::vector<std::shared_ptr<Bernstein> > basis,
              const std::array<double, 2> endpoints )
    : d_coefs( coef )
    , d_basis( basis )
    , d_endpoints( endpoints )
{ 
    // Compute p choose a.
    d_local = true;
}

//---------------------------------------------------------------------------//
// Get the value of u at a point x.
double Shape::evaluate( const double& coord )
{   
    double val = 0;

    for ( std::size_t i = 0; i < d_coefs.size(); ++i )
        val += d_coefs[i] * d_basis[i]->evaluate( coord );

    return val
}

//---------------------------------------------------------------------------//
// Get the value of u at a point x.
double Shape::evaluateDeriv( const double& coord )
{   
    double val = 0;

    for ( std::size_t i = 0; i < d_coefs.size(); ++i )
        val += d_coefs[i] * d_basis[i]->evaluateDeriv( coord );

    return val
}

//---------------------------------------------------------------------------//
// Get the value of u at a point x.
bool Shape::checkLocal()
{   
    return d_local;
}

//---------------------------------------------------------------------------//

} 

//---------------------------------------------------------------------------//


//---------------------------------------------------------------------------//
// end Shape.cc
//---------------------------------------------------------------------------//
