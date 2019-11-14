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
Shape::Shape( std::vector<double> coef,
              std::vector<std::shared_ptr<FEA::Bernstein> > basis,
              std::array<double, 2> endpoints )
    : d_coefs( coef )
    , d_basis( basis )
    , d_endpoints( endpoints )
{ /********/ }

//---------------------------------------------------------------------------//
// Get the value of shape func at a point x.
double Shape::evaluate( double& coord )
{   
    double val = 0;

    // Linear combination of bernstein basis and coefs.
    for ( std::size_t i = 0; i < d_coefs.size(); ++i )
        val += d_coefs[i] * d_basis[i]->evaluate( coord );

    return val;
}

//---------------------------------------------------------------------------//
// Get the derivative of shape func at a point x.
double Shape::evaluateDeriv( double& coord )
{   
    double val = 0;

    // Linear combination of bernstein basis and coefs.
    for ( std::size_t i = 0; i < d_coefs.size(); ++i )
        val += d_coefs[i] * d_basis[i]->evaluateDeriv( coord );

    // Scale by inverse of mapping Jacobian.
    double scale = 2./(d_endpoints[1]-d_endpoints[0]);
    return val*scale;
}

//---------------------------------------------------------------------------//
// Get the second derivative of shape func at a point x.
double Shape::evaluateSecondDeriv( double& coord )
{   
    double val = 0;

    // Linear combination of bernstein basis and coefs.
    for ( std::size_t i = 0; i < d_coefs.size(); ++i )
        val += d_coefs[i] * d_basis[i]->evaluateSecondDeriv( coord );

    // Scale by inverse of mapping Jacobian squared.
    double scale = pow(2./(d_endpoints[1]-d_endpoints[0]),2);
    return val*scale;
}

//---------------------------------------------------------------------------//
// Get the value of u at a point x.
std::array<double,2> Shape::getEndpoints() const
{   
    return d_endpoints;
}

//---------------------------------------------------------------------------//

} 

//---------------------------------------------------------------------------//


//---------------------------------------------------------------------------//
// end Shape.cc
//---------------------------------------------------------------------------//
