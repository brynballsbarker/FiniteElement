//---------------------------------------------------------------------------//
/*!
 * \file Shape.hh
 */
//---------------------------------------------------------------------------//

#ifndef FEA_SHAPE_HH
#define FEA_SHAPE_HH

#include "Bernstein.hh"

#include <array>
#include <vector>

namespace FEA
{

//---------------------------------------------------------------------------//
/*!
 * \class Shape
 *
 * \brief Shape function for one element.
 */
class Shape
{
  public:

    // Constructor.
    Shape( std::vector<double> coef,
           std::vector<std::shared_ptr<FEA::Bernstein> > basis,
           std::array<double, 2> endpoints );

    // Get the value of shape func at a point x.
    double evaluate( double& coord );

    // Get the derivate at point x.
    double evaluateDeriv( double& coord );

    // Get the second derivate at point x.
    double evaluateSecondDeriv( double& coord );

    // Get the endpoints of the correpsonding element.
    std::array<double,2> getEndpoints() const;

  private:

    // C matrix for bezier extraction.
    std::vector<double> d_coefs;

    // Bernstein basis
    std::vector<std::shared_ptr<FEA::Bernstein> > d_basis;

    // Endpoints for corresponding element.
    std::array<double,2> d_endpoints;

};

} 

//---------------------------------------------------------------------------//

#endif // end FEA_SHAPE_HH

//---------------------------------------------------------------------------//
// end Shape.hh
//---------------------------------------------------------------------------//
