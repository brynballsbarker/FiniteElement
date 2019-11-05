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
 * \brief True solution to system.
 */
class Shape
{
  public:

    // Constructor.
    Shape( const std::vector<double> coef,
           std::vector<std::shared_ptr<Bernstein> > basis,
           const std::array<double, 2> endpoints );

    // Get the value of shpae func at a point x.
    double evaluate( const double& coord );

    // Get the derivate at point x.
    double evaluateDeriv( const double& coord );

    // 
    bool checkLocal();

  private:

    // Order of f.
    std::vector<double> d_coefs;

    // Bernstein basis
    std::vector<std::shared_ptr<Bernstein> > d_basis;

    // Defined on local
    bool d_local;

    std::array<double,2> d_endpoints;

};

} 

//---------------------------------------------------------------------------//

#endif // end FEA_SHAPE_HH

//---------------------------------------------------------------------------//
// end Shape.hh
//---------------------------------------------------------------------------//
