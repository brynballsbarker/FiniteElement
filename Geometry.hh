//---------------------------------------------------------------------------//
/*!
 * \file Geometry.hh
 */
//---------------------------------------------------------------------------//

#ifndef FEA_GEOMETRY_HH
#define FEA_GEOMETRY_HH

#include "Element.hh"

#include <array>
#include <functional>

namespace FEA
{
//---------------------------------------------------------------------------//
/*!
 * \class Geometry
 *
 * \brief Geometry interface for the creation of initial particle state.
 */
class Geometry
{
  public:

    // Constructor.
    Geometry( const std::array<double,2> bounds );

    double leftBoundary() const;

    double rightBoundary() const;

    // Determine if a particle is in the geometry.
    bool valueInGeometry( const double& x ) const;

  private:

    // Left boundary.
    double d_left;

    // Right boundary.
    double d_right;
};

//---------------------------------------------------------------------------//

} // end namespace FEA

#endif // end FEA_GEOMETRY_HH

//---------------------------------------------------------------------------//
// end Geometry.hh
//---------------------------------------------------------------------------//
