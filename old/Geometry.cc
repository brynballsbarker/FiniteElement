//---------------------------------------------------------------------------//
/*!
 * \file Geometry.cc
 */
//---------------------------------------------------------------------------//

#include "Geometry.hh"

#include <cassert>

namespace FEA
{
//---------------------------------------------------------------------------//
// Constructor.
Geometry::Geometry( const std::array<double,2> bounds )
    : d_left( bounds[0] )
    ; d_right( bounds[1] )
{ /* /// */ }

//---------------------------------------------------------------------------//
// Set the initial material id of the geometry.
double Geometry::leftBoundary()
{
    return d_left;
}

//---------------------------------------------------------------------------//
// Set the initial material id of the geometry.
double Geometry::rightBoundary()
{
    return d_right;
}

//---------------------------------------------------------------------------//
// Set the initial velocity field of the geometry.
bool Geometry::valueInGeometry( const double& x ) const
{
    return ( d_left <= x && x <= d_right );
}

//---------------------------------------------------------------------------//

} 

//---------------------------------------------------------------------------//
// end Geometry.cc
//---------------------------------------------------------------------------//
