//---------------------------------------------------------------------------//
/*!
 * \file Force.hh
 */
//---------------------------------------------------------------------------//

#ifndef FEA_FORCE_HH
#define FEA_FORCE_HH

#include "Function.hh"

#include <array>
#include <vector>

namespace FEA
{

//---------------------------------------------------------------------------//
/*!
 * \class Force
 *
 * \brief Force function. 
 */
class Force : public Function
{
  public:

    // Constructor.
    Force( const int f_order,
           const double c );

    // Get the value of f at a point x.
    double evaluate( const double& coord ) const override;
         
  private:

    // Material id.
    int d_f_order;

    // Color
    int d_c;

    int d_u;

    int d_f;

};

} 

//---------------------------------------------------------------------------//

#endif // end FEA_FORCE_HH

//---------------------------------------------------------------------------//
// end Force.hh
//---------------------------------------------------------------------------//
