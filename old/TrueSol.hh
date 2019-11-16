//---------------------------------------------------------------------------//
/*!
 * \file TrueSol.hh
 */
//---------------------------------------------------------------------------//

#ifndef FEA_TRUESOL_HH
#define FEA_TRUESOL_HH

#include "Function.hh"

#include <array>
#include <vector>

namespace FEA
{

//---------------------------------------------------------------------------//
/*!
 * \class TrueSol
 *
 * \brief True solution to system.
 */
class TrueSol : public Function
{
  public:

    // Constructor.
    TrueSol( const int f_order,
             const double c );

    // Get the value of u at a point x.
    double evaluate( const double& coord ) const override;
    
  private:

    // Material id.
    int d_f_order;

    // Color
    int d_c;

};

} 

//---------------------------------------------------------------------------//

#endif // end FEA_TRUESOL_HH

//---------------------------------------------------------------------------//
// end TrueSol.hh
//---------------------------------------------------------------------------//
