//---------------------------------------------------------------------------//
/*!
 * \file Function.hh
 */
//---------------------------------------------------------------------------//

#ifndef FEA_FUNCTION_HH
#define FEA_FUNCTION_HH


#include <array>
#include <vector>

namespace FEA
{

//---------------------------------------------------------------------------//
/*!
 * \class Function
 *
 * \brief Function function. 
 */
class Function 
{
  public:

    // Constructor.
    Function( const int f_order,
              const double c,
              const bool u,
              const bool f,
              const bool beam,
              const std::array<double,5> consts );

    // Get the value of func at a point x.
    double evaluate( const double& coord );
         
  private:

    // For u_xx = f, f order.
    int d_f_order;

    // For u_xx = c, c value.
    int d_c;

    // Is this function u?
    bool d_u;

    // Is this function f?
    bool d_f;

    // Is this the beam problem?
    bool d_beam;

    // M val for beam problem.
    double d_M;

    // Q val for beam problem.
    double d_Q;

    // f val for beam problem.
    double d_r;

    // I val for beam problem.
    double d_I;

    // E val for beam problem.
    double d_E;

};

} 

//---------------------------------------------------------------------------//

#endif // end FEA_FUNCTION_HH

//---------------------------------------------------------------------------//
// end Function.hh
//---------------------------------------------------------------------------//
