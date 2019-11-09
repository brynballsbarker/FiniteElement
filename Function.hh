//---------------------------------------------------------------------------//
/*!
 * \file Function.hh
 */
//---------------------------------------------------------------------------//

#ifndef FEA_FUNCTION_HH
#define FEA_FUNCTION_HH


namespace FEA
{
//---------------------------------------------------------------------------//
/*!
 * \class Function
 *
 * \brief Function interface for the creation of initial particle state.
 */
class Function
{
  public:

    // Constructor.
    Function(const int f_order,
             const double c,
             const bool u,
             const bool f );

    // Get the value of f at a point x.
    double evaluate( const double& coord );

  private:

    // Material id.
    int d_f_order;

    // Color
    int d_c;

    bool d_u;

    bool d_f;
};

//---------------------------------------------------------------------------//

} // end namespace FEA

#endif // end FEA_FUNCTION_HH

//---------------------------------------------------------------------------//
// end Function.hh
//---------------------------------------------------------------------------//
