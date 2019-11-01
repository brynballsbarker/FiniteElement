//---------------------------------------------------------------------------//
/*!
 * \file Element.hh
 */
//---------------------------------------------------------------------------//

#ifndef FEA_ELEMENT_HH
#define FEA_ELEMENT_HH

#include <array>

namespace FEA
{

//---------------------------------------------------------------------------//
/*!
 * \class Element
 */
class Element
{
  public:

    //@{
    //! ID.
    int id;

    //! Left endpoint.
    double n1;

    //! Right endpoint.
    double n2;

    //! Left node id.
    int ind1;

    //! Right node id.
    int ind2;

    //! Element basis functions.
    std::array<double,2> basis_values;

    //! Element basis derivatives.
    std::array<double,2> basis_derivatives;;

    //! Element stiffness matrix. 
    std::array<std::array<double,2>,2> k;
    //@}

    // Constructor.
    Element()
    { /* //// */ }
};

//---------------------------------------------------------------------------//

} 

#endif // end FEA_ELEMENT_HH

//---------------------------------------------------------------------------//
// end Element.hh
//---------------------------------------------------------------------------//
