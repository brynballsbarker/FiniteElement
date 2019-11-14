//---------------------------------------------------------------------------//
/*!
 * \file Mesh.cc
 */
//---------------------------------------------------------------------------//

#include "Mesh.hh"
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>

namespace FEA
{

//---------------------------------------------------------------------------//
// Constructor.
Mesh::Mesh( const int n, const int p )
    : d_num_nodes( n+1 )
    , d_num_elements( n )
    , d_element_width( 1./n )
    , d_shape_order( p )
{ 
    // Create the boundary nodes.
    d_boundary_nodes.resize(1);
    d_boundary_nodes[0] = d_num_nodes;

    // Create bernstein basis.
    d_bernstein.resize( d_shape_order+1 );
    for ( int a = 1; a <= d_shape_order + 1; ++a )
        d_bernstein[a-1] = std::make_shared<Bernstein>( d_shape_order, a );

/*
    // Write bernstein basis to file to check plot.
    std::string filename = "toplotbern.csv";
    std::ofstream file( filename );

    // Add column names.
    file << "x,y" << std::endl;

    // Create domain.
    std::array<double,100> dom;
    int dom_len = dom.size();
    double h_step = 2. / ( dom_len - 1. );
    for ( int i = 0; i < dom_len; ++i )
        dom[i] = -1. + i * h_step;

    // Write basis values for each basis function to file.
    for ( auto& b : d_bernstein  )
        for ( int i = 0; i < dom_len; ++i )
            file << dom[i] << "," << b->evaluate(dom[i]) << std::endl;
*/

    // Allocate storage for C matrices
    d_coefs.resize( d_num_elements );

    // C0 basis.
    if ( d_shape_order == 1 )
    { 
        // All elements have same C matrix.
        for ( int e = 0; e < d_num_elements; ++e )
        {
            d_coefs[e].resize( d_shape_order + 1 );

            std::vector<double> vect{ 1., 0 };
            d_coefs[e][0] = vect;

            vect = { 0., 1.};
            d_coefs[e][1] = vect;

        }
    }
    // C1 basis
    if ( d_shape_order == 2 )
    {
        // First and last elements have unique C.
        for ( int e = 0; e < d_num_elements; ++e )
        {
            d_coefs[e].resize( d_shape_order + 1 );
            if ( e == 0 )
            {
                std::vector<double> vect{ 1., 0., 0. };
                d_coefs[e][0] = vect;

                vect = { 0., 1., .5 };
                d_coefs[e][1] = vect;

                vect = {0.,0.,.5};
                d_coefs[e][2] = vect;
            }
            else if ( e == ( d_num_elements - 1 ) )
            {
                std::vector<double> vect{ .5, 0., 0. };
                d_coefs[e][0] = vect;

                vect = { .5, 1., 0. };
                d_coefs[e][1] = vect;

                vect = {0.,0.,1.};
                d_coefs[e][2] = vect;
            }
            else
            {
                std::vector<double> vect{ .5, 0., 0. };
                d_coefs[e][0] = vect;

                vect = { .5, 1., .5 };
                d_coefs[e][1] = vect;

                vect = {0.,0.,.5};
                d_coefs[e][2] = vect;
            }
        }
    }
    // C2 basis.
    if ( d_shape_order == 3 )
    {
        // First and last two have unique C.
        for ( int e = 0; e < d_num_elements; ++e )
        {
            d_coefs[e].resize( d_shape_order + 1 );
            if ( e == 0 )
            {
                std::vector<double> vect{ 1., 0., 0., 0. };
                d_coefs[e][0] = vect;

                vect = { 0., 1., .5, .25 };
                d_coefs[e][1] = vect;

                vect = {0.,0.,.5, 7./12.};
                d_coefs[e][2] = vect;

                vect = {0.,0.,0.,1./6.};
                d_coefs[e][3] = vect;
            }
            else if ( e == 1 )
            {
                std::vector<double> vect{ .25, 0., 0., 0. };
                d_coefs[e][0] = vect;

                vect = { 7./12., 2./3., 1./3., 1./6. };
                d_coefs[e][1] = vect;

                vect = {1./6.,1./3.,2./3., 2./3.};
                d_coefs[e][2] = vect;

                vect = {0.,0.,0.,1./6.};
                d_coefs[e][3] = vect;
            }
            else if ( e == ( d_num_elements - 2 ) )
            {
                std::vector<double> vect{ 1./6., 0., 0.,0. };
                d_coefs[e][0] = vect;

                vect = {2./3.,2./3.,1./3.,1./6.};
                d_coefs[e][1] = vect;

                vect = {1./6.,1./3.,2./3.,7./12.};
                d_coefs[e][2] = vect;

                vect = {0.,0.,0.,1./4.};
                d_coefs[e][3] = vect;
            }
            else if ( e == ( d_num_elements - 1 ) )
            {
                std::vector<double> vect{ 1./6.,0., 0., 0. };
                d_coefs[e][0] = vect;

                vect = {7./12.,1./2.,0.,0.};
                d_coefs[e][1] = vect;

                vect = {.25,.5,1.,0.};
                d_coefs[e][2] = vect;

                vect = {0.,0.,0.,1.};
                d_coefs[e][3] = vect;
            }
            else
            {
                std::vector<double> vect{ 1./6., 0., 0.,0. };
                d_coefs[e][0] = vect;

                vect = {2./3.,2./3.,1./3.,1./6.};
                d_coefs[e][1] = vect;

                vect = {1./6.,1./3.,2./3.,2./3.};
                d_coefs[e][2] = vect;

                vect = {0.,0.,0.,1./6.};
                d_coefs[e][3] = vect;
            }
        }
    }

}

//---------------------------------------------------------------------------//
// Get the total number of elements in the mesh.
int Mesh::totalNumElements() const
{
    return d_num_elements;
}

//---------------------------------------------------------------------------//
// Get the total number of nodes in the mesh.
int Mesh::totalNumNodes() const
{
    return d_num_nodes;
}

//---------------------------------------------------------------------------//
// Get the width of each element in the mesh.
double Mesh::getElementWidth() const
{
    return d_element_width;
}

//---------------------------------------------------------------------------//
// Get the Jacobian for the global to local mapping.
double Mesh::getJacobian() const
{
    return d_element_width/2.;
}

//---------------------------------------------------------------------------//
// Given a node id get its coordinate.
void Mesh::nodeCoordinate( const int node_id,
                            double& coord ) const
{
    assert( node_id < totalNumNodes() );

    // Get the coordinate.
    coord = node_id * d_element_width;
}

//---------------------------------------------------------------------------//
// Given an element id get its left and right endpoints. 
void Mesh::elementCoordinates( const int element_id,
                              std::array<double,2>& coords ) const
{
    assert( element_id < d_num_elements );

    // Initialize the coordinates.
    double n1;
    double n2;

    // Get the node coordinates.
    nodeCoordinate(element_id, n1);
    nodeCoordinate(element_id+1, n2);

    // Assign them to coords.
    coords[0] = n1;
    coords[1] = n2;
}

//---------------------------------------------------------------------------//
// Given a cardinal element id intitalize the element.
void Mesh::initializeElement(
    const int element_id,
    Element& element ) const
{
    
    // Save node and element ids
    std::array<double,2> coords;
    elementCoordinates( element_id, coords );

    // Establish basic characteristics.
    element.id = element_id;
    element.n1 = coords[0];
    element.n2 = coords[1];
    element.ind1 = element_id;
    element.ind2 = element_id + 1;
    
    // Initialize shape functions.
    element.shape_functions.resize( d_shape_order+1 );
    for ( int a = 1; a <= d_shape_order+1; ++a )
        element.shape_functions[a-1] = std::make_shared<FEA::Shape>(
                d_coefs[element_id][a-1], d_bernstein, coords );

}

//---------------------------------------------------------------------------//
// Given a point in the domain, locate the cardinal index of the element
// in which it is located
void Mesh::locateX( const double& x,
                    int& element_id ) const
{
    // Mesh spans (0.0, 1.0)
    element_id = std::floor( x / d_element_width );

    // Check for last element.
    if ( x == 1. )
        element_id = d_num_elements - 1;

    assert( 0 <= element_id );
    assert( element_id < d_num_elements );
}

//---------------------------------------------------------------------------//
// Map the value x from the global frame to the local frame of the element
// in which it is located.
void Mesh::mapGlobalToLocalFrame(
    const double& x,
    const int& element_id,
    double& ref_coord ) const
{
    assert( element_id == std::floor( x / d_element_width ) || element_id == (d_num_elements-1));

    // The reference cell spans -1 to 1 in the x direction 
    // and 0 to 1 in the y direction.
    //
    //     (-1,1)------------(1,1) 
    //       |                |
    //       |                |
    //       |                |
    //       |                |
    //     (-1,0)-----------(1,0)  

    ref_coord = ( x / d_element_width - element_id ) * 2. - 1.;

    assert( -1.0 <= ref_coord && ref_coord <= 1.0 );
}

//---------------------------------------------------------------------------//
// Map a reference coordinate in an element from its local frame to 
// the global frame.
void Mesh::mapLocalToGlobalFrame(
    const double& ref_coord,
    const double& n1, const double& n2,
    double& x ) const
{
    x = ( ( 1. - ref_coord ) * n1 + ( 1. + ref_coord ) * n2 ) / 2.;
}
//---------------------------------------------------------------------------//

}

//---------------------------------------------------------------------------//
// end Mesh.cc
//---------------------------------------------------------------------------//
