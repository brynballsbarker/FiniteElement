//---------------------------------------------------------------------------//
/*!
 * \file Solver.hh
 */
//---------------------------------------------------------------------------//

#ifndef FEA_SOLVER_HH
#define FEA_SOLVER_HH

#include "Mesh.hh"
#include "Element.hh"
#include "Function.hh"

#include <vector>
#include <array>

namespace FEA
{
//---------------------------------------------------------------------------//
/*!
 * \class Solver
 * \brief Simulation manager.
 */
class Solver
{
  public:

    // Constructor.
    Solver( const int num_elements,
            const int f_order,
            const double c,
            const int p );

    // Initialize the problem by initializing the elements.
    void initialize();

    // Solve for the FEA Approximate Solution.
    void solve();

  private:

    // Update the global stiffness matrix.
    void updateK( std::vector<std::vector<double> >& K,
                  FEA::Element& e );

    // Update the global force vector.
    void updateF( std::vector<double>& F,
                  FEA::Element& e);

    // Compute the inner products (N_1,f) and (N_2,f).
    void integrateRHS( const FEA::Element& e,
                       std::array<double,2>& f );

    // Use Gauss Elimination to solve the system Kd = F.
    void gauss( std::vector<double>& d,
                std::vector<std::vector<double> > K );

    // Create the domain which consists of each node 
    // and each element midpoint.
    void createDomain( std::vector<double>& domain );

    // Evaluate the FEA solution over the domain.
    void evaluateApprox( const std::vector<double>& domain,
                         const std::vector<double>& d,
                         std::vector<double>& uh );

    // Evaluate the true solution over the domain.
    void evaluateTrue( const std::vector<double>& domain,
                         std::vector<double>& u );

    // Compute the L2 norm of |u - uh|.
    void computeError( double& error,
                       const std::vector<double>& d );

    // Compute the L2 norm of |u - uh|.
    double innerProduct( std::shared_ptr<FEA::Shape> func1, 
                         std::shared_ptr<FEA::Function> func2 );
   
    // Compute the L2 norm of |u - uh|.
    double energyInnerProduct( std::shared_ptr<FEA::Shape> func1, 
                               std::shared_ptr<FEA::Shape> func2 );

    // Look at behavior on the nodes and on the interior.
    void checkBehavior( double& nodes,
                            double& midpoints,
                            const std::vector<double>& u,
                            const std::vector<double>& uh );

    // Write the solutions to output files for plotting..
    void plot( const std::vector<double>& domain,
               const std::vector<double>& uh,
               const std::vector<double>& u );

  private:

    // n value.
    int d_n;

    // Order of f.
    int d_f_order;

    // Shape funciton order.
    int d_p;

    // Focre function;
    std::shared_ptr<Function> d_f;

    // True solution.
    std::shared_ptr<Function> d_u;

    // Mesh.
    std::shared_ptr<Mesh> d_mesh;

    // Elements.
    std::vector<Element> d_elements;

};

//---------------------------------------------------------------------------//

}

#endif 

//---------------------------------------------------------------------------//
// end Solver.hh
//---------------------------------------------------------------------------//
