//-----------------------------------------------------------------------------
// Module C
//
// Description - Contains the diagonalise function which passes the Hamiltonian
//               matrix for diagonalisation to obtain the eigenvectors and
//               eigenvalues. The vector containing the occupation numbers are
//               also filled.
//-----------------------------------------------------------------------------
#ifndef MODULEC_HPP
#define MODULEC_HPP

#include "../atom/atom.h"
#include "armadillo"
#include "../atom/vect3.h"
#include "../sim/unitcell.h"

using namespace moldycat;
using namespace arma;

// Tolerance to determine whether levels are degenerate
const double eval_tol = 1e-10;

//----------------------------------------------------------------------------------
// Function diagonalise
//----------------------------------------------------------------------------------
// Returns: energy - Total electronic energy
//
// Params: H (mat&)                  - Matrix containing the Hamiltonian
//         evals (vec&)              - Vector containing the eigenvalues
//         evects (mat&)             - Matrix containing the eigenvectors
//         occupation (vec&)         - Vector containing the occupation numbers
//         max_occupied_index (int&) - Index of highest level that is occupied
//
// Operation: Diagonalises the hamiltonian to populate the eigenvalues and
//            eigenvectors. Also takes into account degeneracy to populate the 
//            occupation vector (length 4N), which will be required to work out 
//            the correct hamiltonian forces.
//----------------------------------------------------------------------------------
double diagonalise( mat& H, vec& evals, mat& evects, vec& occupation,int& max_occupied_index)
{
  //Diagonalise the Hamiltonian, get eigenvals and evects
  eig_sym(evals,evects,H); 
  int i,n=occupation.n_elem;
  double energy = 0.0;
    
  //Sum occupied energy levels and reset occupations to 2
  //below fermi level (n/2) and 0.0 above
  for(i=0;i<n/2;i++){
    energy += 2.0*evals(i);
  }

  int level=n/2;
  int degeneracy_above = 0;
  int degeneracy_below = 1;
  double fermi_energy = evals(n/2-1);

  //Look for degeneracy
  while ((level<n) && (abs(evals(level) - fermi_energy) < eval_tol) ) {
    degeneracy_above += 1;
    level += 1;
  }
  if (degeneracy_above != 0) {
    level = n/2 - 2;
    while ((level >= 0) && (abs(evals(level) - fermi_energy) < eval_tol) ) {
      degeneracy_below += 1;
      level -= 1;
    }
  }

  max_occupied_index = n/2 + degeneracy_above;
    
  //--------------------------------------------------------------------------
  // Populate occupation vector.
  //
  // E.g. if the following levels are degenerate
  //
  //   Level index
  //-------> n/2 + 2
  //-------> n/2 + 1
  //-------> n/2
  //-------> n/2 - 1 (highest occupied level without degeneracy)
  //-------> n/2 - 2
  //
  // Then degeneracy_above = 3
  //      degeneracy_below = 2 (total degenerate levels = 5)
  // So for energy levels below n/2 - 2 ( = n/2 - degeneracy_below)
  //   --> occupation number - 1
  // For energy levels between n/2 - 2 and n/2 + 2 inclusive
  //   --> occupation number 
  //        = 2/5 ( = degeneracy_below / (degeneracy_above + degeneracy_below)
  // For energy levels above n/2 + 2 ( = n/2-1 + degeneracy_above)
  //   --> occupation number = 0
  //--------------------------------------------------------------------------
  for(int j=0;j<n;j++){
    if (degeneracy_above != 0) {
      if (j < n/2-degeneracy_below) {
        occupation(j) = 1.0;
      }
      else if (j > n/2-1 + degeneracy_above) {
        occupation(j) = 0.0;
      }
      else {
        occupation(j) = 1.0*degeneracy_below / (degeneracy_above + degeneracy_below);
      }
    }
    else {
      if(j<n/2){
        occupation(j)=1.0;
      }
      else{
        occupation(j)=0.0;
      }
    }
  }
  return energy;
}

//----------------------------------------------------------------------------------
// Function diagonalise
//----------------------------------------------------------------------------------
// Returns: energy - Total electronic energy
//
// Params: H (mat&)                  - Matrix containing the Hamiltonian
//         evals (vec&)              - Vector containing the eigenvalues
//         evects (mat&)             - Matrix containing the eigenvectors
//         occupation (vec&)         - Vector containing the occupation numbers
//
// Operation: Same as above but without max_occupied_index. Used when only the 
//            energy is needed.
//----------------------------------------------------------------------------------
double diagonalise( mat& H, vec& evals, mat& evects, vec& occupation)
{
  int gubbins;
  return diagonalise(H,evals,evects,occupation,gubbins);
}

#endif //MODULEC_HPP

