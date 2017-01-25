//-----------------------------------------------------------------------------
// Module B
//
// Description - Contains all the "matrixy" functions and sub-routines:
//               scalingFunction
//
//               --> Used in the Xu et al. model to scale the elements
//
//               scalingfunctionDeriv
//               --> Derivative of the scalingFunction
//
//               buildHamiltonian
//               --> Builds the TB Hamiltonian matrix
//
//               buildHamiltonianDeriv
//               --> Builds the derivative of the TB Hamiltonian, one matrix
//                   for each component
//
//               buildDensityMatrix
//               --> Builds the density matrix (parallelised)
//
//               buildDensityMatrixSlowly
//               --> Builds the density matrix (unparallelised)
//-----------------------------------------------------------------------------
#ifndef MODULEB_HPP
#define MODULEB_HPP

#include "../atom/atom.h"
#include "armadillo"
#include "../atom/vect3.h"
#include "../sim/unitcell.h"
#include "omp.h"

using namespace moldycat;
using namespace arma;

// Actual constants from Xu 1992 paper.
const double r_m = 2.6;
const double r_0 = 1.536329;
const double r_c = 2.18;
const double r_1 = 2.45;
const double V_ss = -5.0;
const double V_sp = 4.7;
const double V_pp_sigma = 5.5;
const double V_pp_pi = -1.55;
const double n = 2.0;
const double n_c = 6.5;
const double E_s = -2.99;
const double E_p = 3.71;

//Polynomial coefficients, also from Xu 1992.
const double c0 = 6.7392620074314E-3;
const double c1 = -8.1885359517898E-2;
const double c2 = 0.1932365259144;
const double c3 = 0.3542874332380;

//-----------------------------------------------------------------------------
// Function scalingFunction
//
// Returns: Scaling factor
//
// Params: r - Distance between atoms, unit cell 
//
// Operation - Gives the scaling factor for calculation of Hamiltonian element
//-----------------------------------------------------------------------------
double scalingFunction(double r) 
{
  //Condition to use tail function
  if (r > r_1) {
    double r_2 = r - r_1;
    return c0 + r_2*(c1 + r_2*(c2 + r_2*(c3))); 
  }

  // Else use normal scaling function
  else {
    return pow(r_0/r,n)*exp(n*(-pow(r/r_c,n_c) + pow(r_0/r_c,n_c)));
  }
} 

//-----------------------------------------------------------------------------
// Function scalingFunctionDeriv
//
// Returns: Derivative of scaling factor at r
//
// Params: r - Distance between atoms, unit cell 
//
// Operation - Gives the scaling factor derivative for calculation of 
//             Hamiltonian derivative
//-----------------------------------------------------------------------------
double scalingFunctionDeriv(double r) 
{
  // Use tail function derivative if r > r_1
  if (r > r_1) {
    double r_2 = r - r_1;
    return c1 + 2*c2*r_2 + 3*c3*r_2*r_2; 
  }

  // Else use normal scaling function derivative
  else {
    return -(n/r)*(1+n_c*pow(r/r_c,n_c))*scalingFunction(r); 
  }
} 

//-----------------------------------------------------------------------------
// Function buildHamiltonian
//
// Returns: None
//
// Params: atoms (atom_list) - List of atoms
//         unitcell (cell)   - Dimensions of unit cell 
//         H (matrix)        - The Hamiltonian
//
// Operation: Builds the Hamiltonian.
//
// Notes - Can cut cost by half using symmetry.
//-----------------------------------------------------------------------------
void buildHamiltonian(const atom_list& atoms, unitcell& cell, mat& H)
{   
  const int n=atoms.count();

  //V_pp_delta used to simplify calculations/code (see below).
  const double V_pp_delta = V_pp_sigma - V_pp_pi;

  //Loop over all atoms (blocks of 4 rows in the Hamiltonian)
  for (int i=0;i<n;i++) {

    // direction_array is an array of 4 elements containing '1' and the 3
    // direction cosines. 
    double direction_array[4];
    direction_array[0] = 1;  

    //Loop over all atoms (blocks of 4 columns in the Hamiltonian)
    for (int j=i;j<n;j++) {
      
      //For 4x4 blocks down the main diagonal. Only the diagonal elements need
      // to be filled.
      if (i==j) {
        H(4*i, 4*j) = E_s;
        H(4*i+1, 4*j+1) = E_p;
        H(4*i+2, 4*j+2) = E_p;
        H(4*i+3, 4*j+3) = E_p;
        continue;
      }
      
      //Unit vector from atom i to atom j (=the direction cosines)
      vect3 displacement = atoms[j].pos.diff_bc(atoms[i].pos,cell);
      
      //Distance from atom i to atom j
      double r = displacement.norm();

      //We are now reusing a normalised displacement where we had dc before
      displacement/=r;
    
      //Scale factor for distance r
      double sf = scalingFunction(r);
          
      //Loop over all orbitals (single rows in each block)
      for (int a=0;a<4;a++) {
       
        //Loop over all orbitals (single columns in each block)
        for (int b=a;b<4;b++) {
          int x = 4*i + a;
          int y = 4*j + b;
          
          //Diagonally opposite elements within the same 4x4 block
          // used later to save computational time.
          int x_opp = 4*i + b;
          int y_opp = 4*j + a;
          
          //Fill in array
          if (r < r_m) {

            //Fill in direction_array (0,l,m,n)
            for (int k=0;k<3;k++) {
              direction_array[k+1] = displacement(k);
            }
            
            //-----------------------------------------------------------------
            // Diagram of what general 4x4 block (i!=j) looks like 
            // Note: Each elements needs to be multiplied by the scale factor.
            //       DF = direction factor --> see below.
            //       V_delta = V_pp_sigma - V_pp_pi
            //       
            //       The 4x4 block is almost symmetric except for the 
            //       first row and first column
            //-----------------------------------------------------------------
 //==========================================================================
 //       ||                                a                               |
 //       ||-----------------------------------------------------------------
 //       ||       0        |       1       |       2       |        3      | 
 //===|===||=================================================================
 //   |   ||                |               |               |               |
 //   | 0 ||    DF*V_ss     | DF*V_sp_sigma | DF*V_sp_sigma | DF*V_sp_sigma |
 //   |   ||                |               |               |               |
 //   |---||-----------------------------------------------------------------
 //   |   ||                |  DF*V_delta   |               |               |
 //   | 1 || -DF*V_sp_sigma |   + V_pp_pi   |   DF*V_delta  |  DF*V_delta   |
 //   |   ||                |               |               |               |
 // b |---||-----------------------------------------------------------------
 //   |   ||                |               |   DF*V_delta  |               |
 //   | 2 || -DF*V_sp_sigma |  DF*V_delta   |   + V_pp_pi   |  DF*V_delta   |
 //   |   ||                |               |               |               |
 //   |---||-----------------------------------------------------------------
 //   |   ||                |               |               |  DF*V_delta   |
 //   | 3 || -DF*V_sp_sigma |  DF*V_delta   |   DF*V_delta  |   + V_pp_pi   |
 //   |   ||                |               |               |               |
 //---|---||-----------------------------------------------------------------

            // direction_factor is used to simplify the code:
            // This is the table of direction factors (DF)
            //------------------------------------
            //\     ||      |      |      |      |
            //  \ a || 0(1) | 1(l) | 2(m) | 3(n) |
            //  b \ ||      |      |      |      |
            //====================================
            //      ||      |      |      |      |
            // 0(1) || 1*1  | l*1  | m*1  | n*1  |
            //      ||      |      |      |      |
            //------------------------------------
            //      ||      |      |      |      |
            // 1(l) || 1*l  | l*l  | m*l  | n*l  |
            //      ||      |      |      |      |
            //------------------------------------
            //      ||      |      |      |      |
            // 2(m) || 1*m  | l*m  | m*m  | n*m  |
            //      ||      |      |      |      |
            //------------------------------------
            //      ||      |      |      |      |
            // 3(n) || 1*n  | l*n  | m*n  | n*n  |
            //      ||      |      |      |      |
            //------------------------------------
            double direction_factor = direction_array[a]*direction_array[b];

            if ((a == 0) && (b == 0)) {
              H(x,y) = sf*V_ss; 
            }
            
            else if (a == b) {
              H(x,y) = sf*(direction_factor*V_pp_delta + V_pp_pi);
            }
            
            else if (a == 0) {
              H(x,y) = sf*direction_factor*V_sp; 
              //Use anti-symmetry of first row with first column
              H(x_opp,y_opp) = -H(x,y);
            }
            
            else {
              H(x,y) = sf*direction_factor*V_pp_delta; 
              H(x_opp,y_opp) = H(x,y);
            }
          }
        }
      }
    } 
  }
  //Make fill in lower half of matrix via symmetry
  H=symmatu(H);
  return;
}

//-----------------------------------------------------------------------------
// Function buildHamiltonianDeriv
//
// Returns: None
//
// Params: atoms(atom_list) - List of atoms 
//         cell(unitcell) - Unit cell
//         dH(cube) - 4*n x 4*n x 3 matrix containing the elements of the 
//                    derivative in either the x,y or z directions - each 
//                    4*n x 4*n slice contains the derivative elements in a
//                    specific elements. 
//                    e.g. dH(x,y,1) is the derivative of the (x,y) element of
//                         H resolved in the y direction.
//         
// Operation: Builds the Hamiltonian derivative cube, dH, composed of 3
//            matrices one representing the derivatives in each component.
//            The i^{th} block of 4 rows in matrix c contains the derivatives 
//            of the corresponding block of the Hamiltonian with respect to the 
//            c^{th} component of r_i. 
//-----------------------------------------------------------------------------
void buildHamiltonianDeriv(const atom_list& atoms,const unitcell& cell, cube& dH)
{   
  const int n=atoms.count();
  const double V_pp_delta = V_pp_sigma - V_pp_pi;

#pragma omp parallel for
  //Loop over all atoms (blocks of 4 rows in the Hamiltonian)
  for (int i=0;i<n;i++) {
    
    //Loop over all atoms (blocks of 4 columns in the Hamiltonian)
    for (int j=i+1;j<n;j++) {
      
      //Loop over all orbitals (single rows in each block)
      for (int a=0;a<4;a++) {
       
        //Loop over all orbitals (single columns in each block)
        for (int b=0;b<4;b++) {
          int x = 4*i + a;
          int y = 4*j + b;
		      
          //Vector from atom i to atom j
          vect3 displacement = atoms[j].pos.diff_bc(atoms[i].pos,cell);
         
          //Direction cosines (l,m,n)
          vect3 dc=displacement.getUnit();

          //Distance from atom i to atom j
          double r = displacement.norm();

          //Derivative of the scale factor
          double sf_deriv = scalingFunctionDeriv(r);
          
          //Scale factor
          double sf = scalingFunction(r);
          
          if (r < r_m) {
            
            //Loop over each component
            for (int comp=0;comp<3;comp++) {
              
              //====================== For each 4x4x3 block=====================
              
              //---------------Covers (0,0,0), (0,0,1) and (0,0,2)--------------
              if ((a == 0) && (b == 0)) {
                dH(x,y,comp) = sf_deriv*V_ss*dc(comp); 
              }
              //----------------------------------------------------------------

              //---------Other elements in the diagonal of the 4x4 block--------
              else if (a == b) {
                
                //Covers (1,1,0), (2,2,1) & (3,3,2)------------------------
                if ((a-1)==comp) {
                  dH(x,y,comp) = ( 2*dc(comp)*(1-dc(comp)*dc(comp))/r*sf \
                                   + dc(comp)*dc(comp)*dc(comp)*sf_deriv \
                                 ) *V_pp_delta \
                                 + sf_deriv*dc(comp)*V_pp_pi;
                }                                                           
                //---------------------------------------------------------
                
                //Covers (1,1,1),(1,1,2),(2,2,0),(2,2,2),(3,3,0),(3,3,1)
                else {
                  dH(x,y,comp) = ( -2*dc(a-1)*dc(a-1)*dc(comp)/r*sf \
                                   + dc(a-1)*dc(a-1)*dc(comp)*sf_deriv \
                                 ) *V_pp_delta \
                                 + sf_deriv*dc(comp)*V_pp_pi;
                }
                //---------------------------------------------------------
              }
              //----------------------------------------------------------------

              //---------Other elements in the first row of the 4x4 block-------
              else if (a == 0) {

                //Covers (0,1,0),(0,2,1),(0,3,2)---------------------------
                if ((b-1)==comp) {
                  dH(x,y,comp) = ( (1-dc(comp)*dc(comp))/r*sf \
                                   + dc(comp)*dc(comp)*sf_deriv \
                                 ) *V_sp;
                }
                //---------------------------------------------------------

                //Covers (0,1,1),(0,1,2),(0,2,0),(0,2,2),(0,3,0),(0,3,1)---
                else {
                  dH(x,y,comp) = ( -dc(b-1)*dc(comp)/r*sf \
                                   + dc(b-1)*dc(comp)*sf_deriv \
                                 ) *V_sp;
                }
                //---------------------------------------------------------
              }
              
              //-------Other elements in the first column of the 4x4 block------
              else if (b == 0) {

                //Covers (1,0,0),(2,0,1),(3,0,2)---------------------------
                if ((a-1)==comp) {
                  dH(x,y,comp) = -( (1-dc(comp)*dc(comp))/r*sf \
                                    + dc(comp)*dc(comp)*sf_deriv \
                                  ) *V_sp;
                }
                //---------------------------------------------------------
                
                //Covers (1,0,1),(1,0,2),(2,0,0),(2,0,2),(3,0,0),(3,0,1)---
                else {
                  dH(x,y,comp) = -( -dc(a-1)*dc(comp)/r*sf \
                                    + dc(a-1)*dc(comp)*sf_deriv \
                                  ) *V_sp;
                }
                //---------------------------------------------------------
              }
              //----------------------------------------------------------------
              

              //----------------------All other elements------------------------
              else {
                
                //Covers (1,2,0),(1,3,0),(2,3,1),(2,1,1),(3,1,2),(3,2,2)---
                if ((a-1)==comp) {
                  dH(x,y,comp) = ( (1-2*dc(a-1)*dc(a-1))*dc(b-1)/r*sf \
                                   +dc(a-1)*dc(b-1)*dc(comp)*sf_deriv \
                                 ) *V_pp_delta;
                }
                //---------------------------------------------------------

                //Covers (1,2,1),(1,3,2),(2,1,0),(2,3,2),(3,1,0),(3,2,1)---
                else if ((b-1)==comp) {
                  dH(x,y,comp) = ( (1-2*dc(b-1)*dc(b-1))*dc(a-1)/r*sf \
                                   + dc(a-1)*dc(b-1)*dc(comp)*sf_deriv \
                                 ) *V_pp_delta;
                }
                //---------------------------------------------------------

                //Covers (1,2,2),(1,3,1),(2,1,2),(2,3,0),(3,1,1),(3,2,0)---
                else {
                  dH(x,y,comp) = ( -2*dc(a-1)*dc(b-1)*dc(comp)/r*sf \
                                   +dc(a-1)*dc(b-1)*dc(comp)*sf_deriv \
                                 ) *V_pp_delta;
                }
                //---------------------------------------------------------
              }
              //----------------------------------------------------------------

              //To make dH antisymmetric
              dH(y,x,comp) = -dH(x,y,comp);
            }
          }
        }
      } 
    }
  }
  return;
}


//-----------------------------------------------------------------------------
// Function buildDensityMatrix
//
// Returns: None
//
// Params: atoms(atom_list) - List of atoms 
//         cell(unitcell)   - Unit cell
//         evals(vec)       - Vector of eigenvalues (ascending) 
//         evects(mat)      - Matrix of corresponding eigenvectors
//         occupation(vec)  - Vector of occupation numbers
//         max_occupied_index(int) - index of Fermi level within energy spectrum
//          (typically allows computational effort to be halved as top N/2 
//              states unfilled)
//         
// Operation: Builds the density matrix used to calculate forces 
//-----------------------------------------------------------------------------
void buildDensityMatrix(const atom_list& atoms, unitcell& cell, vec& evals, mat& evects, mat& DM, const vec& occupation, int max_occupied_index) 
{
  int x,y,n=atoms.count();
  double elem;

  //First do elements in diagonal 4x4 blocks to prevent double counting 
  //Loop over all atoms (blocks of 4 rows/columns in the Hamiltonian)
  #pragma omp parallel for private(x,y,elem)
  for (int i=0;i<n;i++) {

    //Loop over all orbitals (single rows in each 4x4 block)
    for (int a=0;a<4;a++) {       
      //Loop over all orbitals (single columns in each 4x4 block)
      for (int b=0;b<4;b++) {
        x = 4*i + a;
        y = 4*i + b;

        //Add element to density matrix
        for (int level=0;level<max_occupied_index;level++) {
          DM(x,y) += occupation(level)*evects(x,level)*evects(y,level);
        }
      }
    }

    //Then do the elements in the upper off-diagonal 4x4 blocks 
    // and use use symmetry (rho_ij = rho_ji)
    
    //Loop over all atoms (blocks of 4 columns in the Hamiltonian)
    for (int j=i+1;j<n;j++) {  

      //Loop over all orbitals (single rows in each 4x4 block)
      for (int a=0;a<4;a++) {
       
        //Loop over all orbitals (single columns in each 4x4 block)
        for (int b=0;b<4;b++) {
          x = 4*i + a;
          y = 4*j + b;

          //Add element to density matrix
          for (int level=0;level<max_occupied_index;level++) {
            elem =  occupation(level)*evects(x,level)*evects(y,level);
            DM(x,y) += elem;
            DM(y,x) += elem;
          }
        }
      }
    } 
  }
}

//-----------------------------------------------------------------------------
// Function buildDensityMatrixSlowly
//
// Returns: None
//
// Params: atoms(atom_list) - List of atoms 
//         cell(unitcell)   - Unit cell
//         evals(vec)       - Vector of eigenvalues (ascending) 
//         evects(mat)      - Matrix of corresponding eigenvectors
//         occupation(vec)  - Vector of occupation numbers
//         max_occupied_index(int) - index of Fermi level within energy spectrum
//          (typically allows computational effort to be halved as top N/2 
//              states unfilled)
//         
// Operation: Builds the density matrix used to calculate forces 
//            Same as buildDensity (above) but used in phonon code where 
//            parallelisation occurs at top level and too many threads would 
//            slow down computation
//-----------------------------------------------------------------------------
void buildDensityMatrixSlowly(const atom_list& atoms, unitcell& cell, vec& evals, mat& evects, mat& DM, const vec& occupation, int max_occupied_index) 
{
  int x,y,n=atoms.count();
  double elem;
  
  //First do elements in diagonal 4x4 blocks to prevent double counting 
  //Loop over all atoms (blocks of 4 rows/columns in the Hamiltonian)
  for (int i=0;i<n;i++) {

    //Loop over all orbitals (single rows in each 4x4 block)
    for (int a=0;a<4;a++) {       

      //Loop over all orbitals (single columns in each 4x4 block)
      for (int b=0;b<4;b++) {
        x = 4*i + a;
        y = 4*i + b;

        //Add element to density matrix
        for (int level=0;level<max_occupied_index;level++) {
          DM(x,y) += occupation(level)*evects(x,level)*evects(y,level);
        }
      }
    }

    //Then do the elements in the upper off-diagonal 4x4 blocks 
    // and use use symmetry (rho_ij = rho_ji)

    //Loop over all atoms (blocks of 4 columns in the Hamiltonian)
    for (int j=i+1;j<n;j++) {  
      //Loop over all orbitals (single rows in each 4x4 block)
      for (int a=0;a<4;a++) {
       
        //Loop over all orbitals (single columns in each 4x4 block)
        for (int b=0;b<4;b++) {
          x = 4*i + a;
          y = 4*j + b;

          //Add element to density matrix
          for (int level=0;level<max_occupied_index;level++) {
            elem =  occupation(level)*evects(x,level)*evects(y,level);
            DM(x,y) += elem;
            DM(y,x) += elem;
          }
        }
      }
    } 
  }
}
#endif //MODULEB_HPP
