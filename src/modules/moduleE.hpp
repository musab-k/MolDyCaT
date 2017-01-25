//-----------------------------------------------------------------------------
// Module E
//
// Description - Contains functions used to calculate the Hellmann-Feynman and
//               repulsive forces (separately).
//-----------------------------------------------------------------------------
#ifndef MODULEE_HPP
#define MODULEE_HPP
#define TINY 1e-14

#include "../atom/atom.h"
#include "armadillo"
#include <cmath>
#include "../atom/vect3.h"
#include "../sim/unitcell.h"
#include "omp.h"

using namespace moldycat;
using namespace arma;

//-----------------------------------------------------------------------------
// Function TBforce 
//
// Returns: vector of force on a-th orbital of atom i due to b-th orbital of
//          atom j
//
// Params: i - Atom number
//         j - Atom number
//         a - Orbital number 
//         b - Orbtial number
//        dH - Deriative of Hamiltonian (cube)
//        DM - Density matrix
//
// Operation - Calculates Hellman-Feynman force on a-th orbital of atom i due 
//             to b-th orbital of atom j
//-----------------------------------------------------------------------------
vect3 TBforce(int i, int j, int a, int b,const cube& dH,const mat& DM,const vect3 displacement) 
{
  vect3 force(0,0,0);
  int x = 4*i + a;
  int y = 4*j + b;
  
  force.x() = -2*DM(x,y)*dH(y,x,0);
  if (abs(force.x()) < TINY) {
    force.x() = 0.0;
  }
  
  force.y() = -2*DM(x,y)*dH(y,x,1);
  if (abs(force.y()) < TINY) {
    force.y() = 0.0;
  }

  force.z() = -2*DM(x,y)*dH(y,x,2);
  if (abs(force.z()) < TINY) {
    force.z() = 0.0;
  }

  return force;
}

//-----------------------------------------------------------------------------
// Function buildTBforces 
//
// Returns: None 
//
// Params (type): atoms(atom_list) - Positions of atoms
//                cell(unitcell)   - Unit cell
//                dH(cube)         - The "Hamiltonian derivative" cube
//                DM(mat)          - Density matrix
//                forces(vec_list) - Vector list of Hellmann-Feynamnn forces
//
// Operation - Constructs vectors list of the Hellmann-Feynman forces on each
//             atom.
//-----------------------------------------------------------------------------
void buildTBforces(const atom_list& atoms, unitcell& cell,cube& dH, mat& DM, vec_list& forces)
{
  int n=atoms.count(),i,j,a,b;
  vect3 displacement;
  vec_list test(n);
  vect3 f;
  
  //Loop over atoms that the force is to be calculated on
  for (i=0; i<n; i++) {

    //Loop over atoms that the force is due to
    for (j=0;j<n;j++) {
      
      // No force contribution from own atom
      if (j==i) { 
        continue; 
      }
      
      else {
        displacement = atoms[j].pos.diff_bc(atoms[i].pos,cell);

        //Optional if statement (since Hamiltonian element should be 0 if r > r_m anyway)
        if(displacement.norm()<=r_m) {
          
          //Loop over orbitals on atom i
          for(a=0;a<4;a++){
            
            //Loop over orbtials on atom j
            for(b=0;b<4;b++){
              f=TBforce(i,j,a,b,dH,DM,displacement);
              (forces[i]) += f;
            }
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function buildRepulsiveForces
//
// Returns: None 
//
// Params (type): atoms(atom_list)   - Positions of atoms
//                cell(unitcell)     - Unit cell
//                repForces(veclist) - List of repulsive forces on each atom
//
// Operation - Constructs list of repulsive forces. 
//-----------------------------------------------------------------------------
void buildRepulsiveForces(const atom_list& atoms, const unitcell& cell, vec_list& repForces ) {
  int n=atoms.count(),i,j; 
  double r;
  double x;
  vect3 y;
  vect3 f;
  vect3 force,displacement;
  
  for(i=0;i<n;i++){
    x=0.0;
    y=vect3(0.0,0.0,0.0);
    f=vect3(0.0,0.0,0.0);
    
    for(j=0;j<n;j++){
      
      if (i==j) { 
        continue; 
      }
      
      else {
        displacement=atoms[j].pos.diff_bc(atoms[i].pos, cell);
        r=displacement.norm();
        
        if (r<=r_m) {
          displacement.normalise();
          x+=phi(r);
          y+=phiDeriv(r)*displacement;
        }
      }
    }
   
   f.x()+=y.x()*(f1+(x*(2*f2+x*(3*f3+x*(4*f4)))));
   if (abs(f.x())<TINY) {
     f.x()=0;
   }

   f.y()+=y.y()*(f1+(x*(2*f2+x*(3*f3+x*(4*f4)))));
   if (abs(f.y())<TINY) {
     f.y()=0;
   }   
   
   f.z()+=y.z()*(f1+(x*(2*f2+x*(3*f3+x*(4*f4)))));
   if (abs(f.x())<TINY) {
     f.x()=0;
   }   
   
   repForces[i]+=f;
  }
}

#endif //MODULEE_HPP
