//-----------------------------------------------------------------------------
// Module D
//
// Description - Contains functions used to calculate the repulsive energy
//-----------------------------------------------------------------------------
#ifndef MODULEDNEW_HPP
#define MODULEDNEW_HPP

#include "../atom/atom.h"
#include "armadillo"
#include <cmath>
#include "../atom/vect3.h"
#include "../sim/unitcell.h"
#include "../sim/simulation.hpp"

using namespace moldycat;
using namespace arma;

//Define parameters as taken from Xu 1992 paper
const double phi0=8.18555;
const double d0=1.64;
const double d1=2.57;
const double m=3.30304;
const double mc=8.6655;
const double dc=2.1052;
const double c0r=2.2504290109E-8;
const double c1r=-1.4408640561E-6;
const double c2r=2.1043303374E-5;
const double c3r=6.6024390226E-5;
const double f0=-2.5909765118191;
const double f1=0.5721151498619;
const double f2=-1.7896349903996E-3;
const double f3=2.3539221516757E-5;
const double f4=-1.24251169551587E-7;

//-------------------------------------------------------------
// Function phi
//-------------------------------------------------------------
// Returns: phi(r) (from Xu et al.)
//
// Params: r - Distance between atoms 
//
// Operation: Calculates phi(r) used to calculate the repulsive
//            energy and the repulsive force
//-------------------------------------------------------------
double phi(double r) { 

  if (r>d1) {
    double r_2=r-d1;
    return c0r + r_2*(c1r + r_2*(c2r + r_2*(c3r)));
  }
  
  else {
    return phi0*pow((d0/r),m)*exp(m*(pow(d0/dc,mc)-pow(r/dc,mc)));
  }

}

//------------------------------------------------------------------------------
// Function phiDeriv
//------------------------------------------------------------------------------
// Returns: Phi'(r)
//
// Params: r - Distance between atoms 
//
// Operation: Calculates the derivative of Phi(r) at r, used to
//            calcaulte the repulsive forces
//------------------------------------------------------------------------------
double phiDeriv(double r) {
  
  if (r>d1) {
    double r_2=r-d1;
    return r_2*c1r*1.0+r_2*(c2r*2.0+r_2*(c3r*3.0));
  }
  
  else {
    return -m*phi(r)*((1.0/r)+(mc/r)*pow((r/dc),mc)); 
  }

}

//------------------------------------------------------------------------------
// Function repulsiveEnergy
//------------------------------------------------------------------------------
// Returns: E_rep - Total repulsive energy
//
// Params: atoms(atom_list) - Positions of atoms 
//         cell(unitcell)   - Unit cell
//
// Operation: works out displacement betweens atoms and returns total repulsive 
//            energy accordingly.
//------------------------------------------------------------------------------
double repulsiveEnergy(const atom_list& atoms, unitcell cell) {

  int n=atoms.count(),i,j;
  double x,R,E_rep=0.0,dif;
  vec E=zeros<vec>(n);
  vect3 r;

  for (i=0;i<n;i++) {
    x=0.0;
    
    for (j=0;j<n;j++) {
      
      if(j==i) { 
        continue;
      }
      
      else {
        r=atoms[j].pos.diff_bc(atoms[i].pos,cell);
        R=r.norm();
        if (R<r_m) {
          x+=phi(R);
        }
      }
    
    }
    
    dif=f0+x*(f1+x*(f2+x*(f3+x*f4)));
    E(i)=dif;
    E_rep+=dif;
  
  }
  return E_rep;
}

#endif
