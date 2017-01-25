//--------------------------------------------------------------------
// Module G 
//
// Description: Contains functions to write positions, velocities
//              etc. to output files.
//--------------------------------------------------------------------
#ifndef MODULEG_HPP
#define MODULEG_HPP

#include "../atom/atom.h"
#include "armadillo"
#include "../atom/vect3.h"
#include "../sim/unitcell.h"
#include <fstream>
#include "string.h"
#include <cmath>

using namespace moldycat;
using namespace arma;

//--------------------------------------------------------------------
// Function writePositions
//--------------------------------------------------------------------
// Returns: nothing
//
// Params: atom_list, append
//
// Operation: writes positions of atoms to file for visualisation 
// using jMOL/VMD. If append=1, it appends the positions, otherwise it 
// creates a new file. If called after every time step jMOL/VMD can 
// animate the positions.
//--------------------------------------------------------------------
void writePositions(const atom_list& atoms, int append, string filename = "outputs/positions.xyz"){
  ofstream posfile;
  
  if (append==1)
    posfile.open(filename.c_str(),ios::out | ios::app);
  else
    posfile.open(filename.c_str());
  
  int n=atoms.count();
  posfile << n << endl << endl;
  posfile.precision(10);
  
  for (int i=0; i<n; i++){
    posfile << "C \t" << atoms[i].pos.x() << "\t" << atoms[i].pos.y();
    posfile << "\t" << atoms[i].pos.z() << endl;
  }

  posfile.close();
} 

//--------------------------------------------------------------------
// Function writeVelocities
//--------------------------------------------------------------------
// Returns: nothing
//
// Params: atom_list, append, filename
//
// Operation: Writes velocities to output file
//--------------------------------------------------------------------
void writeVelocities(const atom_list& atoms, int append, string filename = "outputs/velocities.out"){
  ofstream posfile;

  if (append==1)
    posfile.open(filename.c_str(),ios::out | ios::app);
  else
    posfile.open(filename.c_str());
    
  int n=atoms.count();
  posfile << n << endl << endl;
  posfile.precision(10);
    
  for (int i=0; i<n; i++){
    posfile << "C \t" << atoms[i].vel.x() << "\t" << atoms[i].vel.y();
    posfile << "\t" << atoms[i].vel.z() << endl;
  }

  posfile.close();
}

//--------------------------------------------------------------------
// Function writeVelocities
//--------------------------------------------------------------------
// Returns: nothing
//
// Params: atom_list, append, filename
//
// Operation: Writes velocities to output file
//--------------------------------------------------------------------
void writeUnitcell(const simulation& sim, string filename = "outputs/unitcell.out"){
  ofstream ufile;
  ufile.open(filename.c_str());
  vect3 length = (sim.cell).len;
  vect3 angles = (sim.cell).ang;
  const bctype* boundary = (sim.cell).bc;
  
  ufile.precision(10);
  ufile << "LENGTH_X_1 " << length.x() << endl;
  ufile << "LENGTH_X_2 " << length.y() << endl;
  ufile << "LENGTH_X_3 " << length.z() << endl;
  ufile << "ALPHA " << angles[0]*180.0/M_PI << endl;
  ufile << "BETA " << angles[1]*180.0/M_PI << endl;
  ufile << "GAMMA " << angles[2]*180.0/M_PI << endl;
  ufile << "BOUND_X_1 " << boundary[0] << endl;
  ufile << "BOUND_X_2 " << boundary[1] << endl;
  ufile << "BOUND_X_3 " << boundary[2] << endl;

}

//--------------------------------------------------------------------
// Function writePositionsandCell
//--------------------------------------------------------------------
// Returns: nothing
//
// Params: atom_list, append, filename
//
// Operation: Meow 
//--------------------------------------------------------------------
void writePostionsandCell(const atom_list& atoms, const unitcell& cell, \
                          vect3 Lx, vect3 Ly, vect3 Lz, \
                          string filename = "outputs/positions.xyz") {
  
  ofstream posfile;
  posfile.open(filename.c_str(),ios::out | ios::app);
  int n=atoms.count();
  atom_list newatoms(n);
  int i,j,k;
  posfile << n+60 << endl << endl;
  posfile.precision(10);

  for (i=0; i<n; i++){
    posfile << "C \t" << atoms[i].pos.x() << "\t" << atoms[i].pos.y();
    posfile << "\t" << atoms[i].pos.z() << endl;
  }

  vect3 x;
  mat rot;

  vect3 A,B,C;

  double ca,sa;

  ca = cos(cell.ang.x());
  sa = sin(cell.ang.x());

  A = vect3(ca*cell.len.x()-sa*cell.len.y(),sa*cell.len.x()+ca*cell.len.y(),0.0);
  B = vect3(ca*cell.len.x()-sa*cell.len.y(),sa*cell.len.x()+ca*cell.len.y(),0.0);
  C = vect3(0.0,0.0,cell.len.z());

  A = Lx;
  B = Ly;
  C = Lz;

  k=0;
  
  for (j=0; j<19; j+=2){
    newatoms[j].pos = vect3(A*0.1*k);
    newatoms[j+1].pos = vect3(B+A*0.1*k);
    k+=1;
  }

  k=0;

  for (j=20; j<39; j+=2){
    newatoms[j].pos = vect3(B*0.1*k);
    newatoms[j+1].pos = vect3(A+B*0.1*k);
    k+=1;
  }
  
  k=0;

  for (j=40; j<59; j+=2){
    newatoms[j].pos = vect3(C*0.1*k);
    newatoms[j+1].pos = vect3(A+C*0.1*k);
    k+=1;
  }

  for (i=0; i<n; i++){
    posfile << "S \t" << newatoms[i].pos.x() << "\t" << newatoms[i].pos.y();
    posfile << "\t" << newatoms[i].pos.z() << endl;
  }
    
  posfile.close();
} 

#endif
