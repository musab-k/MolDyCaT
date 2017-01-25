#ifndef EVSTRETCHER_HPP 
#define EVSTRETCHER_HPP

#include <cmath>
#include <vector>
#include <fstream>
#include "armadillo"
#include <iostream>
#include "../atom/vect3.h"
#include "../atom/atom.h"
#include "../sim/unitcell.h"
#include "../sim/simulation.hpp"

//-----------------------------------------------------------------------------------
// Function evStretcher
//
// Params:      - atom list, unit cell, scaling factor, filename (optional) 
// Returns:     - nothing
// Operation:   - Works out the energy/atom-nearest neighbour distance curve by scaling
//              the cell (and supersizing it if required).
//              - Writes to outputs/energy-volume by default, or filename specified.
//-----------------------------------------------------------------------------------
void evStretcher(atom_list& atoms_init, unitcell& cell_init, int scaleDir[], string filename="outputs/energy-volume"){
    int n=atoms_init.count(), i, j, nsuper,sx,sy,sz;
    double Lx, Ly, Lz, x, xMax, dx;
    double Lx_init, Ly_init, Lz_init;
    double tb_energy, repulsive_energy;
    mat H, evects;
    vec evals, occupation;
    ofstream evCurves;
    evCurves.open(filename.c_str());

    double nearestNeighbour;

    nearestNeighbour = atoms_init.nearestNeighbour(cell_init);

    // Set up initial cell of n atoms.
    atom_list atoms_frac(n);
    atom_list atoms_current(n);
    Lx_init=cell_init.len.x();
    Ly_init=cell_init.len.y();
    Lz_init=cell_init.len.z();
    unitcell cell=cell_init;
    sx = scaleDir[0];
    sy = scaleDir[1];
    sz = scaleDir[2];

    // Obtain fractional co-ordinates.
    for ( i=0 ; i<n ; i++ ){
        atoms_frac[i].pos.x() = atoms_init[i].pos.x()/Lx_init;
        atoms_frac[i].pos.y() = atoms_init[i].pos.y()/Ly_init;
        atoms_frac[i].pos.z() = atoms_init[i].pos.z()/Lz_init;
    }

    x=-0.15;
    xMax=0.15;
    dx=(xMax-x)/15.0;
    for ( j=0 ; x<xMax ; j++ ){
        // Set the new cell size
        Lx=Lx_init*(1.0+x);
        Ly=Ly_init*(1.0+x);
        Lz=Lz_init*(1.0+x);
        cell.len=vect3(Lx,Ly,Lz);
      //  cout << cell.volume() << endl;
        // Calculate new co-ordinates by scaling fractional coordinates.
        for ( i=0 ; i<n ; i++ ){
            atoms_current[i].pos.x() = atoms_frac[i].pos.x()*Lx;
            atoms_current[i].pos.y() = atoms_frac[i].pos.y()*Ly;
            atoms_current[i].pos.z() = atoms_frac[i].pos.z()*Lz;
        }

        // Super-size the cell 
        atoms_current.scale_supercell(sx,sy,sz,cell);
        nsuper=atoms_current.count();

        cout << nsuper << endl;
        // Work out the energies
        H = zeros<mat>(4*nsuper,4*nsuper);
        evects = zeros<mat>(4*nsuper,4*nsuper);
        evals = zeros<vec>(4*nsuper);
        occupation = zeros<vec>(4*nsuper);

        buildHamiltonian(atoms_current, cell, H);                                                                 
        tb_energy = diagonalise(H,evals,evects,occupation);
        repulsive_energy = repulsiveEnergy(atoms_current, cell);
        
        // Write to file
        evCurves << nearestNeighbour*(1.0+x) << "\t" << (tb_energy+repulsive_energy)/(nsuper) << endl;

        // Re-initialise atoms. 
        atoms_current = atoms_init;
        x+=dx;
    }
    cout << nsuper << endl;
    evCurves.close();
}
#endif
