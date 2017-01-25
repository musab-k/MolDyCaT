#ifndef ELASTICCONSTANTS_HPP
#define ELASTICCONSTANTS_HPP

#include <cmath>
#include <vector>
#include <fstream>
#include "armadillo"
#include <iostream>
#include "../atom/vect3.h"
#include "../atom/atom.h"
#include "../sim/unitcell.h"
#include "../sim/simulation.hpp"

//------------------------------------------------------------------------
// Function elasticConstant 
//
// Params: atomlist, unit cell, filename (optional)
// Operation: Works out the elastic constants by supersizing the unit
// cell (by an amount specified in scaling factor) and writes them to 
// file.
//------------------------------------------------------------------------
void elasticConstants(atom_list& atoms_init, unitcell& cell_init, string filename="outputs/elasticConst"){
    int n=atoms_init.count(), i, j;
    double Lx, Ly, Lz, x, xMax, dx;
    double Lx_init, Ly_init, Lz_init;
    double tb_energy, repulsive_energy;
    mat H, evects;
    vec evals, occupation;
    ofstream c11c12;
    string c1filename = filename + "-c11c12.txt";
    c11c12.open(c1filename.c_str());

    // Set up initial cell of n atoms.
    atom_list atoms_frac(n);
    atom_list atoms_current(n);
    Lx_init=cell_init.len.x();
    Ly_init=cell_init.len.y();
    Lz_init=cell_init.len.z();
    unitcell cell=cell_init;

    // Obtain fractional co-ordinates.
    for ( i=0 ; i<n ; i++ ){
        atoms_frac[i].pos.x() = atoms_init[i].pos.x()/Lx_init;
        atoms_frac[i].pos.y() = atoms_init[i].pos.y()/Ly_init;
        atoms_frac[i].pos.z() = atoms_init[i].pos.z()/Lz_init;
    }

    x=-0.05;
    xMax=0.05;
    dx=(xMax-x)/25.0;
        //writePositions(atoms_current, 0); 
    // Work out C11-C12 (only uses extensions/compression - no shear)
    for ( j=0 ; x<xMax ; j++ ){
        // Set the new cell size
        Lx=Lx_init*(1.0+x);
        Ly=Ly_init*(1.0-x);
        Lz=Lz_init*(1.0+(x*x)/(1.0-x*x));
        cell.len=vect3(Lx,Ly,Lz);
      
        // Calculate new co-ordinates by scaling fractional coordinates.
        for ( i=0 ; i<n ; i++ ){
            atoms_current[i].pos.x() = atoms_frac[i].pos.x()*Lx;
            atoms_current[i].pos.y() = atoms_frac[i].pos.y()*Ly;
            atoms_current[i].pos.z() = atoms_frac[i].pos.z()*Lz;
        }

        // Work out the energies
        H = zeros<mat>(4*n,4*n);
        evects = zeros<mat>(4*n,4*n);
        evals = zeros<vec>(4*n);
        occupation = zeros<vec>(4*n);

        buildHamiltonian(atoms_current, cell, H);                                                                 
        tb_energy = diagonalise(H,evals,evects,occupation);
        repulsive_energy = repulsiveEnergy(atoms_current, cell);
        
        // Write to file
        c11c12 << x << "\t" << (tb_energy+repulsive_energy)/(cell.volume()) << endl;

        // Re-initialise atoms. 
        atoms_current = atoms_init;
        x+=dx;
    }
    c11c12.close();
}

//------------------------------------------------------------------------
// Function elasticConstant-c44
//
// Params:      - atom list, unit cell, filename (optional).
// Operation:   - Works out the energy as a function of a monoclinic strain
//              and writes it to file. This can then be analysed to obtain 
//              c44.
//------------------------------------------------------------------------
void elasticConstants_c44(atom_list& atoms_init, unitcell& cell_init, string filename="outputs/elasticConst"){
    int n=atoms_init.count(), i, j;
    double Lx, Ly, Lz, x, xMax, dx, alpha;
    double tb_energy, repulsive_energy; 
    mat H, evects, strain;
    vec evals, occupation, lx_curr(3), ly_curr(3), lz_curr(3);
    ofstream c44;
    string c4filename = filename + "-c44.txt";
    c44.open(c4filename.c_str());

    // Set up initial cell of n atoms.
    atom_list atoms_current(n);
    atom_list atoms_frac(n);

    Lx=cell_init.len.x();
    Ly=cell_init.len.y();
    Lz=cell_init.len.z();
 
    unitcell cell=cell_init;
    
    x=-0.05;
    xMax=-1.0*x;
    dx=(xMax-x)/25.0;

    strain = zeros<mat>(3,3);

    writePositions(atoms_current, 0); 

    // Work out c44
    for ( j=0 ; x<xMax ; j++ ){
        // Define strain matrices
        strain(0,0)=1.0;
        strain(0,1)=0.5*x;

        strain(1,0)=0.5*x;
        strain(1,1)=1.0;

        strain(2,2)=1.0+x*x/(4.0-x*x);

        // Apply strain to length vectors.
        lx_curr=strain*vect3(Lx,0.0,0.0);
        ly_curr=strain*vect3(0.0,Ly,0.0);
        lz_curr=strain*vect3(0.0,0.0,Lz);

        alpha = acos(norm_dot(lx_curr,ly_curr));
        cell.len=vect3(norm(lx_curr,2),norm(ly_curr,2),norm(lz_curr,2));
        cell.ang = vect3(alpha, M_PI/2.0, M_PI/2.0);

        cell.calcS();

        // Calculate new co-ordinates by scaling fractional coordinates.
        for ( i=0 ; i<n; i++ ){
              atoms_current[i].pos = vect3(strain*atoms_init[i].pos); 
              atoms_current[i].pos = vect3(cell.Sinv*atoms_init[i].pos); 
        }


        // Supersize the cell 
        writePositions(atoms_current, 1); 

        // Work out the energies
        H = zeros<mat>(4*n,4*n);
        evects = zeros<mat>(4*n,4*n);
        evals = zeros<vec>(4*n);
        occupation = zeros<vec>(4*n);

        buildHamiltonian(atoms_current, cell, H);                                                                 
        tb_energy = diagonalise(H,evals,evects,occupation);
        repulsive_energy = repulsiveEnergy(atoms_current, cell);
        
        // Write to file
        c44 << x << "\t" << (tb_energy+repulsive_energy)/(cell.volume()) << endl;

        // Re-initialise atoms. 
        atoms_current = atoms_init;
        x+=dx;
    }
    c44.close();
}


//------------------------------------------------------------------------
// Function bulkModulus
//
// Params: atomlist, unit cell, filename (optional)
// Operation: Works out the bulk modulus by supersizing the unit
// cell (by an amount specified in scaling factor) and writes the energy 
// as a function of volume  to file so that the bulk modulus can be fitted.
//------------------------------------------------------------------------
void bulkModulus(atom_list& atoms_init, unitcell& cell_init, string filename="outputs/atoms"){
    int n=atoms_init.count(), i, j;
    double Lx, Ly, Lz, x, xMax, dx;
    double Lx_init, Ly_init, Lz_init;
    double tb_energy, repulsive_energy;
    mat H, evects;
    vec evals, occupation;
    ofstream bulk;
    string bfilename = filename + "-bulk.txt";
    bulk.open(bfilename.c_str());

    // Set up initial cell of n atoms.
    atom_list atoms_frac(n);
    atom_list atoms_current(n);
    Lx_init=cell_init.len.x();
    Ly_init=cell_init.len.y();
    Lz_init=cell_init.len.z();
    unitcell cell=cell_init;

    // Obtain fractional co-ordinates.
    for ( i=0 ; i<n ; i++ ){
        atoms_frac[i].pos.x() = atoms_init[i].pos.x()/Lx_init;
        atoms_frac[i].pos.y() = atoms_init[i].pos.y()/Ly_init;
        atoms_frac[i].pos.z() = atoms_init[i].pos.z()/Lz_init;
    }

    x=-0.05;
    xMax=0.05;
    dx=(xMax-x)/15.0;
        //writePositions(atoms_current, 0); 
    // Work out C11-C12 (only uses extensions/compression - no shear)
    for ( j=0 ; x<xMax ; j++ ){
        // Set the new cell size
        Lx=Lx_init*(1.0+x);
        Ly=Ly_init*(1.0+x);
        Lz=Lz_init*(1.0+x);
        cell.len=vect3(Lx,Ly,Lz);
        
        // Calculate new co-ordinates by scaling fractional coordinates.
        for ( i=0 ; i<n ; i++ ){
            atoms_current[i].pos.x() = atoms_frac[i].pos.x()*Lx;
            atoms_current[i].pos.y() = atoms_frac[i].pos.y()*Ly;
            atoms_current[i].pos.z() = atoms_frac[i].pos.z()*Lz;
        }

        // Work out the energies
        H = zeros<mat>(4*n,4*n);
        evects = zeros<mat>(4*n,4*n);
        evals = zeros<vec>(4*n);
        occupation = zeros<vec>(4*n);

        buildHamiltonian(atoms_current, cell, H);                                                                 
        tb_energy = diagonalise(H,evals,evects,occupation);
        repulsive_energy = repulsiveEnergy(atoms_current, cell);
        
        // Write to file
        bulk << cell.volume() << "\t" << (tb_energy+repulsive_energy) << endl;

        // Re-initialise atoms. 
        atoms_current = atoms_init;
        x+=dx;
    }
    bulk.close();
}
#endif
