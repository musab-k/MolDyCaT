#ifndef RADIAL_HPP
#define RADIAL_HPP

#include <vector>
#include <cmath>
#include <vector>
#include "armadillo"
#include "string.h"
#include <iostream>
#include "../atom/vect3.h"
#include "../atom/atom.h"
#include "../sim/unitcell.h"
#include "../sim/simulation.hpp"


 void binsetup (int noBins, double cutoff, colvec& r, colvec& bins){
     double width = cutoff/noBins;
     r = zeros<colvec>(noBins);
     bins = zeros<colvec>(noBins);
     for(int i=0;i<noBins;i++){
         r(i)=(width*i)+width;
     }
 }

 void addtobins(atom_list& atoms, const unitcell& cell, colvec& bins, double cutoff, int& counter){
        counter++;
        unitcell newcell = cell;
        atom_list newatoms = atoms;
        newatoms.scale_supercell(3,3,3,newcell);
     for(int i=0;i<newatoms.count();i++){
         for(int j=i+1;j<newatoms.count();j++){
             double dist = newatoms[j].pos.dist_bc(newatoms[i].pos,newcell);
             if(dist<cutoff){
                 int bin = int(ceil((dist/cutoff)*bins.size()))-1;
                 bins(bin)+=2;
                 }
         }
     }
 }

 void printradial(colvec& r, colvec& bins, int& counter, const unitcell& cell, int n, string filename){
    ofstream outfile;
    outfile.open(filename.c_str());
    n*=27;
    double volume=27*cell.volume();
    double rho=n/volume;
    double width = r(1)-r(0);
    double R;

    for(unsigned int i=0;i<bins.size();i++){
         R=r(i)-(width/2);
         bins(i)/=(R*R*width*4*M_PI*n*rho);
         bins(i)/=counter;
         outfile<<R<<"\t"<<bins(i)<<endl;
    }
     outfile.close();
 }

             











#endif
