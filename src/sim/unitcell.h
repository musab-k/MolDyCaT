#ifndef MOLDYCAT_SIM_UNITCELL_H
#define MOLDYCAT_SIM_UNITCELL_H

#include <vector>
#include <cmath>
#include <vector>
#include "armadillo"
#include <iostream>
#include "../atom/vect3.h"

using namespace std;
using namespace arma;


namespace moldycat
{
    enum bctype{ Open = 0, Periodic = 1};
    
//    class vect3;

 //-----------------------------------------------------------------------------
// Non-orthorhombic unit cell. Constructor specified using 3 lattice vector lengths
// and angles between using standard crystallographic conventions and a list of boundary conditions. 
// These can be specified using bctype enum. e.g 
// bctype bcs[3] = {Periodic, Open, Open} specifies periodic BCs along 1st and 3rd directions and
// open boundaries along 2nd lattice vector direction. Precalculate 
// transformation matrix S  which transforms non-orthorhombic vectors into orthorhombic
// cell with lattice length A,B and C. Also calcualte and store inverse Sinv.
// Apply BCs in transformed, orthorhombic space and then transform back
// Method based on that of http://zeolites.cqe.northwestern.edu/Music/nonorthorhombic.pdf
//-----------------------------------------------------------------------------
	class unitcell
    {
        public:
            unitcell(double A, double B, double C, double alpha, double beta, double gamma,bctype* bcs);
            //Fields to store transformation matrix, inverse, lattice sides, angles between and boundary contid
            mat S;
            mat Sinv;
            vect3 len;
            vect3 ang;
            bctype bc[3];

            void calcS();
        
            double dot(const vect3 v1,const vect3 v2) const;
            double norm(vect3 v) const;  
         
            vect3 diff_bc(const vect3& v1,const vect3& v2) const;
             double dist_bc(const vect3& v1,const vect3& v2) const;

            void add_bc(vect3& v1,const vect3& v2) const;
            double volume() const;      

    };

}
#endif MOLDYCAT_SIM_UNITCELL_H


