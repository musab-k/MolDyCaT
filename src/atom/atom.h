#ifndef MOLDYCAT_ATOM_ATOM_H
#define MOLDYCAT_ATOM_ATOM_H

#include <vector>
#include <cmath>
#include "armadillo"
#include "vect3.h"
#include "../sim/simulation.hpp"
#include <iostream>

using namespace std;
using namespace arma;

namespace moldycat
{
//-----------------------------------------------------------------------------
// Atoms object
//
// Basically just two vect3 objects for position (pos) and velocity (vel) respectively
//-----------------------------------------------------------------------------
    class atom
    {
    public:
    
        atom();
		atom(vect3 position, vect3 velocity);
		atom(double x, double y, double z, double vx, double vy, double vz);

        vect3 pos;
        vect3 vel;

        friend ostream& operator<< (ostream &out, atom& a);

    };

//-----------------------------------------------------------------------------
// atom list -  just STL vector of atoms 
//
// you can use all the functions usual STL vectors support (google: c++ vector stl)
// or something like it for reference
//
// have used wrapper as we may want to add 
//-----------------------------------------------------------------------------

    class atom_list : public vector<atom>
    {
        public:
			atom_list( int count);

            const int count() const;
            bool check_dist(const unitcell& cell) const;
            void scale_supercell(int sx, int sy, int sz, unitcell& cell);
            double nearestNeighbour(const unitcell& cell) const;

    };


}
#endif //MOLDYCAT_ATOM_ATOM_H

