#include <vector>
#include <cmath>
#include "armadillo"
#include "vect3.h"
#include <iostream>
#include "atom.h"
#include "../sim/unitcell.h"

#include "iostream"

using namespace std;
using namespace arma;

namespace moldycat
{
//-----------------------------------------------------------------------------
// Atoms object
//
// Basically just two vect3 objects for position (pos) and velocity (vel) respectively
//-----------------------------------------------------------------------------


        atom::atom() : pos(), vel(){};
		atom::atom(vect3 position, vect3 velocity) : pos(position), vel(velocity){};
		atom::atom(double x, double y, double z, double vx, double vy, double vz) : pos(x,y,z) , vel(vx,vy,vz){};

        ostream& operator<< (ostream &out, atom& a) 
        {
            out << "pos = " << a.pos << "\tvel = " << a.vel;
            return out;
        }

//-----------------------------------------------------------------------------
// atom list -  just STL vector of atoms 
//
// you can use all the functions usual STL vectors support (google: c++ vector stl)
// or something like it for reference
//
// have used wrapper as we may want to add 
//-----------------------------------------------------------------------------

    	atom_list::atom_list( int count) : vector<atom>(count){}



        const int atom_list::count() const
        {
            return vector<atom>::size();
        }

//-----------------------------------------------------------------------------
// Rescale atom class for supercell calculations with si repeating units 
// along direction i. Unitcell lengths are scaled accordingly
// Iterate over all output cells and copy atomic positions from unit cell
// provided into super cell array
//-----------------------------------------------------------------------------

        void atom_list::scale_supercell(int sx, int sy, int sz, unitcell& cell)
        {

            int index = 0;
            int old_count = size();
            (*this).resize(size()*sx*sy*sz);

            for(int nx = 0; nx < sx; nx++)
            {
                for(int ny = 0; ny < sy; ny++)//corrected to ny
                {
                    for(int nz = 0; nz < sz; nz++)//corrected to nz
                    {
                        for(int i = 0; i < old_count; i++)
                        {
                            (*this)[index++].pos = static_cast<vect3>((*this)[i].pos +  vect3(nx*cell.len.x(),ny*cell.len.y(),nz*cell.len.z()));
                            //(*this)[index].pos = vect3((*this)[i].pos);
                        }
                    }
                }
            }

            cell.len.x() *= sx;
            cell.len.y() *= sy;
            cell.len.z() *= sz;

            const double min_cell = 5.2;

          if(cell.len[0] < min_cell || cell.len[1] < min_cell || cell.len[2] < min_cell)
           {
               cerr << "Unit cell lengths must be greater than the potential cutoff radius (5.2A)" << endl;
               exit(1);
           }

        }

//-----------------------------------------------------------------------------
// Calculate nearest neihbour distance between atoms (only useful for periodic
// carbon structures)
//-----------------------------------------------------------------------------
        
        double atom_list::nearestNeighbour(const unitcell& cell) const
        {
            int k, n=(*this).count();
            double dmin=1.0E6, dcurrent;
            for( k=1; k<n; k++ ){
                    dcurrent=(*this)[k].pos.dist_bc((*this)[0].pos, cell);
                    if (dcurrent < dmin) {
                        dmin=dcurrent;
                    }
            }
            if(dmin<1.0E6){ return dmin; }
            else{ 
                cerr << "nearest neighbour distance incorrect!" << endl;
                return 1;
            }
                    
        }
//-----------------------------------------------------------------------------
// Check that all atoms are at least 0.8A apart. If atoms are much closer together
// they will have unphysically high energies (corresponds to 50% bond deformation)
//-----------------------------------------------------------------------------
        bool atom_list::check_dist(const unitcell& cell) const
        {
           const double minsep = 0.8;

            // Double loop over atom pairs to check if atoms are too close
            // (after application of PBCs)
            for(unsigned int i = 0; i < size(); i++)
            {
                for(unsigned int j = i+1; j < size(); j++)
                {

                    double dist = (*this)[i].pos.dist_bc((*this)[j].pos,cell);

                    if(dist < minsep)
                    {
                        cerr << "Atoms " << i+1 << " and " << j+1 << " are unphysically close together (" << dist << " A)" << endl;
                       
                        return false;
                    }
                }
            }
            return true;
        }


}

