#ifndef MOLDYCAT_THERMO_ANDERSEN_CPP
#define MOLDYCAT_THERMO_ANDERSEN_CPP


#include <cmath>
#include "../atom/vect3.h"
#include "../atom/atom.h"
#include "thermobase.cpp"
#include "rand.cpp"


using namespace std;


namespace moldycat
{
//-----------------------------------------------------------------------------
//
//  Andersen thermostat - will destroy your dynamics but at least it
//  is consistent with the canonical ensemble 
//
//-----------------------------------------------------------------------------
    class andersen :public thermobase 
    {
    public:
    
        andersen(double T, double Ttol,double timestep) : thermobase(T,Ttol), nrng(0,sqrt(kB*T/mass)), dt(timestep), nu(0.2)
        {
        }

        //Rplace current velocities with velocities
        //taken from Maxwell-Boltzmann distribution at given temperature
        void apply(atom_list& atoms)
        {
            // current T function defined in thermobase
            double impliedT = T(atoms);

            if(abs(impliedT - T0) > Ttol*T0)
            {
                for(int i = 0 ; i < atoms.count(); i++)
                {
                    if(nrng.next() < nu * dt)
                    for(unsigned int j = 0 ; j < 3; j++)
                    {
                        atoms[i].vel[j] = nrng.norm();
                    }
                }
            }
        }

    private:

        rangauss nrng;
        double dt;
        double nu;
    };


}

#endif // MOLDYCAT_THERMO_ANDERSEN_CPP

