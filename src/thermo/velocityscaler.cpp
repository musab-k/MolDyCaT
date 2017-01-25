#ifndef MOLDYCAT_THERMO_VELOCITYSCALER_CPP
#define MOLDYCAT_THERMO_VELOCITYSCALER_CPP


#include <cmath>
#include "../atom/vect3.h"
#include "../atom/atom.h"
#include "thermobase.cpp"



using namespace std;


namespace moldycat
{
//-----------------------------------------------------------------------------
// 
// Velocity scaler class
// 
// Example of object orientated base-derived pattern to use as model
// when coding further thermostats. NOTE: Code is untested and is intended
// more as an example for object orientated programming. Will test once i
// get a working system to play with
//-----------------------------------------------------------------------------
    class velocityscaler :public thermobase
    {
    public:
    
        velocityscaler(double T, double Ttol) : thermobase(T,Ttol)
        {
        }

      // www.mpip-mainz.mpg.de/~andrienk/journal_club/thermostats.pdf
        void apply(atom_list& atoms)
        {
            // current T function defined in thermobase
            double impliedT = T(atoms);

            if(abs(impliedT - T0) > Ttol*T0)
            {
                double scalefac = sqrt(T0/impliedT);
                for(int i = 0; i < atoms.count(); i++)
                    atoms[i].vel *= scalefac;
            }
        }
    };


}

#endif // MOLDYCAT_THERMO_VELOCITYSCALER_CPP

