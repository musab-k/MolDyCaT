#ifndef MOLDYCAT_THERMO_NOTHERMO_CPP
#define MOLDYCAT_THERMO_NOTHERMO_CPP


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
// Do nothing thermostat
//
//-----------------------------------------------------------------------------
    class nothermo :public thermobase 
    {
    public:
    
        nothermo() : thermobase(0.0,0.0)        
        {
        }

        void apply(atom_list& atoms)
        {
        }
    };


}

#endif // MOLDYCAT_THERMO_NOTHERMO_CPP

