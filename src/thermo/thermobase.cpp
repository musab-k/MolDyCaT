#ifndef MOLDYCAT_THERMO_THERMOBASE_HPP
#define MOLDYCAT_THERMO_THERMOBASE_HPP

#include <vector>
#include <cmath>
#include "armadillo"
#include "../atom/atom.h"
#include "assert.h"

using namespace std;


namespace moldycat
{
//-----------------------------------------------------------------------------
// Base object for all thermostats - uses same convention as stepper objects
// Will explain carefully on wednesday
//
// Provided functions which may be useful for derived thermostats
// e.g calculating KE and instantaneous T as well as kB constant
//------------------------------------------------------------------------------
    class thermobase
    {
    public:
  
        thermobase(double temp0, double temptol) : T0(temp0), Ttol(temptol),
        kB(8.6173E-5), mass(1244.82175)
        {}

        virtual void apply( atom_list& in) = 0;


        double T0;
        double Ttol;

        double calcKE(const atom_list& atoms) const
        {
            double ke = 0.0;
            for(int i = 0 ; i <atoms.count(); i++)
            {
                ke += atoms[i].vel.normsq();
            }
            return 0.5 * ke * mass;
        }

        const double kB; // Need to figure out kB in our dubious units!!!!!!!!
        const double mass;


         double T(const atom_list& atoms) const
        {
            int dof = atoms.count() * 3 - 3;
            /// set number of degrees of freedom to 3n-3

            return calcKE(atoms) * 2 / (dof * kB);        
        }

    };
}
#endif MOLDYCAT_THERMO_THERMOBASE_HPP

