#ifndef MOLDYCAT_STEPPER_SPRK_CPP
#define MOLDYCAT_STEPPER_SPRK_CPP 


#include <vector>
#include <cmath>
#include "armadillo"
#include "../atom/vect3.h"
#include "../atom/atom.h"
#include "stepperbase.hpp"
#include "../sim/forceeval.hpp"

using namespace std;


namespace moldycat
{
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
    template<int S>
       class sprk :public stepperbase
    {
    public:
    
        sprk(forceeval forcer, const atom_list& atoms,const simulation& sim, double* A, double* B) 
        :   stepperbase(forcer,sim), force(atoms.count())
            {
                for(int i=0; i < S; i++)
                {
                    a[i] = A[i];
                    b[i] = B[i];
                }
            }

        void step(atom_list& atoms, forceeval& forcer, const simulation& sim)
        {
            step_count++;
            double dt = sim.dt;
            int n = atoms.count();
            for( int stage = 0; stage < S; stage++)
            {
                // Don't recalculate forces if forces aren't required for this stage
                if(b[stage] != 0)
                {
                    if(stage > 0)
                     if( a[stage-1] != 0.0)
                        force = forcer.calc_forces(atoms);

                } // this needs some carefuly attention!

                for(int i = 0; i < n; i++)
                {
                    atoms[i].vel += (b[stage]*dt/mass) * force[i];
                    atoms[i].pos.add_bc(static_cast<vect3>((a[stage]*dt) * atoms[i].vel),sim.cell);
                }
            }
        }
        
    private:
        double h[S]; 
        double a[S];
        double b[S];
        vec_list force;
};
}

#endif //MOLDYCAT_STEPPER_SPRK_HPP
