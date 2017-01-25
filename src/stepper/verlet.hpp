#ifndef MOLDYCAT_STEPPER_VERLET_CPP
#define MOLDYCAT_STEPPER_VERLET_CPP 


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
// Verlet class
//
// Basically just module F copied and pasted and shoehorned into the stepperbase schema
//
// NOT TESTED - DO NOT TRUST - NEEDS ATTENTION OF MODULE F AUTHOR
//
//-----------------------------------------------------------------------------
       class verlet :public stepperbase
    {
    public:
    
        verlet(forceeval forcer, const atom_list& atoms,const simulation& sim) 
        :   stepperbase(forcer,sim), 
            old_force(forcer.calc_forces(atoms))
            {
                }

        void step(atom_list& atoms, forceeval& forcer, const simulation& sim)
        {
            step_count++;
            double mass=1244.82175;
            double dt = sim.dt;
            int n = atoms.count();
            
            unitcell cell = sim.cell;

            for (int i=0; i<n; i++)
            {
                   vect3 newpos(atoms[i].vel * dt +(0.5*dt*dt/mass)*old_force[i]);
                    atoms[i].pos.add_bc(newpos,cell);

            } 
           vec_list new_force = forcer.calc_forces(atoms);
           for (int i=0; i<n; i++){
               atoms[i].vel+=0.5*(old_force[i]+new_force[i])*(dt/mass);
            }
            
            old_force = new_force;
        }
        
    private:
        vec_list old_force;
    
        
    };


}

#endif //MOLDYCAT_STEPPER_VERLET_HPP
