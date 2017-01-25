#ifndef SIM_FORCEEVAL
#define SIM_FORCEEVAL

#include <vector>
#include <cmath>
#include <vector>
#include "armadillo"
#include <iostream>
#include "../atom/vect3.h"
#include "../atom/atom.h"
#include "unitcell.h"
#include "simulation.hpp"


#include "../modules/moduleB.hpp"
#include "../modules/moduleC.hpp"
#include "../modules/moduleDnew.hpp"
#include "../modules/moduleE.hpp"


using namespace std;
using namespace arma;

namespace moldycat
{
//-----------------------------------------------------------------------------
// force eval object
//
// allows all the technicalities of calculating forces and energies to be encapsulated
// into a single type. Also allows resuse of Hamiltonians and variables to store
// derivatives and other goodies
//
// I think the calc_forces method is the next thing to get working assuming
// modules B, C, D and E are working
// Just change to what you will need to return a vector of forces (tight binding + repulsive)
// you have access to the simulation parameters using sim.[whatever you will need - cell presumably in the short term]
//-----------------------------------------------------------------------------
    class forceeval
    {
        public:
         //Allocate objects once and recycle in successive force evaluations. 
            forceeval(const atom_list& atoms,const simulation& simul) :
             // atoms(init_atoms), 
             sim(simul), 
              tb_energy(0.0),repulsive_energy(0.0), eval_count(0),
              forces(atoms.count())
               {}


vec_list force(const atom_list atoms,const unitcell& cell)
{
    double eq = 1.0;
    vect3 d  = atoms[1].pos.diff_bc(atoms[0].pos,cell);
    double r  = d.norm();
    d.normalise();

    vec_list f(2);

    f[0] = static_cast<vect3>(10000*(r - eq) * d);
    f[1] = static_cast<vect3>(10000*(eq - r) * d);

    return f;
}



const vec_list& calc_forces(const atom_list& atoms)
{
    forces.fill(0.0);
                   
    for(int i = 0; i < atoms.count(); i++)
    {
        for(int j = 0; j < atoms.count(); j++)
        {
            if(i==j) continue;

            vect3 r_v = atoms[i].pos.diff_bc(atoms[j].pos,sim.cell);
            double r2 = r_v.normsq();


            double ri6 = r2*sqrt(r2);
            ri6 *= ri6;
            ri6 = 1.0/ri6;
             forces[i] += static_cast<vect3>((48.0/r2) * (ri6 * ri6 - 0.5*ri6) * r_v);
        }

    }
    calc_energy(atoms);
    return forces;
  }
const vec_list& calc_phonon_forces(const atom_list& atoms)
{
    return calc_forces(atoms);
}
  

  double calc_energy(const atom_list& atoms)
  {
    double energy = 0.0;

    for(int i = 0; i < atoms.count(); i++)
    {
        for(int j = 0; j < i; j++)
        {
            double r = atoms[i].pos.dist_bc(atoms[j].pos,sim.cell);
            double ri6 = r*r*r;
            ri6 *= ri6;
            ri6 = 1.0/ri6;

            energy += (ri6 * ri6 - ri6);
        }
    }

    tb_energy = 4.0*energy;
    return tb_energy;
  }

            double get_repulsive_energy() const
            {
                return repulsive_energy;
            }

            double get_tb_energy() const
            {
                return tb_energy;
            }

            double get_total_energy() const
            {
                return repulsive_energy + tb_energy;
            }

            double get_eval_count() const
            {
                return  eval_count;
            }


            simulation sim;
        private:
            //atom_list atoms;

            vec_list forces;

            double repulsive_energy;
            double tb_energy;
            int eval_count;


            
    };

}
#endif //SIM_FORCEEVAL

