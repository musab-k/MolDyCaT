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
              evects(zeros<mat>(4*atoms.count(),4*atoms.count())),
              evals(zeros<vec>(4 *atoms.count())),
              occupation(zeros<vec>(4*atoms.count())),
              H(zeros<mat>(atoms.count()*4, atoms.count()*4)),
              dH(zeros<cube>(atoms.count()*4, atoms.count()*4,3)),
              DM(zeros<mat>(atoms.count()*4, atoms.count()*4)),
              forces(atoms.count())
               {}


vec_list force(const atom_list atoms,const unitcell& cell)
{
    double eq = 1.0;
    vect3 d  = atoms[1].pos.diff_bc(atoms[0].pos,cell);
    double r  = d.norm();
    d.normalise();

    vec_list f(2);

    f[0] = static_cast<vect3>(1.0*(r - eq) * d);
    f[1] = static_cast<vect3>(1.0*(eq - r) * d);

    return f;
}

double calc_energy(const atom_list& atoms)
{
    return 0.0;
}

double energy(const atom_list atoms,const unitcell& cell)
{
    double eq = 1.0;
    vect3 d  = atoms[1].pos.diff_bc(atoms[0].pos,cell);
    double r  = d.norm();

    return 0.5 * (r-eq)*(r-eq);

}


            const vec_list& calc_forces(const atom_list& atoms)
            {
               
                eval_count++;                 
                tb_energy = energy(atoms,sim.cell);
                 forces =  force(atoms,sim.cell);
                 return forces;

            }
            const vec_list& calc_phonon_forces(const atom_list& atoms)
            {
               
                eval_count++;                 
                tb_energy = energy(atoms,sim.cell);
                 forces =  force(atoms,sim.cell);
                 return forces;

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


            mat H;
            cube dH;
            mat evects;
            vec evals;
            vec occupation;
            mat DM;

            vec_list forces;

            double repulsive_energy;
            double tb_energy;
            int eval_count;


            
    };

}
#endif //SIM_FORCEEVAL

