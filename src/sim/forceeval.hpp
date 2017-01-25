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
             sim(simul), 
              H(zeros<mat>(atoms.count()*4, atoms.count()*4)),
              dH(zeros<cube>(atoms.count()*4, atoms.count()*4,3)),
              evects(zeros<mat>(4*atoms.count(),4*atoms.count())),
              evals(zeros<vec>(4 *atoms.count())),
              occupation(zeros<vec>(4*atoms.count())),
              DM(zeros<mat>(atoms.count()*4, atoms.count()*4)),
              forces(atoms.count()),
              repulsive_energy(0.0),
              tb_energy(0.0), eval_count(0),max_occupied_index(0)
                             {}

            // Calculates the tb and repulsive ENERGIES only, and returns the total.
            // computational expense of force calculation is saved 
            double calc_energy(const atom_list& atoms){
                //Reset matrices to zero as not sure if your codes assumes these to be initialised to zero or not
                 H.fill(0.0);
                 unitcell cell = sim.cell;
                 buildHamiltonian(atoms, cell, H);
                 tb_energy = diagonalise(H,evals,evects,occupation,max_occupied_index);
                 repulsive_energy = repulsiveEnergy(atoms, cell);
                 return tb_energy+repulsive_energy;
            }
           
            
            vec_list& calc_phonon_forces(const atom_list& atoms)
            {                //Reset matrices to zero as not sure if your codes assumes these to be initialised to zero or not
                 H.fill(0.0);
                 dH.fill(0.0);
                 DM.fill(0.0);
                 forces.fill(0.0);
                 unitcell cell = sim.cell;
                 buildHamiltonian(atoms, cell, H);
 
                 buildHamiltonianDeriv(atoms, cell,dH);
//                 cout << "class" << endl;
                 tb_energy = diagonalise(H,evals,evects,occupation,max_occupied_index);
                 buildDensityMatrixSlowly(atoms,cell,evals,evects,DM,occupation,max_occupied_index);

                 //set tb_energy accordingly

                 repulsive_energy = repulsiveEnergy(atoms, cell);

              //   cout <<get_repulsive_energy()  << endl << get_tb_energy() << endl;



                buildRepulsiveForces(atoms,cell,forces);//cout<<forces[5]<<endl;
               buildTBforces(atoms,cell,dH,DM,forces);//cout<<forces[5]<<endl;
                
                eval_count++;     
                
                  

                 return forces;

            }

            vec_list& calc_forces(const atom_list& atoms)
            {
                //Reset matrices to zero as not sure if your codes assumes these to be initialised to zero or not
                 H.fill(0.0);
                 dH.fill(0.0);
                 DM.fill(0.0);
                 forces.fill(0.0);
                 unitcell cell = sim.cell;
                 buildHamiltonian(atoms, cell, H);
 
                 buildHamiltonianDeriv(atoms, cell,dH);
                 tb_energy = diagonalise(H,evals,evects,occupation,max_occupied_index);
                 buildDensityMatrix(atoms,cell,evals,evects,DM,occupation,max_occupied_index);

                 repulsive_energy = repulsiveEnergy(atoms, cell);
                buildRepulsiveForces(atoms,cell,forces);//cout<<forces[5]<<endl;
               buildTBforces(atoms,cell,dH,DM,forces);//cout<<forces[5]<<endl;
                
                eval_count++;     
                
                  

                 return forces;

            }

            double get_repulsive_energy() const
            {
                 // This actually only returns the total repulsive energy - might want to change.
                return repulsive_energy;
            }

            double get_tb_energy() const
            {
                 // This actually only returns the total tb energy - might want to change.
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
            int max_occupied_index;

            
    };

}
#endif //SIM_FORCEEVAL

