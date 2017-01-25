#ifndef MOLDYCAT_OPT_SD_CPP
#define MOLDYCAT_OPT_SD_CPP 


#include <vector>
#include <cmath>
#include "armadillo"
#include "../atom/vect3.h"
#include "../atom/atom.h"
#include "optbase.hpp"
#include "../sim/forceeval.hpp"
#include <iostream>

using namespace std;

using namespace arma;
namespace moldycat
{
//-----------------------------------------------------------------------------
// Steepest descents optimiser
// Very basic SD algorithm without adaptive step size control
//-----------------------------------------------------------------------------
       class sd :public optbase
    {
    public:
    
        sd(forceeval forcer, const atom_list& atoms,int max_steps, double force_tol,double alpha_step = 1.0e-5,bool auto_rescale_step = false) :   optbase(forcer,max_steps, force_tol), 
        force(atoms.count()),
        alpha(alpha_step),
        auto_alpha(auto_rescale_step),
        prev_energy(0.0)
             {
             }

        // Calculate gradients and move atoms by a constant (alpha) in the opposite direction to the gradient vector.
        // This method is likely to get stuck around the minimum if the step size selected is too large.
        bool step(atom_list& atoms, forceeval& forcer)
        {
             unitcell cell(forcer.sim.cell);
             int n = atoms.count();


            if(step_count >= maxsteps)
            {
                cerr << "Relaxation steps exceeded MAX_RELAXATION_STEPS = " << maxsteps <<  " before reaching specified force tolerance" << endl;
                return false;
            }
            
            force = forcer.calc_forces(atoms);
            double energy = forcer.get_total_energy();
            max_force = force.max_norm();
            if(max_force < force_tol)
            {
                converged = true;
                return false;
            }

            for(int i = 0; i < n; i++)
            {
                atoms[i].pos += static_cast<vect3>(alpha*force[i]);
            }
                      
                      
                      // optionally, if enabled
                      // adjust  step size automatically
                    //  cout << "alpha " << alpha << endl;
            if(auto_alpha && step_count > 0)
            {
                if(energy > prev_energy)
                    alpha *= 0.6;
                else
                    alpha *= 1.1;
            }

            prev_energy = energy;
                        
            step_count++;          
           return true;
        }
        
    private:
        vec_list force;
        double alpha;
        bool auto_alpha;
        double prev_energy;
    };


}

#endif //OPT_SD_HPP
