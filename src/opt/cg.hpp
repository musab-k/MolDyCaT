#ifndef MOLDYCAT_OPT_CG_CPP
#define MOLDYCAT_OPT_CG_CPP 


#include <vector>
#include <cmath>
#include "armadillo"
#include "../atom/vect3.h"
#include "../atom/atom.h"
#include "optbase.hpp"
#include "../sim/forceeval.hpp"
#include <iostream>

using namespace std;


namespace moldycat
{
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
    class cg :public optbase
    {
    public:
    
        cg(forceeval forcer, const atom_list& atoms,int max_steps, double force_tol,double alpha = 0.001) :   optbase(forcer,max_steps, force_tol), 
        force(atoms.count()),
        d(atoms.count()),
        d_old(atoms.count()),
        k(alpha),
        force_old(atoms.count())
        {
        }


        bool step(atom_list& atoms, forceeval& forcer)
        {
             unitcell cell(forcer.sim.cell);
             int n = atoms.count();


            if(step_count >= maxsteps)
            {
                cerr << "Relaxation steps exceeded MAX_RELAXATION_STEPS = " << maxsteps <<  " before reaching specified force tolerance" << endl;
                return false;
            }
            
            if(step_count ==0) 
            {
                d_old.fill(0.0);
            }
            force = forcer.calc_forces(atoms);

            max_force = force.max_norm();
            cout << atoms.count() << endl;
            if(max_force < force_tol)
            {
                converged = true;
                return false;
            }

           
            double gam = 0.0;
            double f_f = 0.0;
            double fold_fold = 0.0;
            

            // Calculate gamma
            for(int i = 0; i < n; i++)
            {
                for(unsigned int j = 0; j < 3; j++)
                {                      
                    f_f += force[i][j]*force[i][j];
                    fold_fold += force_old[i][j]*force_old[i][j];              
                   // h[i][j] = (force[i][j] + gamma[i][j] * h_old[i][j]);
                }
       //         cout << atoms_old[i].pos << endl;
         //       cout << k*h[i] << endl << endl;
            }
           
          
            gam = (step_count == 0 )  ? 0.0 : f_f / fold_fold;
          //  cout << "gamma = " << gam << endl;

            //calculate new d
           for(int i = 0; i < n; i++)
           {
                d[i] = static_cast<vect3>(force[i] + gam * d_old[i]);
            //    cout << "d " << i << " " << d[i] <<endl;
           }

            //update positions
           for(int i = 0; i < n; i++)
            {
                atoms[i].pos.add_bc(static_cast<vect3>(k*d[i]),cell);
            }           
            d_old = d;
            force_old = force;
            
                 
            step_count++;          
           return true;
            
        }
        
    private:
        vec_list force;
        vec_list d;
        vec_list d_old;
        double k;

        vec_list force_old;
        };


}

#endif //OPT_CG_HPP

