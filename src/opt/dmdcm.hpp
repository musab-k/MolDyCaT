#ifndef MOLDYCAT_OPT_DMD_CPP
#define MOLDYCAT_OPT_DMD_CPP 


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
// DMD (coupled modes) algorithm implemented as desribed in http://www.sciencedirect.com/science/article/pii/S0021999103003085
//
// On initialisation and after gamma_recal_steps have relapsed, the damping constants, gamma are recalculated
//-----------------------------------------------------------------------------
       class dmdcm :public optbase
    {
    public:
    
        dmdcm(forceeval forcer, const atom_list& atoms,int max_steps, double force_tol, int gamma_recal_steps = 100, double time_step_scaler = 7) :   optbase(forcer,max_steps, force_tol),gamma(0.0), dt(0.1), N(gamma_recal_steps),
        old_force(forcer.calc_forces(atoms)), vel_half(atoms.count()),gamma_steps(gamma_recal_steps),calc_gamma_flag(true), new_force(atoms.count())
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
            
            if(step_count % gamma_steps ==0)
                calc_gamma_flag = true;


            // Need to compute gamma and dt
            if(calc_gamma_flag)
            {      
                //dt=0.1; 
                for (int i=0; i<n; i++)
                {
                    vect3 newpos((0.5*dt*dt/mass)*old_force[i]);
                   atoms[i].pos.add_bc(newpos,cell);
                } 

                int stage = step_count % gamma_steps;
                old_force = forcer.calc_forces(atoms);

                max_force = old_force.max_norm();
                if(max_force < force_tol)                
                {
                    converged = true;
                    return false;
                }                  

                if(old_force.max_norm() < force_tol)
                {
                    converged = true;
                    return false;
                }

                E[stage] =  forcer.get_total_energy();                
                     
                if(stage == 2)
                {
                    gamma = sqrt(0.5*log(abs((E[0]-E[1])/(E[1]-E[2])))) / dt;
                    
                    if(abs((E[0]-E[1])/(E[1]-E[2])) < 1.0)
                    {
                        cerr << "Unable to calculate next damping factor. Consider increasing DMD_RECALC_STEPS above RELAX_STEPS to prevent recalculation of gamma during relaxation. Or start minimisation closer to minimum energy." << endl;
                        return false;
                    }
                    
                //    cout << "E0 - E1 " << (E[0]-E[1]) <<endl;
                //    cout << "E1 - E2 " << (E[1]-E[2]) <<endl;
                //    cout << "gamma " << gamma << endl;
                //    cout << "dt " << dt << endl;  

                   
                  dt = 2*M_PI / (gamma*N);
                    for(int i = 0; i < n; i++)
                    {
                        old_force[i] -= static_cast<vect3>(2*mass*gamma*atoms[i].vel);
                    }
                    calc_gamma_flag = false;
                }
            }
            else
            {

                for(int i = 0; i < n; i++)
                {
                     vect3 newpos(atoms[i].vel * dt +(0.5*dt*dt/mass)*old_force[i]);
                    atoms[i].pos.add_bc(newpos,cell);
                    vel_half[i] = static_cast<vect3>(atoms[i].vel + old_force[i]*(0.5*dt/mass));
                }
                new_force = forcer.calc_forces(atoms);




      //          cout << "Before damping : " << new_force[0] << endl;
                for(int i = 0; i < n; i++)
                {
                    new_force[i] -= static_cast<vect3>((2*gamma*mass/dt)*vel_half[i]);
                    new_force[i] *= (1.0/(1.0+gamma));
                    atoms[i].vel = static_cast<vect3>(vel_half[i]+new_force[i]*(0.5*dt/mass));
                }
      //         cout << "After damping : " << new_force[0] << endl; 
                old_force = new_force;
                
                max_force = new_force.max_norm();
                if(max_force < force_tol)                
                {
                    converged = true;
                    return false;
                }    
            }


            step_count++;
                        
           return true;

        }
        
    private:
        double gamma;
        double dt;
        const double N;    
        double E[3];
        vec_list old_force; 
        
        vec_list vel_half;
        int gamma_steps; 
        bool calc_gamma_flag;
        vec_list new_force;
    };


}

#endif //OPT_DMD_HPP
