#ifndef MOLDYCAT_SIM_SIMULATION_H
#define MOLDYCAT_SIM_SIMULATION_H

#include <vector>
#include <cmath>
#include "armadillo"
#include <iostream>
#include "../atom/vect3.h"
#include "unitcell.h"

using namespace std;
using namespace arma;


namespace moldycat
{
 //-----------------------------------------------------------------------------
// Simulation object contains all the parameter values specified in input file.
// by encapsulating them in a class, we can reduce the number of variable that
// must be passed between functions which can make things unclear.
//
// We also enforce const-correctness by preventing any other part of the code
// from modifying user input values accidentially during the simulation
//
//-----------------------------------------------------------------------------

	class simulation
    {
        public:
        simulation(string SEED, unitcell ucell, int ENERGY_OUT_STEPS, int VELOCITY_OUT_STEPS, int POSITIONS_OUT_STEPS, double TIME_STEP, int MD_STEPS,int RELAX_STEPS, double SKIN_DEPTH, double TEMPERATURE, double TEMP_TOL, string STRUCT_RELAX, string RELAXER, double FORCE_TOL, string THERMOSTAT, string INTEGRATOR, string SUPERCELL, int S_x, int S_y, int S_z, string RESTART,int DMD_GAMMA_RECALC_STEPS,double SD_STEP_SIZE,double DMD_TIMESTEP_SCALE, string TASK,vec_list phonon_q_point_list,double PHONON_DISPLACEMENT,int RDF_BINS,double RDF_CUTOFF,int RDF_STEPS) : 
        seed(SEED),
        cell(ucell), 
        temp(TEMPERATURE), 
        temptol(TEMP_TOL),
        dt(TIME_STEP),
        skin_depth(SKIN_DEPTH),  
        velocity_out_steps(VELOCITY_OUT_STEPS),
        energy_out_steps(ENERGY_OUT_STEPS), 
        positions_out_steps(POSITIONS_OUT_STEPS), 
        sx(S_x), sy(S_y), sz(S_z),
        md_steps(MD_STEPS),
        relax_steps(RELAX_STEPS),  
        struct_relax(STRUCT_RELAX), 
        relaxer(RELAXER),
        force_tol(FORCE_TOL),
        thermostat(THERMOSTAT), 
        integrator(INTEGRATOR), 
        supercell(SUPERCELL),
        restart(RESTART),
        dmd_gamma_recalc_steps(DMD_GAMMA_RECALC_STEPS),
        sd_step_size(SD_STEP_SIZE),
        dmd_timestep_scale(DMD_TIMESTEP_SCALE),
        task(TASK),
        phonon_displacement(PHONON_DISPLACEMENT),
        phonon_q_points(phonon_q_point_list),
        rdf_bins(RDF_BINS),
        rdf_cutoff(RDF_CUTOFF),
        rdf_steps(RDF_STEPS)
               { }


        simulation(const simulation& sim) : 
        seed(sim.seed),
        cell(sim.cell), 
        temp(sim.temp), 
        temptol(sim.temptol),
        dt(sim.dt),
        skin_depth(sim.skin_depth),  
        velocity_out_steps(sim.velocity_out_steps),
        energy_out_steps(sim.energy_out_steps), 
        positions_out_steps(sim.positions_out_steps), 
        sx(sim.sx), sy(sim.sy), sz(sim.sz),
        md_steps(sim.md_steps),
        relax_steps(sim.relax_steps),  
        struct_relax(sim.struct_relax), 
        relaxer(sim.relaxer),
        force_tol(sim.force_tol),
        thermostat(sim.thermostat), 
        integrator(sim.integrator), 
        supercell(sim.supercell),
        restart(sim.restart),
        dmd_gamma_recalc_steps(sim.dmd_gamma_recalc_steps),
        sd_step_size(sim.sd_step_size),
        dmd_timestep_scale(sim.dmd_timestep_scale),
        task(sim.task),
        phonon_displacement(sim.phonon_displacement),
        phonon_q_points(sim.phonon_q_points),
        rdf_bins(sim.rdf_bins),
        rdf_cutoff(sim.rdf_cutoff),
        rdf_steps(sim.rdf_steps)
               { }

            const string seed;
            const unitcell cell;

            const double temp;
            const double temptol;
            const double dt;
            const double skin_depth;
            const int velocity_out_steps;
            const int energy_out_steps;
            const int positions_out_steps;
            const int sx;
            const int sy;
            const int sz;
            const int md_steps;
            const int relax_steps;
            const string struct_relax;
            const string relaxer;
            const double force_tol;
            const string thermostat;
            const string integrator;
            const string supercell;
            const string restart;
            const int dmd_gamma_recalc_steps;
            const double sd_step_size;
            const double dmd_timestep_scale;
            const string task;
            const double phonon_displacement;
            const vec_list phonon_q_points;
            const int rdf_bins;
            const int rdf_cutoff;
            const int rdf_steps;        
    };

}
#endif MOLDYCAT_SIM_SIMULATION_H



