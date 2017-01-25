// 	          \`*-.                    
// 	           )  _`-.                 
// 	          .  : `. .                
// 	          : _   '  '
// 	          ; *` _.   `*-._          
// 	          `-.-'          `-.       
// 	            ;       `       `.     
// 	            :.       .        '
// 	            . \  .   :   .-'   .   
// 	            '  `+.;  ;  '      :   
// 	            :  '  |    ;       ;-. 
// 	            ; '   : :`-:     _.`* ;
// 	 [MolDy] .*' /  .*' ; .*`- +'  `*' 
// 	         `*-*   `*-*  `*-*'        


#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "atom/atom.h"
#include "atom/vect3.h"
#include "sim/unitcell.h"
#include "sim/simulation.hpp"

#include "sim/forceeval.hpp"
#include <string>
// Uncomment this line to use lennard jones
// rather than TB
//#include "sim/forceeval_lj.hpp"


#include "modules/moduleA.hpp"
#include "modules/moduleG.hpp"

#include "stepper/stepperbase.hpp"
#include "stepper/verlet.hpp"
#include "stepper/sprk.hpp"
#include "stepper/factory.hpp"

#include "thermo/thermofactory.hpp"
#include "thermo/thermobase.cpp"
#include "thermo/velocityscaler.cpp"
#include "thermo/andersen.cpp"
#include "thermo/rand.cpp"
#include "anal/radial.hpp"
#include "anal/elasticConstants.hpp"
#include "anal/graphiteGenerator.hpp"
#include "anal/evStretcher.hpp"
#include "opt/optbase.hpp"
#include "opt/factory.hpp"

#include "sim/phonon.hpp"

using namespace moldycat;
using namespace arma;


void calc_phonons( atom_list& atoms, forceeval& forcer, const simulation& sim);


int main(int argc, char* argv[])
{
    cout << setprecision(9);


    // Allocate pointer to hold simulation object which will contain
    // all the parameters specified in the user input files
    simulation* sim_ptr;
    // An atoms list objet is created from the user specified positions
    // file and the pointer set to point at a newly allocated simulation object
    atom_list atoms = read_input_file(argc,argv,sim_ptr);
    // This simulation is dereferenced and a local copy is made to be passed 
    // to other parts of the code
    simulation sim=*sim_ptr; 

    cout << "Simulation initialised with " << atoms.count() << " atoms." <<endl; 
    writePositions(atoms,0,("outputs/" + sim.seed +"_positions.xyz").c_str());

    // Bins are initialised to for sampling the RDF
    int step=0;
    int counter=0;
    colvec r, bins;
    binsetup(sim.rdf_bins,sim.rdf_cutoff,r,bins);
    
    // With the simulation (and therefore unticell) known
    // and the atom positions specified a forcer object
    // may be allocated and a local copy made
    // for 'syntactic sugar'
    forceeval* forcer_ptr = new forceeval(atoms,sim);
    forceeval forcer = *forcer_ptr;

    // If a structural relaxation is requested do this first
    // otherwise skip this section
    if(sim.struct_relax=="ON"){
      // open a filestream to output forces and energies as the relaxation progresses
      // and another to output the updated atomic positions
      ofstream relaxfile;
      relaxfile.open(("outputs/" + sim.seed + "_relax.out").c_str());

      // write initial positions to input file
      writePositions(atoms,0,("outputs/" + sim.seed +"_opt.xyz").c_str());
      // Instantialise a new optbase type to the derived specified in the input file
      // the optFactory switches the sim.relaxer flag to determine if
      // a supported optimisation algorithm is specified
      optbase* relaxer = optFactory(forcer,atoms,sim);
     
        // The relaxation algorithm is stepped until
        // the maximum number of relaxation steps is reached,
        // the maximum force on any atoms is lower that the 
        // specified force tolerance or an error occours in
        // which case an error is output to the user
        while(relaxer->step(atoms,forcer))
        {
            writePositions(atoms,1,("outputs/" + sim.seed +"opt.xyz").c_str());
            cout << "Step\t" << step++
                 << "\nEnergy\t" << forcer.get_total_energy() << " eV"
                << "\nMax. Force\t" << relaxer->get_max_force() << " eV/A" << endl;
           relaxfile << step << '\t' << forcer.get_total_energy() << '\t' << relaxer->get_max_force() << endl;
        }
        relaxfile.close();
    
            // Test if the stop condition is the result of successful convergence
            // If so continue without stopping. Otherwise prompt user to see
            // whether or not to continue with the simulation
            if(relaxer->is_converged()){
                cout << "Geometery relaxation successful. Forces converged within " << sim.force_tol << " eV/A" << endl;
            }
            else{
                cout << "Forces not fully converged. Continue with "<<sim.task<<" simulation? (Y or N)" << endl;
                char s;
                cin >> s;
                if(s != 'Y' && s != 'y'){ cout << "Aborting\n"; return 0;};
            }
            // free memory associated with relaxer
            delete relaxer;
     }



    // See what type of simulation requested (ENERGY, 
    // PHONON or MD) and determine program flow accordingly
    if(sim.task == "PHONON" || sim.task == "PHONONS")
    {
        cerr << "Phonons feature is deprecated." << endl;
        exit(1);
        calc_phonons(atoms,forcer,sim);
        return 0;
    }

    if(sim.task == "ENERGY" || sim.task == "SINGLE_POINT")
    {   
        cout << "System energy " << forcer.calc_energy(atoms) << " eV" << endl;
        return 0;
    }

    // Instantiate integrators and thermostat using factory methods
    // as defined in user input 
    stepperbase* integrator =  integratorFactory(sim,atoms,forcer);
    thermobase* thermostat =   thermoFactory(sim); 


    // Prepare filestreams for simulation output
    ofstream energyfile;
    ofstream tempfile;
    energyfile.open(("outputs/" + sim.seed +"_energy.out").c_str());
    tempfile.open(("outputs/" + sim.seed +"_temp.out").c_str());
    
    writeVelocities(atoms,0,("outputs/" + sim.seed +"_velocities.out").c_str());
    
    for(step = 0; step < integrator->maxsteps; step++)
        {
       // Iterate integrator by one step using the forces returned by
       // forcer at the current atomic positions
            integrator->step(atoms,forcer,sim);
       //Get energies and T from appropriate simulation objects and output to user and write to a file     
            double total_energy =  forcer.get_total_energy() + thermostat->calcKE(atoms);
                      cout << "Step " << step << "\n\tTB energy\t" << forcer.get_tb_energy() << endl
                                    << "\tRep energy\t" << forcer.get_repulsive_energy() << endl
                                    << "\tTot energy\t" << total_energy << endl
                                    << "\tT\t" << thermostat->T(atoms) << endl;

            tempfile << step << " " << thermostat->T(atoms) << endl;

            // Apply thermostat to atomic velocities
            thermostat->apply(atoms);
            
            // Write details of current simulation state to output files if requested
            // (i.e sim.X_steps == 0) every sim.X_steps 
            if(sim.energy_out_steps != 0 && step%sim.energy_out_steps==0){
                  energyfile << step  << " " << forcer.get_total_energy()
                                << " " << total_energy << endl;
                                writePositions(atoms,0,("outputs/" + sim.seed +"_positions.end").c_str());
                                writeVelocities(atoms,0,("outputs/" + sim.seed +"_velocities.end").c_str());
                                writeUnitcell(sim,("outputs/" + sim.seed +"_unitcell.out").c_str());
            }
            if(sim.positions_out_steps != 0 && step%sim.positions_out_steps==0){
                writePositions(atoms, 1,("outputs/" + sim.seed +"_positions.xyz").c_str());
                }

            if(sim.rdf_steps != 0 && step%sim.rdf_steps==0){
                addtobins(atoms,sim.cell,bins,sim.rdf_cutoff,counter);
                }

            if(sim.velocity_out_steps != 0 && step%sim.velocity_out_steps==0){
                writeVelocities(atoms, 1,("outputs/" + sim.seed +"_velocities.out").c_str());
                }
                            
   }
   // Write final positions and velocities to a file for restart

   writePositions(atoms,0,("outputs/" + sim.seed +"_positions.end").c_str());
   writeVelocities(atoms,0,("outputs/" + sim.seed +"_velocities.end").c_str());
   writeUnitcell(sim,("outputs/" + sim.seed +"_unitcell.out").c_str());
   energyfile.close();
   tempfile.close();

   // If requested, print radial distribution function to a file
   if(sim.rdf_steps > 0){ printradial(r,bins,counter,sim.cell,atoms.count(),("outputs/" + sim.seed +"_rdf.dat").c_str());}
    // free dynamically allocated memory

    delete integrator;
    delete thermostat;

    delete forcer_ptr;
    delete sim_ptr;

   return 0;

}
// phonon calculation was sadly deprecated
void calc_phonons( atom_list& atoms, forceeval& forcer, const simulation& sim)
{/*
   ofstream phonfile;
   phonfile.open(("outputs/" + sim.seed +"_phonon.freqs").c_str());

   int n = atoms.count()/(sim.sx*sim.sy*sim.sz);
   cout << n << endl;
   cx_mat DyMat = zeros<cx_mat>(3*n,3*n);
   for(int i = 0; i < sim.phonon_q_points.count(); i++)
   {
        vect3 q_pt = sim.phonon_q_points[i];
        cx_vec evals(3*n);
        cx_mat evects(3*n,3*n);
        buildDynamicalMatrix(atoms,evects,evals, forcer, q_pt,sim); 
        //cout << evals;
        //cout << evects;
        phonfile << "q-point:\t" <<  q_pt.x() << '\t' << q_pt.y() << '\t' << q_pt.z() << endl;
        for(int j = 0; j < 3 * n; j++)
        {
            phonfile << real(evects[j]) << endl;
        }
   }
   phonfile.close();*/
} 
