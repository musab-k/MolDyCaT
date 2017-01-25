//-----------------------------------------------------------------------------
// Module A
//
// Description - Does the parsing of the input file and reads the positions
//               file
//-----------------------------------------------------------------------------
#ifndef MODULEA_HPP
#define MODULEA_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <cmath>

#include "../atom/atom.h"
#include "armadillo"
#include "../atom/vect3.h"
#include "../sim/unitcell.h"
#include "../sim/simulation.hpp"

using namespace moldycat;
using namespace arma;

// Declare the function "to_upper"
string to_upper(string);

// Declare function to read in positions and return an atom_list
atom_list readpositions(string,const unitcell&);

// Declare function to read in velocities and return a vect_list
vec_list readvelocities(string);

// Declare a function to read in the unitcell
unitcell readunitcell(string);

//-----------------------------------------------------------------------------
// Function read_input_file 
//
// Returns: atom_list containing the initial positions of the atoms
//
// Params: argc(int)     - The number of arguments from the command line
//         argv[](char*) - The arguments from the command line
//         sim_ptr(simulation*&) - The pointer to a simulation object
// 
// Operation - atom_list is filled in by calling the function readpositions.
//             The unitcell of the positions is required to perform a check
//             that the positions are not unphysically close together.
//-----------------------------------------------------------------------------
atom_list read_input_file(int argc, char* argv[], simulation*& sim_ptr)
{
  string seed;
  // Read in the seed name.
  if (argc==1){
    cout << "Please enter a seed name. The name of the input file should be {seed}_params.in" << endl;
    cin >> seed;
  }
  else// assume that the second arguement is the string
      // this may change to something like --seed {seed}
  {
    seed = argv[1];

#ifdef DBG
    cout << "argv[1]=" << argv[1] << endl;
#endif

  }
  // Temporarily output the name of the seed

#ifdef DBG
  cout << "Outputting the name of the seed." << endl;
  cout << seed << endl;
#endif

  // -----------------------------------------------------------
  // Start reading in the input file
  // It is taken to called {seed}_params.in
  // -----------------------------------------------------------
  ifstream in;
  ofstream outfile;
  string paramfile = "inputs/" + seed + "_params.in";
  string paramsout = "outputs/" + seed + "_params.out";
  in.open(paramfile.c_str());
  outfile.open(paramsout.c_str());
  if(!in){//if input file not opened, error message
    cerr << "ERROR: The input parameters file could not be opened. \n";
    exit (1);
  }
  else{

#ifdef DBG
  cout << "The parameters file has been opened." << endl;
#endif

  }
  string name, NAME, ATOMS_FILE_ADDRESS,first_letter,letter,bound_x_1,bound_x_2,bound_x_3, STRUCT_RELAX, THERMOSTAT, INTEGRATOR, SUPERCELL,RELAXER,RESTART, TASK;
  int ENERGY_OUT_STEPS, VELOCITY_OUT_STEPS,POSITIONS_OUT_STEPS,MD_STEPS,RELAX_STEPS, S_x, S_y, S_z,DMD_GAMMA_RECALC_STEPS,RDF_BINS, RDF_STEPS;
  int maxsize = 1024; // The maximum expected line length (in chars)
  char skip[maxsize];
  double x_1, x_2, x_3, TIME_STEP, SKIN_DEPTH, TEMPERATURE, TEMP_TOL,alpha,beta,gamma,FORCE_TOL,SD_STEP_SIZE,DMD_TIMESTEP_SCALE,PHONON_DISPLACEMENT,RDF_CUTOFF;
  int i=0,IMAX=150;

  
  bool in_phonon_q_list = false;
  vec_list phonon_q_points(0);

  // Set the default values

  ENERGY_OUT_STEPS=100;
  VELOCITY_OUT_STEPS=100;
  POSITIONS_OUT_STEPS=100;
  TIME_STEP=0.1;
  MD_STEPS=1000;
  RELAX_STEPS=100;
  x_1=6.0;// x_1 length of unit cell
  x_2=6.0;// x_2 length of unit cell
  x_3=6.0;// x_3 length of unit cell
  bound_x_1=to_upper("Open");
  bound_x_2=to_upper("Open");
  bound_x_3=to_upper("Open");
  alpha=M_PI/2;
  beta=M_PI/2;
  gamma=M_PI/2;
  SKIN_DEPTH=1.0;
  TEMPERATURE=350;
  TEMP_TOL=0.01;
  FORCE_TOL=1e-2;
  PHONON_DISPLACEMENT=1e-5;
  ATOMS_FILE_ADDRESS=seed+"_atoms.in";
  STRUCT_RELAX=to_upper("Off");
  RELAXER=to_upper("SD");
  RESTART="OFF";
  TASK="ENERGY";
  THERMOSTAT=to_upper("NOTHERMO");
  INTEGRATOR=to_upper("VERLET");
  SUPERCELL=to_upper("N");
  RDF_BINS=100;
  RDF_CUTOFF=3.0;
  RDF_STEPS=100;
  S_x=1;
  S_y=1;
  S_z=1;

  DMD_GAMMA_RECALC_STEPS=100000;
  SD_STEP_SIZE = 0.001;
  DMD_TIMESTEP_SCALE = 7;

  while(!in.eof())
  {

    i++;
    if (i>IMAX){

#ifdef VERBOSE
      cerr << "The while loop has been broken by IMAX." << endl;
#endif
     
      cerr << "Problem reading the input file. Problem occured directly after reading the word " << name << endl;
      exit(1);
    }

    if(in_phonon_q_list) {
      string token;
      in >> token;
      token = to_upper(token);
      
      if (token == "END"){
        in >> token;
        in_phonon_q_list = false; continue; 
      }

      double valx,valy,valz;
      valx = atof(token.c_str());
      in >> valy;
      in >> valz;
      vect3 q_pt(valx,valy,valz);
      phonon_q_points.push_back(q_pt);
      continue;
    }


    in >> name;
    first_letter = name[0];//first letter of string
 
    NAME=to_upper(name);

    // If there is a comment we skip it
    if (first_letter == "#" || first_letter == "!"){
      in.getline(skip,maxsize);
      continue;
    }

    // Otherwise see if we recognise token and try to parse result
    // If the >> operator fails and we continue to iterate through
    // until i_max lines of code reached at which point we exit the
    // loop and assume there was an issue and complain to user
    else if (NAME =="ENERGY_OUT_STEPS"){
      in >> ENERGY_OUT_STEPS;//set Energy to that value
    }

    else if (NAME =="VELOCITY_OUT_STEPS" ){
      in >> VELOCITY_OUT_STEPS;
    }

    else if (NAME == "POSITIONS_OUT_STEPS"){
      in >> POSITIONS_OUT_STEPS;
    }

    else if (NAME == "TIME_STEP"){
      in >> TIME_STEP;
    } 

    else if (NAME == "RELAX_STEPS"){
      in >> RELAX_STEPS;
    } 

    else if (NAME == "MD_STEPS"){
      in >> MD_STEPS;
    }

    else if (NAME == "CELL_DIMENSIONS"){
      in >> x_1;
      in >> x_2;
      in >> x_3;
    }

    else if (NAME == "BOUNDARY_TYPE"){
      in >> bound_x_1;
      in >> bound_x_2;
      in >> bound_x_3;
      bound_x_1 = to_upper(bound_x_1);
      bound_x_2 = to_upper(bound_x_2);
      bound_x_3 = to_upper(bound_x_3);
    }

    else if (NAME == "CELL_ANGLES"){
      in >> alpha;
      in >> beta;
      in >> gamma;
    }

    else if (NAME == "SKIN_DEPTH" ){
      in >> SKIN_DEPTH;
    }

    else if (NAME == "TEMPERATURE"){
      in >> TEMPERATURE;
    }

    else if (NAME == "TEMP_TOL"){
      in >> TEMP_TOL;
    }
    
    else if (NAME == "FORCE_TOL"){
      in >> FORCE_TOL;
    }

    else if (NAME == "ATOMS_FILE_ADDRESS"){
      in >> ATOMS_FILE_ADDRESS;
    }

    else if (NAME == "STRUCT_RELAX"){
      in >> STRUCT_RELAX;
      STRUCT_RELAX = to_upper(STRUCT_RELAX);
    }

    else if (NAME == "RELAXER"){
      in >> RELAXER;
      RELAXER = to_upper(RELAXER);
    }

    else if (NAME == "THERMOSTAT"){
      in >> THERMOSTAT;
      THERMOSTAT = to_upper(THERMOSTAT);
    }

    else if (NAME == "INTEGRATOR"){
      in >> INTEGRATOR;
      INTEGRATOR = to_upper(INTEGRATOR);
    }

    else if (NAME == "RESTART"){
      in >> RESTART;
      RESTART = to_upper(RESTART);
    }

    else if (NAME == "TASK"){
      in >> TASK;
      TASK = to_upper(TASK);
    }

    else if (NAME == "SUPERCELL"){
      in >> SUPERCELL;
      SUPERCELL=to_upper(SUPERCELL);
      if (SUPERCELL!="Y" && SUPERCELL!="N"){
        SUPERCELL="N";// i.e. if not specified properly turn off.
      }
    }

    else if (NAME == "MULTIPLY_CELL_DIRN"){
      in >> S_x;
      in >> S_y;
      in >> S_z;
    }

    else if (NAME == "DMD_GAMMA_RECALC_STEPS"){
      in >> DMD_GAMMA_RECALC_STEPS;
    } 

    else if (NAME == "SD_STEP_SIZE"){
      in >> SD_STEP_SIZE;
    } 
    
    else if (NAME == "DMD_TIMESTEP_SCALE"){
      in >> DMD_TIMESTEP_SCALE;
    }
    
    else if (NAME == "RDF_BINS"){
      in >> RDF_BINS;
    }
    
    else if (NAME == "RDF_CUTOFF"){
      in >> RDF_CUTOFF;
    }
    
    else if (NAME == "RDF_STEPS"){
      in >> RDF_STEPS;
    }
    
    else if (NAME == "PHONON_DISPLACEMENT"){
      in >> PHONON_DISPLACEMENT;
    }
    
    else if (NAME == "PHONON_Q_POINTS"){
      in_phonon_q_list = true;
      in >> name;
      in >> name;
    }

    else {
      cerr << "Error on line " << i << ". Invalid token '" << name << "'" << endl;
      exit(1);
    }
  }

  if(in_phonon_q_list)
  {
    cerr << "Expected END_PHONON_Q_POINT_LIST token." << endl;
    exit(1); 
  }  

  //Print out to params.out depending on type of simulation specified
  outfile << "TASK" << '\t' << TASK << '\n';
  outfile << "RESTART" << '\t' << RESTART << '\n';  
  outfile << "CELL_DIMENSIONS" << '\t' << x_1 << '\t' << x_2 << '\t' << x_3 << '\n';
  outfile << "BOUNDARY_TYPE" << '\t' << bound_x_1 << '\t' << bound_x_2 << '\t' << bound_x_3 << '\n';
  outfile << "CELL_ANGLES" << '\t' << alpha << '\t' << beta << '\t' << gamma <<'\n';
  outfile << "ATOMS_FILE_ADDRESS" << '\t' << ATOMS_FILE_ADDRESS << '\n';
  outfile << "SUPERCELL" << '\t' << SUPERCELL << '\n';
  outfile << "MULTIPLY_CELL_DIRN" << '\t' << S_x << '\t' << S_y << '\t' << S_z << '\n';

  if(TASK == "MD" || TASK == "MOLECULAR_DYNAMICS")
  {
    outfile << "THERMOSTAT" << '\t' << THERMOSTAT << '\n';
    outfile << "INTEGRATOR" << '\t' << INTEGRATOR << '\n';
    outfile << "ENERGY_OUT_STEPS" << '\t' << ENERGY_OUT_STEPS << '\n';
    outfile << "VELOCITY_OUT_STEPS" << '\t' << VELOCITY_OUT_STEPS << '\n';
    outfile << "POSITIONS_OUT_STEPS" << '\t' << POSITIONS_OUT_STEPS << '\n';
    outfile << "TIME_STEP" << '\t' << TIME_STEP << '\n';
    outfile << "MD_STEPS" << '\t' << MD_STEPS << '\n';
    outfile << "SKIN_DEPTH" << '\t' << SKIN_DEPTH << '\n';
    outfile << "TEMPERATURE" << '\t' << TEMPERATURE << '\n';
    outfile << "TEMP_TOL" << '\t' << TEMP_TOL << '\n';
    outfile << "RDF_BINS" << '\t' << RDF_BINS << '\n';
    outfile << "RDF_CUTOFF" << '\t' << RDF_CUTOFF << '\n';
    outfile << "RDF_STEPS" << '\t' << RDF_STEPS << '\n';
  }
  else if(TASK == "PHONON" || TASK == "PHONONS") {
    outfile << "PHONON_Q_POINTS\nBEGIN LIST" << '\n';
    
    for(int w = 0; w < phonon_q_points.count(); w++)
    {
      outfile << phonon_q_points[w].x() << '\t' << phonon_q_points[w].y() << '\t' << phonon_q_points[w].z() << '\n';
    }
    outfile << "END LIST" << "\n#\n";
    outfile << "PHONON_DISPLACEMENT" << '\t' << PHONON_DISPLACEMENT << '\n'; 
  }

  outfile << "STRUCT_RELAX" << '\t' << STRUCT_RELAX << '\n';
  
  if(STRUCT_RELAX == "YES" || STRUCT_RELAX == "ON")
  {
    outfile << "RELAXER" << '\t' << RELAXER << '\n';
    outfile << "DMD_GAMMA_RECALC_STEPS" << '\t' << DMD_GAMMA_RECALC_STEPS << '\n';
    outfile << "SD_STEP_SIZE" << '\t' << SD_STEP_SIZE << '\n';
    outfile << "DMD_TIMESTEP_SCALE" << '\t' << DMD_TIMESTEP_SCALE << '\n';
    outfile << "RELAX_STEP" << '\t' << TIME_STEP << '\n'; 
    outfile << "FORCE_TOL" << '\t' << FORCE_TOL << '\n';   
  }

  in.close();
  outfile.close();
  bctype boundary[3]={bound_x_1 == "OPEN" ? Open : Periodic,
                      bound_x_2 == "OPEN" ? Open : Periodic,
                      bound_x_3 == "OPEN" ? Open : Periodic};
  
  if((TASK!="ENERGY" || TASK!="SINGLE_POINT") && (alpha!=90.0 || beta!=90.0 || gamma!=90.0))
  {
      cerr << "Non-orthogonal unit cells only supported for energy calculations. " << endl;
      cerr << '\'' << TASK << '\'' << endl;
      //exit(1);
  }
  alpha *= M_PI /180.0;
  beta *= M_PI / 180.0;
  gamma *= M_PI / 180.0;
  ATOMS_FILE_ADDRESS = "inputs/" + ATOMS_FILE_ADDRESS;

  // -----------------------------------------------------------
  // Decide what to return and return it.
  // -----------------------------------------------------------
  
  // If the restart flag is on
  if (RESTART=="ON"){

#ifdef DBG
    cout << "The restart flag is on." << endl;
    cout << "Reading in the unitcell from outputs/seed_unitcell.out" << endl;
#endif

    unitcell cell = readunitcell(("outputs/" + seed +"_unitcell.out").c_str());

#ifdef DBG
    cout << "Reading in positions from outputs/seed_positions.end" << endl;
#endif
    
    atom_list atoms = readpositions(("outputs/" + seed +"_positions.end").c_str(),cell);
    // Check that the positions are OK
    if (atoms.check_dist(cell)==false){
      exit(1);
    }
    sim_ptr =  new simulation(seed,cell,ENERGY_OUT_STEPS,VELOCITY_OUT_STEPS,POSITIONS_OUT_STEPS,TIME_STEP,MD_STEPS,RELAX_STEPS,SKIN_DEPTH,TEMPERATURE,TEMP_TOL,STRUCT_RELAX,RELAXER,FORCE_TOL,THERMOSTAT,INTEGRATOR,SUPERCELL,S_x,S_y,S_z,RESTART,DMD_GAMMA_RECALC_STEPS,SD_STEP_SIZE,DMD_TIMESTEP_SCALE,TASK,phonon_q_points,PHONON_DISPLACEMENT,RDF_BINS,RDF_CUTOFF,RDF_STEPS);
    // Read in the velocities.
#ifdef DBG
    cout << "Reading in velocities from outputs/seed_velocities.end" << endl;
#endif
    vec_list velocities = readvelocities(("outputs/" + seed +"_velocities.end").c_str());
    for (i=0;i<atoms.count();i++){
      atoms[i].vel = velocities[i];
    }
    return atoms;
  }
  // if supercell is on
  else if (SUPERCELL=="Y"){
#ifdef DBG
    cout << "The supercell option is turned on." << endl;
#endif
    unitcell cell(x_1,x_2,x_3,alpha,beta,gamma,boundary);
    atom_list atoms = readpositions(ATOMS_FILE_ADDRESS,cell);
#ifdef DBG
    cout << "Atoms have been read in. Now scaling them into the supercell." << endl;
#endif
    atoms.scale_supercell(S_x,S_y,S_z,cell);


    if (atoms.check_dist(cell)==false){
      exit(1);
    }
    // create the simulation object, after modifying the cell
    sim_ptr =  new simulation(seed,cell,ENERGY_OUT_STEPS,VELOCITY_OUT_STEPS,POSITIONS_OUT_STEPS,TIME_STEP,MD_STEPS,RELAX_STEPS,SKIN_DEPTH,TEMPERATURE,TEMP_TOL,STRUCT_RELAX,RELAXER,FORCE_TOL,THERMOSTAT,INTEGRATOR,SUPERCELL,S_x,S_y,S_z,RESTART,DMD_GAMMA_RECALC_STEPS,SD_STEP_SIZE,DMD_TIMESTEP_SCALE,TASK,phonon_q_points,PHONON_DISPLACEMENT,RDF_BINS,RDF_CUTOFF,RDF_STEPS);
    return atoms;
  }
  else{
       S_x = 1; S_y = 1; S_z = 1;
  // if the supercell is off and restart is off, define the simulation object
    unitcell cell(x_1,x_2,x_3,alpha,beta,gamma,boundary);
    sim_ptr =  new simulation(seed,cell,ENERGY_OUT_STEPS,VELOCITY_OUT_STEPS,POSITIONS_OUT_STEPS,TIME_STEP,MD_STEPS,RELAX_STEPS,SKIN_DEPTH,TEMPERATURE,TEMP_TOL,STRUCT_RELAX,RELAXER,FORCE_TOL,THERMOSTAT,INTEGRATOR,SUPERCELL,S_x,S_y,S_z,RESTART,DMD_GAMMA_RECALC_STEPS,SD_STEP_SIZE,DMD_TIMESTEP_SCALE,TASK,phonon_q_points,PHONON_DISPLACEMENT,RDF_BINS,RDF_CUTOFF,RDF_STEPS);
  // to get the sim object elsewhere dereference by simulation sim = *sim_ptr;
  // if supercell turned on supercell the atoms
  // otherwise return without supercell
#ifdef DBG
    cout << "The supercell option is turned off. The restart option is turned off. Reading in positions from file: " << ATOMS_FILE_ADDRESS << endl;
#endif
    return readpositions(ATOMS_FILE_ADDRESS,sim_ptr->cell);
  }
}

 // -----------------------------------------------------------
 // Start reading in the positions file
 // -----------------------------------------------------------
atom_list readpositions(string position_file, const unitcell& cell)
{
  cout << position_file << endl;
  ifstream in;
  in.open(position_file.c_str());
  if (!in){//A test to see if the file was opened
    cerr << "ERROR: The positions file could not be opened. \n";
    exit(1);
  }
  string line, letter;
  // re-initialise skipline and i.
  int i=0;
  int N;
  double posx, posy, posz;
// Get the number of atoms in the file (i.e. the first thing in the file).
  in >> N;
// Create the new list of atoms.
  atom_list atoms(N);
  for (i=0; i<N; i++){
	  in >> line; //gets rid of C
	  in >> posx;
	  in >> posy;
	  in >> posz;
      atoms[i].pos.x()=posx;
      atoms[i].pos.y()=posy;
      atoms[i].pos.z()=posz;
  }
  in.close();

    cout << "Positions file '" << position_file << "' successfully read." <<endl;

    if (!atoms.check_dist(cell)){
        exit(1);
    }
    return atoms;
}


 // -----------------------------------------------------------
 // Start reading in the velocities file
 // -----------------------------------------------------------
vec_list readvelocities(string velocity_file)
{
  cout << velocity_file << endl;
  ifstream in;
  in.open(velocity_file.c_str());
  if (!in){//A test to see if the file was opened
    cerr << "ERROR: The velocities file could not be opened. \n";
    exit(1);
  }
  string line, letter;
  int i=0;
  int N;
  double velx, vely, velz;
// Get the number of atoms in the file (i.e. the first thing in the file).
  in >> N;
// Create the new list of atoms.
  vec_list velocities(N);
  for (i=0; i<N; i++){
	  in >> line; //gets rid of C
	  in >> velx;
	  in >> vely;
	  in >> velz;
      velocities[i].x()=velx;
      velocities[i].y()=vely;
      velocities[i].z()=velz;
  }
  in.close();
    
    return velocities;
}

unitcell readunitcell(string filename){
  ifstream ufile;
  ufile.open(filename.c_str());
  double x_1, x_2, x_3, alpha, beta, gamma;
  int bound_x_1, bound_x_2, bound_x_3;
  int i=0, IMAX=100;
  string name, NAME;

if(!ufile){//if input file not opened, error message
    cerr << "ERROR: " << filename << " could not be opened for restart.\n";
    exit (1);
  }


  while (!ufile.eof())
  {
    i++;
    if (i>IMAX){
#ifdef VERBOSE
     cerr << "The while loop has been broken by IMAX." << endl;
#endif
    }
    ufile >> name;
    NAME = to_upper(name);
    if (NAME == "LENGTH_X_1"){
      ufile >> x_1;
    }
    else if (NAME == "LENGTH_X_2"){
      ufile >> x_2;
    }
    else if (NAME == "LENGTH_X_3"){
      ufile >> x_3;
    }

    else if (NAME == "ALPHA"){
      ufile >> alpha;
    }
    else if (NAME == "BETA"){
      ufile >> beta;
    }
    else if (NAME == "GAMMA"){
      ufile >> gamma;
    }

    else if (NAME == "BOUND_X_1"){
      ufile >> bound_x_1;
    }
    else if (NAME == "BOUND_X_2"){
      ufile >> bound_x_2;
    }
    else if (NAME == "BOUND_X_3"){
      ufile >> bound_x_3;
    }
    else
    {
         cerr << "Problem reading the input file. Problem occured directly after reading the word " << name << endl;
     exit(1);
    }
  }

  // Check to make sure that the boundaries are reasonable.
  if (abs(bound_x_1)>1 || abs(bound_x_2)>1 || abs(bound_x_3)>1){
    cerr << "ERROR: error occured while reading in the boundary conditions from file " << filename << endl;
    exit(1);
  }
  // create the boundary, 0 -> Open; 1-> Periodic
  bctype boundary[3]={bound_x_1 == 0 ? Open : Periodic,
                      bound_x_2 == 0 ? Open : Periodic,
                      bound_x_3 == 0 ? Open : Periodic};

  // convert into radians
  alpha *= M_PI /180.0;
  beta *= M_PI / 180.0;
  gamma *= M_PI / 180.0;
  // make unitcell
  unitcell cell(x_1,x_2,x_3,alpha,beta,gamma,boundary);
  return cell;
}

// Write out the function "to_upper" explicitly
string to_upper(string str)
{
  std::transform(str.begin(), str.end(),str.begin(), ::toupper);
  return str;
}

#endif //MODULEA_HPP
