#include <iostream>
#include <cmath>
#include <complex>

#include "../atom/atom.h"
#include "../atom/vect3.h"
#include "../sim/unitcell.h"
#include "../modules/moduleE.hpp"
#include "../modules/moduleB.hpp"
#include "../sim/forceeval.hpp"
#include <omp.h>

typedef complex<double> cx;

void buildDynamicalMatrix(atom_list& atoms_orig,cx_mat& evects, cx_vec& evals, forceeval forcer_orig, vect3 q, const simulation& sim) {
  const unitcell cell = sim.cell;
  const int sx = sim.sx;
  const int sy = sim.sy;
  const int sz = sim.sz;
  const int n = atoms_orig.count()/(sx*sy*sz);
  cx_mat DyMat = zeros<cx_mat>(3*n,3*n);

  int x,y,j_index,cell_count;
  const cx ii(0.0,1.0); 
  const double h = sim.phonon_displacement;
  const double mass_const = 1244.82175;
  double force_const;

  vec_list F_orig(atoms_orig.count());
  F_orig = forcer_orig.calc_forces(atoms_orig);

    //Get number of threads we will be using and allocate arrays to hold pointers to
    //forcer, atoms and F_new objects to use in each thread
  //int thread_count = omp_get_num_threads();
/*  forceeval** forcers = new  forceeval*[thread_count];
  atom_list** atomss =  new atom_list*[thread_count];
  vec_list** F_news = new vec_list*[thread_count];

 // copy our forcer and atom_list to new thread safe versions
  for(int t = 0; t < thread_count; t++)
  {
    forcers[t] = new forceeval(atoms_orig,sim);
    //cout << &(forcer.calc_forces(atoms)) << "\t" <<  (forcers[t]->calc_forces(atoms))[0] << endl;
    //forcers[t]->calc_forces(atoms);
    atomss[t] = new atom_list(atoms_orig);
    F_news[t] = new vec_list(atoms_orig.count());
  }*/



#pragma omp parallel private(x,y,j_index,cell_count,force_const) num_threads(4)
{
forceeval* forcer_ptr = new forceeval(atoms_orig,sim);
atom_list* atoms_ptr = new atom_list(atoms_orig);
vec_list* F_new_ptr = new vec_list(atoms_orig.count());

  forceeval& forcer = *forcer_ptr;
  atom_list& atoms = *atoms_ptr;
  vec_list& F_new = *F_new_ptr;



#pragma omp for 
for (int j=0;j<n;j++) {       
      
    //figure out which thread we are in (id) 
    // and replace each call to forcer with forcers[id]->
    // and each call to atoms with (*atoms[id])
    // and each call to F_new with (*F_new[id]) 
     int id = omp_get_thread_num();
  
//       vec_list F_newn = (*forcers[id]).calc_forces(atoms);
      
    for (int j_comp=0;j_comp<3;j_comp++) { 
      y = 3*j + j_comp;
      //cout << "Row " << y+1 << " out of " << n*3 << endl;
      //Sum over atom j in all unit cells
      cell_count = 0;
      for (int ux=0;ux<sx;ux++) {
        for (int uy=0;uy<sy;uy++) {
          for (int uz=0;uz<sz;uz++) {
            //Determine index of atom j in unit cell (ux,uy,uz) in expanded atom list after super cell scaling
            
            //j_index = (4*ux + 2*uy + uz)*n + j;
            //cout << cell_count << endl;
            j_index = cell_count*n + j;
            cell_count += 1;
            //cout << "--Unit cell " << cell_count << " out of " << sx*sy*sz << endl;
            //cout << "j_index = " << j_index << ", pos=" << atoms[j_index].pos << endl;
            //cout << "j=" << j << ", j_comp=" << j_comp << ", ux=" << ux << ", uy=" << uy << ", uz=" << uz << endl;
            //Displace that atom in j_comp direction
            double eq_disp = atoms[j_index].pos[j_comp];
            atoms[j_index].pos[j_comp] += h;
            cout << "id : " << id <<  "Displaced pos = " << atoms[j_index].pos << endl;
            //Calculate force derivative on the i_comp of atom i in the (0,0,0) unit cell
            F_new = forcer.calc_phonon_forces(atoms);
            for (int i=0;i<n;i++) {
              //cout << "Orig F = " << F_orig[i] << ", New F = " << F_new[i] << endl;
              for (int i_comp=0;i_comp<3;i_comp++) {
                x = 3*i + i_comp;
                //cout << x << "\t" << y << endl;
                //cout << "i=" << i << ", i_comp=" << i_comp << endl;
                force_const = -(F_new[i][i_comp] - F_orig[i][i_comp])/h;
                cout << "force const " << force_const << endl;
                //Return atom
                atoms[j_index].pos[j_comp] = eq_disp;
                //Displacement from (m,j) to (0,i) - NEED TO CHECK!!
                vect3 R_diff = atoms[i].pos.diff_bc(atoms[j_index].pos,cell);
                //May need 2PI factor?
                cx exp_factor = exp(-ii*q.dot(R_diff));
                //cout << "DyMat(" << x << "," << y <<") += " << force_const*exp_factor << endl;
                //cout << atoms[i].pos << endl;
                //cout << atoms[j_index].pos << endl;
                //cout << R_diff << endl;
                DyMat(x,y) += force_const*exp_factor/mass_const;

                //cout << "-----------------" << endl;
              }
            }
            atoms[j_index].pos[j_comp] = eq_disp;
          }
        }
      }
    }
  }
  delete forcer_ptr;
  delete atoms_ptr;
  delete F_new_ptr;

}

/*  for(int t = 0; t < thread_count; t++)
  {
    delete forcers[t];
    delete atomss[t];
    delete F_news[t];
  }
  delete[] forcers;
  delete[] atomss;
  delete[] F_news;*/
  eig_gen(evals,evects,DyMat);
  evals.print("Evals");
  //evects.print("Evects");
  cx_mat evects_t = strans(evects);
  evects_t.save("outputs/evects.dat",raw_ascii);
  evals.save("Evals.dat",raw_ascii);
  cout << "VC TESt" << endl;
    DyMat.print("DyMat");

  return;
}

