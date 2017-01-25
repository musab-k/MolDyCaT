#ifndef MOLDYCAT_STEPPER_STEPPERBASE_HPP
#define MOLDYCAT_STEPPER_STEPPERBASE_HPP

#include <vector>
#include <cmath>
#include "armadillo"
#include "../atom/vect3.h"
#include "../atom/atom.h"
#include "../sim/simulation.hpp"
#include "../sim/forceeval.hpp"


using namespace std;

//-----------------------------------------------------------------------------
// Stepperbase from which all integration algorithms are derived
//
// constructor takes force evaluator object, time step size anLd number of time steps to run through
//
// This does no clever work and just makes sure that all our integrators
// follow the same conventions:
//-----------------------------------------------------------------------------
namespace moldycat
{
    class stepperbase
    {
    public:
    
        stepperbase(forceeval& forceevaluater, const simulation& sim) :maxsteps(sim.md_steps),  dt(sim.dt), mass(1244.82175), t(0.0), step_count(0) //, forcer(forceevaluater)
        {}
        double getdt() const {return dt;}


        virtual void step(atom_list& atoms,forceeval& forcer, const simulation& sim) = 0;


        const int maxsteps;
    protected:
        const double dt;
        const double mass;
        double t;
        int step_count;
     //   forceeval& forcer;
    };


}
#endif //MOLDYCAT_STEPPER_STEPPERBASE_HPP
