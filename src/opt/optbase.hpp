#ifndef MOLDYCAT_OPT_OPTBASE_HPP
#define MOLDYCAT_OPT_OPTBASE_HPP

#include <vector>
#include <cmath>
#include "armadillo"
#include "../atom/vect3.h"
#include "../atom/atom.h"
#include "../sim/simulation.hpp"
#include "../sim/forceeval.hpp"


using namespace std;

//-----------------------------------------------------------------------------
// Optbase from which all relaxation algorithms are derived
//
// This does no clever work and just makes sure that all our relaxation algorithms
// follow the same pattern
// step() will continue to return true until the relaxer is converged to within the
// tolerance, the maximum number of iterations is reached or an error prevents
// the next step from being calculated.
// To test whether the stop condition was caused by convergence or not
// optbase provides is is_converged() method
//-----------------------------------------------------------------------------
namespace moldycat
{
    class optbase
    {
    public:
    
        optbase(forceeval& forceevaluater, int max_steps,double force_tolerance ) : 
            maxsteps(max_steps), mass(1244.82175),step_count(0), converged(false), force_tol(force_tolerance)//, forcer(forceevaluater)
        {}

        virtual bool step(atom_list& atoms,forceeval& forcer) = 0;


        const int maxsteps;
        const double mass;

         bool is_converged() const
        {
            return converged;
        }

        double get_max_force() const
        {
            return max_force;
        }
    protected:
        int step_count;
        bool converged;
        const double force_tol;
        double max_force;
    };


}
#endif //MOLDYCAT_OPT_OPTBASE_HPP
