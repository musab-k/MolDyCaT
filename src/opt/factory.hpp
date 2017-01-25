#ifndef OPT_OPTFACTORY_HPP
#define OPT_OPTFACTORY_HPP

#include "optbase.hpp"
#include "dmdcm.hpp"
#include "dmdum.hpp"
#include "cg.hpp"
#include "sd.hpp"
#include <cmath>


optbase* optFactory(forceeval forcer, const atom_list& atoms, const simulation& sim)
{
    double forcetol = sim.force_tol; /// Needs to come from simulation object in future!
    optbase* relaxer;
    if(sim.relaxer == "DMD" || sim.relaxer == "DMDCM" || sim.relaxer == "DMD_COUPLED_MODES")
    {
        relaxer = new dmdcm(forcer,atoms,sim.relax_steps,forcetol,sim.dmd_gamma_recalc_steps,sim.dmd_timestep_scale);
    }
    else if(sim.relaxer == "DMDUM" || sim.relaxer == "DMD_UNCOUPLED_MODES")
    {
        relaxer = new dmdum(forcer,atoms,sim.relax_steps,forcetol,sim.dmd_gamma_recalc_steps,sim.dmd_timestep_scale);
    }
    else if(sim.relaxer == "SD" || sim.relaxer == "STEEPEST_DESCENT" || sim.relaxer == "STEEPEST_DESCENTS" || sim.relaxer == "GRADIENT_DESCENT" || sim.relaxer == "GRADIENT_DESCENTS")
    {
        relaxer = new sd(forcer,atoms,sim.relax_steps,forcetol,sim.sd_step_size);
    }
    else if(sim.relaxer == "CG" || sim.relaxer == "CONJUGATE_GRADIENT" || sim.relaxer == "CONJUGATE_GRADIENTS")
    {
         relaxer  = new cg(forcer,atoms,sim.relax_steps,forcetol,sim.sd_step_size);
    }
    else
    {
        cerr << "Invalid relaxation algorithm specified." << endl;
        exit(1);
    }
 return relaxer;
}


#endif //OPT_OPTFACTORY_HPP
