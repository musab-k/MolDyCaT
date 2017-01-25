#ifndef STEPPER_FACTORY_HPP
#define STEPPER_FACTORY_HPP

#include "../stepper/stepperbase.hpp"
#include "../stepper/verlet.hpp"
#include <cmath>


stepperbase* integratorFactory(const simulation& sim,const atom_list& atoms,forceeval& forcer)
{
    stepperbase* integrator;
    if(sim.integrator == "VERLET")
    {
        integrator = new verlet(forcer,atoms,sim);
    }
    else if(sim.integrator == "SI2A")
    {
         double a[2] = {0.5,0.5};
         double b[2] = {0.0,1.0};
         integrator  = new sprk<2>(forcer, atoms,sim, a,b);
    }
    else if(sim.integrator == "SI3A")
    {
        double a[3] = {(2.0/3.0),
                        -(2.0/3.0),
                        1.0};
        double b[3] = {(7.0/24.0),
                        (3.0/4.0),
                        (-1.0/24.0)};
     integrator  = new sprk<3>(forcer, atoms,sim, a,b);

    }
    else if(sim.integrator == "SI3B")
    {
        double a1 = 0.9196615230173999;
        double a2 = 0.25/a1-a1/2;
        double a[3] = {a1,a2,1.0-a1-a2};                   
        double b[3] = {a[2],a[1],a[0]};

        integrator  = new sprk<3>(forcer, atoms,sim, a,b);

    }
    else if(sim.integrator == "SI4A")
    {
        double a1 = (2.0 + pow(2.0,(1.0/3.0)) + pow(2.0,(-1.0/3.0))) / 6.0;
        double a2 = (1.0 + pow(2.0,(1.0/3.0)) - pow(2.0,(-1.0/3.0))) / 6.0;
 
        double a[4] = {a1,a2,a1,a2};
        double b3 = 1.0/(1.0 - pow(2.0,(2.0/3.0)));
        double b2 = 1.0/(2.0 - pow(2.0,(1.0/3.0)));

        double b[4] = {0.0,b2,b3,b2};

        integrator  = new sprk<4>(forcer, atoms,sim, a,b);
    }
    else if(sim.integrator == "SI4B")
    {
            double a[4] = {0.515352837431122936,
                          -0.085782019412973646,
                           0.441583023616466524,
                           0.128846158365384185};
            double b[4] = {0.134496199277431089,
                          -0.224819803079420806, 
                           0.756320000515668291,
                           0.334003603286321425};

       integrator  = new sprk<4>(forcer, atoms,sim, a,b);
    }
    else if(sim.integrator == "SI4C")
    {
        double a[5] = {0.205177661542290, 
                           0.403021281604210,
                          -0.120920876338910, 
                           0.512721933192410,
                           0.000000000000000};

        double b[5] = {0.061758858135626,
                           0.338978026553640,
                           0.614791307175580,
                          -0.140548014659370,
                           0.125019822794530};
             
       integrator  = new sprk<5>(forcer, atoms,sim, a,b);
    }
    else if(sim.integrator == "SI6A")
    {
        double a1 = 0.78451361047756;
        double a2 = 0.23557321335936;
        double a3 =-1.17767998417890;
        double a4 = 1.31518632068390;
        double a[8] = {a1,a2,a3,a4,a3,a2,a1,0.0};

        double b1 = 0.39225680523878;
        double b2 = 0.51004341191846;
        double b3 =-0.47105338540976;
        double b4 = 0.06875316825252;
        double b[8] = {b1,b2,b3,b4,b4,b3,b2,b1};


       integrator  = new sprk<8>(forcer, atoms,sim, a,b);
    }
    else
    {
        cerr << "Specified integration algorithm " << sim.integrator << " invalid." << endl;
        exit(1);
    }
    return integrator;
}


#endif //STEPPER_FACTORY_HPP
