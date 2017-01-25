#ifndef THERMO_THERMOFACTORY_HPP
#define THERMO_THERMOFACTORY_HPP

#include "../thermo/thermobase.cpp"
#include "../thermo/velocityscaler.cpp"
#include "../thermo/andersen.cpp"
#include "../thermo/nothermo.cpp"
#include <cmath>


thermobase* thermoFactory(const simulation& sim)
{
    thermobase* thermostat;
    if(sim.thermostat == "VS" || sim.thermostat == "VELOCITY_SCALING")
    {
        thermostat = new velocityscaler(sim.temp,sim.temptol);
    }
    else if(sim.thermostat == "ANDERSEN" || sim.thermostat == "A")
    {
         thermostat  = new andersen(sim.temp,sim.temptol,sim.dt);
    }
    else if(sim.thermostat == "NOTHERMO" || sim.thermostat == "NO_THERMOSTAT")
    {
         thermostat = new nothermo();
    }
    else
    {
        cerr << "Thermostat specified is not valid" << endl;
        exit(1);
    }
 return thermostat;
}


#endif //THERMO_THERMOFACTORY_HPP
