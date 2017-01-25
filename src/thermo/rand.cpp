#ifndef THERMO_ANDERSEN_CPP
#define THERMO_ANDERSEN_CPP

#include <cmath>
#include <time.h>


namespace moldycat
{
//Fast, relatively high quality PRNG taken from Numerical Recipes 3rd edition
//Slight modification to allow multiple instances to maintain their own state 
//in order to use for multithreaded application
//If no seed is specified, will use current system time
struct ranq1 {
    unsigned long long v;

    ranq1(unsigned long long j = (unsigned) time(0)) : v(4101842887655102017LL)
    {
        v ^= j;
        v = int64(); 
    }


    inline unsigned long long int64() 
    {
        v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
        return v * 2685821657736338717LL;
    }
    inline double next() { return 5.42101086242752217E-20 * int64(); } 
};


struct  rangauss : ranq1
{
    rangauss(double mmu, double ssig, unsigned long long seed = (unsigned) time(0)) : ranq1(seed), mu(mmu), sig(ssig), storedval(0.0) {}



    double norm() 
    {
     double v1,v2,rsq,fac;
        if (storedval == 0.0) 
        {
            do 
            {
                v1=2.0*next()-1.0;
                v2=2.0*next()-1.0;
                rsq=v1*v1+v2*v2;
            } 
            while (rsq >= 1.0 || rsq == 0.0); 

           fac=sqrt(-2.0*log(rsq)/rsq);
           storedval = v1*fac;
           return mu + sig*v2*fac;
        } else 
        {
           fac = storedval;
           storedval = 0.;
           return mu + sig*fac;
        } 
   }

    private:    
    double mu,sig;
    double storedval;

};
}

#endif //THERMO_ANDERSEN_CPP

