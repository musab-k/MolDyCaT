#ifndef MOLDYCAT_ATOM_VECT3_H
#define MOLDYCAT_ATOM_VECT3_H

#include <vector>
#include <cmath>
#include <vector>
#include "armadillo"
#include <iostream>


using namespace std;
using namespace arma;

namespace moldycat
{
    class unitcell;

	class vect3 : public vec::fixed<3>
    {
        public:
            vect3();
			vect3(const vec &v);
            vect3( double x, double y, double z);

			 double&  x();
			 double&  y();
   			 double&  z();

             const double  x() const;
			 const double  y() const;
   			 const double  z() const;
            
            vect3 diff_bc(const vect3& v,const unitcell& cell) const;
            double dist_bc(const vect3& v,const unitcell& cell) const;
            void add_bc(const vect3& v,const unitcell& cell);
 


            double normsq() const;           
            double norm() const;
            double norm_bc(const unitcell& cell) const;
            double dot(const vect3 v,const unitcell& cell) const;
 
            vect3 getUnit() const;
            vect3& normalise();
            vect3 cross(const vect3 v) const;
            double dot(const vect3 v) const;
            double triple_prod(const vect3 a, const vect3 b) const;            
           friend  ostream& operator<< (ostream &out,const vect3& v);


    };

    class vec_list : public vector<vect3>
    {    
        public:
			vec_list( int count);
            const int count() const;
            void fill(double d);
            double max_norm() const;
            double dot(const vec_list& v) const;
            double normsq() const;
    };
}
#endif //MOLDYCAT_ATOM_VECT3_HPP

