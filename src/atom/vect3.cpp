#include <vector>
#include <cmath>
#include <vector>
#include "armadillo"
#include <iostream>
#include "vect3.h"
#include "../sim/unitcell.h"

#include "assert.h"

using namespace std;
using namespace arma;

//-----------------------------------------------------------------------------
// vector object - essentially wrapper around armadillo vec3 object so we get
//                 to make use of all their optimised arithmetic operations
//
//
// Hope naming convection is sufficiently clear so I wont go through each function in detail
// All functions ending in _bc will be corrected for unitcells (suitable for positions 
// while standard +/- are not (suitable for velocities)
// 
// Consider vect3 A, vect3 B, vect3 C
//
//  Things like A += B * 0.5; 
//
// work as expected but things like
//
// C = A + B;
//
// need a little more love. Both
//
// C = static_cast<vect3>(A+B);
// and
// C = vect3(A+B);
//
// should work, although the first is preferable as it does not involve allocating a temporary
//
//-----------------------------------------------------------------------------


 namespace moldycat
 {
            vect3::vect3()
            {
                x() = 0.0; y() = 0.0; z() = 0.0;
            }
			vect3::vect3(const vec &v) : vec::fixed<3>(v){}

            vect3::vect3( double x, double y, double z) : vec::fixed<3>()
            {
                vec::fixed<3>::operator()(0) = x;
                vec::fixed<3>::operator()(1) = y;
                vec::fixed<3>::operator()(2) = z;
            }

			 double&  vect3::x() {return vec::fixed<3>::operator[](0);}
			 double&  vect3::y() {return vec::fixed<3>::operator[](1);}
			 double&  vect3::z() {return vec::fixed<3>::operator[](2);}

			 const double vect3::x() const {return vec::fixed<3>::operator[](0);}
			 const double  vect3::y() const {return vec::fixed<3>::operator[](1);}
			 const double  vect3::z() const {return vec::fixed<3>::operator[](2);}
            

            //Cell corrected vector operations are done in unit cell class
            vect3 vect3::diff_bc(const vect3& v,const unitcell& cell) const 
            {
                return cell.diff_bc(*this,v);             
            }
            void vect3::add_bc(const vect3& v,const unitcell& cell) 
            {
                cell.add_bc(*this,v);
            }
            double vect3::dist_bc(const vect3& v,const unitcell& cell) const 
            {
                return cell.dist_bc(*this,v); 
            
            }


            vect3 vect3::cross(const vect3 v) const
            {
                return vect3(y()*v.z() - z()*v.y(),
                             z()*v.x() - x()*v.z(),
                             x()*v.y() - y()*v.x());
            }

            double vect3::dot(const vect3 v) const
            {
                return x()*v.x()+y()*v.y()+z()*v.z();
            }

           double vect3::normsq() const { return dot(*this);}
           double vect3::norm() const {return sqrt(vect3::normsq());}

            vect3 vect3::getUnit() const 
            {
                double mag = norm();
                return vect3(x() / mag,y() / mag,z() / mag);
            }

            vect3& vect3::normalise() 
            {
                vec::fixed<3>::operator/=(norm());
                return *this;
            }

            double vect3::norm_bc(const unitcell& cell) const
            {
                return cell.norm(*this);
            }
            double vect3::dot(const vect3 v,const unitcell& cell) const
            {
                return cell.dot(*this,v);
            }


            double vect3::triple_prod(const vect3 a, const vect3 b) const
            {
                return b.dot(cross(a));
            }
            
            ostream& operator<< (ostream &out,const vect3& v) 
            {
                out << '(' << v.x() << ',' << v.y() << ',' << v.z() << ')';
                return out;
            }


//-----------------------------------------------------------------------------
// Wrapper around STL vect3 vector
//-----------------------------------------------------------------------------


			vec_list::vec_list( int count) : vector<vect3>(count)
            {
                (*this).fill(0.0);
            }
            const int vec_list::count() const
            {return vector<vect3>::size();}

            // Fill vector to chosen value
            // to make usage consistent with
            // Armadillo types
            void vec_list::fill(double d)
            {
                for(unsigned int i = 0; i < size(); i++)
                {
                    vector<vect3>::operator[](i).x() = d;
                    vector<vect3>::operator[](i).y() = d;
                    vector<vect3>::operator[](i).z() = d;
                }
            }

            // Return the maximum value of the norm
            // of all component vect3s (useful for
            // dertermining if a relaxation algorithm
            // has converged
            double vec_list::max_norm() const
            {
                double max = (*this)[0].norm();
                for(unsigned int i = 1; i < size(); i++)
                {
                    if((*this)[i].norm() > max)
                        max = (*this)[i].norm();
                }
                return max;
            }

            // dot product between two vector lists
            double vec_list::dot(const vec_list& v) const
            {
                double d;
                for(unsigned int i = 1; i < size(); i++)
                {
                    d += (*this)[i].dot(v[i]);
                }
                return d;
            }

            double vec_list::normsq() const
            {
                return (*this).dot(*this);
            }
}
