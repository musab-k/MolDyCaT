
#include <vector>
#include <cmath>
#include <vector>
#include "armadillo"
#include <iostream>
#include "../atom/vect3.h"
#include "unitcell.h"

using namespace std;
using namespace arma;

namespace moldycat
{
            //Initialise unit cell object from specified arguments calculating transformation matrix and inverse and storing BCs internally
            unitcell::unitcell(double A, double B, double C, double alpha, double beta, double gamma,bctype* bcs) :len(A,B,C), ang(alpha,beta,gamma)
            {

    if(alpha < 0.0 | alpha >= 180.0)
    {
        cerr << "Cell angle alpha must lie in range ]0.0,180[" << endl;
        exit(0);
    }

     if(beta < 0.0 | beta >= 180.0)
    {
        cerr << "Cell angle beta must lie in range ]0.0,180[" << endl;
        exit(0);
    }

    if(gamma < 0.0 | gamma >= 180.0)
    {
        cerr << "Cell angle gamma must lie in range ]0.0,180[" << endl;
        exit(0);
    }

            calcS();
           for(int i=0;i<3;i++)
                bc[i] = bcs[i];

            }

            void unitcell::calcS()
            {
                double ca = cos(ang[0]);
                double cb = cos(ang[1]);
                double cg = cos(ang[2]);

                double eps = 1E-15;



               S = zeros<mat>(3,3);
               S(0,0) = 1.0;
               S(0,1) = cg;
               S(0,2) = cb;

               S(1,1) = sqrt(1.0 - cg*cg);
               S(1,2) = (ca-cg*cb)/S(1,1);

               S(2,2) = sqrt(1.0-cb*cb-S(1,2)*S(1,2));

                if(abs((ang[0]/(M_PI/2.0))-1)<eps){
                    if(abs((ang[1]/(M_PI/2.0))-1)<eps){
                        if(abs((ang[2]/(M_PI/2.0))-1)<eps){
                            S= eye<mat>(3,3);
                        }
                    }
                }

               Sinv = inv(S);

            }
            
            // Transform difference of v1 and v2 into s-space
            // Apply regular cartesian BCs and transform back into real space
            //and inverse Sinv
            vect3 unitcell::diff_bc(const vect3& v1,const vect3& v2) const
            {
                vec3 dr = (v1-v2);
                vec3 ds = S*(dr);
               // cout << S << endl;
               // cout << Sinv << endl;
             //   cout << ds << endl;

                //cout << "Before trans\t" << (v1-v2) << endl;
                //cout << "After trans\t" << ds << endl;
           
                for(int i=0; i < 3; i++){
                   //if(bc[i] == Periodic & abs(ds[i]) > 0.5 * len[i]) ds[i] = len[i] - abs(ds[i]);
                   if(bc[i] == Periodic && ds[i]/len[i] > 0.5){
                       ds[i]-=len[i];
                   }
                   else if(bc[i] == Periodic && ds[i]/len[i] < -0.5){
                       ds[i]+=len[i];
                   }

                }
                 
                 //cout << "After BCs\t" <<  ds << endl;
                 //cout << "Trans back\t" << (Sinv*ds) << endl;

                //Returns everything in terms of orthorhombic coordinates
                //return vect3(Sinv*ds);
                return vect3(S*ds);
                //return vect3(Sinv*ds);       
            
           }

            double unitcell::dist_bc(const vect3& v1,const vect3& v2) const
            {
                vec3 ds = S*(v1-v2);


                //cout << "Before trans\t" << (v1-v2) << endl;
               // cout << "After trans\t" << static_cast<vect3>(ds) << endl;
           
                for(int i=0; i < 3; i++){
                   if(bc[i] == Periodic & abs(ds[i]) > 0.5 * len[i]) ds[i] = len[i] - abs(ds[i]);
                }
                 
                return static_cast<vect3>(S*ds).norm();       

                // previously
                // return vect3(Sinv*ds).norm(); 
           }


            void unitcell::add_bc(vect3& v1,const vect3& v2) const
            {
          
                vec3 ds = S*(v1+v2);
           
                for(int i=0; i < 3; i++)
                   if(bc[i] == Periodic)  ds[i] -= floor( ds[i]/len[i]) * len[i];
                
                
                v1 =  static_cast<vect3>(S*ds);       
            
           }

           double unitcell::dot(const vect3 v1,const vect3 v2) const
           {    
                
                vect3 tmp = static_cast<vect3>(S*v1);
                return tmp.dot(static_cast<vect3>(S*v2));
           }

           double unitcell::norm(const vect3 v) const
           {    
                return unitcell::dot(v,v);
           }
            

            // Calculate volume using formula from 
            // http://onlinelibrary.wiley.com/doi/10.1002/9780470727133.app1/pdf
            double unitcell::volume() const
            {
                double ca = cos(ang[0]);
                double cb = cos(ang[1]);
                double cg = cos(ang[2]);

                return len[0]*len[1]*len[2]*sqrt(1.0 - ca*ca - cb*cb - cg*cg + 2*(ca*cb*cg));
               
            }
      
}

