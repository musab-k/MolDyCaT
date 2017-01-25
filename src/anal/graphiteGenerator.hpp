#ifndef GRAPHITEGENERATOR_HPP
#define GRAPHITEGENERATOR_HPP

//-----------------------------------------------------------------------------------
// Function generateGraphite
//
// Params: atom list, unit cell
// Returns: nothing
// Operation:   - fills the variable atoms with the positions of N graphene atoms, 
//              where N is the size of the atom list. 
//              - NB: Make sure N is a multiple of 40 (ideally at least 160) to get 
//              good results.
//              - method will also fill the unit cell with the correct lengths to ensure
//              the PBCs work as they should.
//-----------------------------------------------------------------------------------
void generateGraphite(atom_list& atoms, unitcell& cell)
{
    int N=atoms.count(),hlim=20,klim=N/(2*hlim);
    double a,b,c, Lx, Ly, Lz, x,y,z;
    int i,k;
    a=1.42;
    b=a+a;
    c=a*sqrt(3.0);

    Lz=3.35;
    z=Lz*0.5;
    x=0.0;
    y=a*0.5;


    for (k=0 ; k<klim ; k++){
        // : : : : : 
        for(i=0 ; i<hlim-1; i+=2){
            atoms[k*40+i].pos=vect3(x,y,z);
            atoms[k*40+i+1].pos=vect3(x,y+a,z);
            x+=c;
        }

        // Go up and do : . : . : where . is the top atom from previous loop
        x=a*sqrt(3.0)*0.5;
        y+=1.5*a;
        for(i=20 ; i<hlim*2-1; i+=2){
            atoms[k*40+i].pos=vect3(x,y,z);
            atoms[k*40+i+1].pos=vect3(x,y+a,z);
            x+=c;
        }
        y+=1.5*a;
        Lx=x-c*0.5;
        x=0.0;
    }
    Ly=y-0.5*a;
    cell.len = vect3(Lx,Ly,Lz);
}
#endif
