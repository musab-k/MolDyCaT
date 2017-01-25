#!/usr/bin/python
import sys
import math
import string

def gen_pos(disp,cell_x=3.55648,cell_y=3.55648,cell_z=3.55648):
  #filename = "phonons/dia_gamma_" + str(disp) + ".xyz"
  filename = "phonons_gamma.txt"
  out_f = open(filename,'w')
  #a = cell_x
  disp = float(disp)
  cell_x = float(cell_x)
  cell_y = float(cell_y)
  cell_z = float(cell_z)
  atom_rel_pos = [[0,0,0],[1.0/4,1.0/4,1.0/4],[0,2.0/4,2.0/4],[1.0/4,3.0/4,3.0/4],[2.0/4,0,2.0/4],[3.0/4,1.0/4,3.0/4],[2.0/4,2.0/4,0],[3.0/4,3.0/4,1.0/4]]
  atom_abs_pos = []
  for count in range(len(atom_rel_pos)):
    if (count%2 == 0):
      print atom_rel_pos[count][0]
      atom_abs_pos.append([(atom_rel_pos[count][0]-disp)*cell_x,(atom_rel_pos[count][1]-disp)*cell_y,(atom_rel_pos[count][2]-disp)*cell_z])
    else:
      atom_abs_pos.append([(atom_rel_pos[count][0]+disp)*cell_x,(atom_rel_pos[count][1]+disp)*cell_y,(atom_rel_pos[count][2]+disp)*cell_z])

  out_f.write("8\n\n")
  for atom in range(len(atom_abs_pos)):
    out_f.write("C " + str(atom_abs_pos[atom][0]) + " " + str(atom_abs_pos[atom][1]) + " " + str(atom_abs_pos[atom][2]) + "\n")
  
  return

if __name__ == "__main__":
  gen_pos(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
