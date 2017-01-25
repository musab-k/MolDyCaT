#!/usr/bin/python
import sys
import math
import string

def gen_pos(disp,cell_x=4.257,cell_y=2.4577801,cell_z=10):
  cell_y = cell_x*math.sqrt(3)/3
  filename = "phonons_E2gb.txt"
  out_f = open(filename,'w')
  #a = cell_x
  disp = float(disp)
  atom_rel_pos = [[0,1.0/2,0],[1.0/6,0,0],[1.0/2,0,0],[2.0/3,1.0/2,0]]
  atom_abs_pos = []
  for count in range(len(atom_rel_pos)):
    if ((count == 0) | (count == 2)):
      atom_abs_pos.append([atom_rel_pos[count][0]*cell_x,(atom_rel_pos[count][1]+disp)*cell_y,atom_rel_pos[count][2]*cell_z])
    elif ((count == 1) | (count == 3)): 
      atom_abs_pos.append([atom_rel_pos[count][0]*cell_x,(atom_rel_pos[count][1]-disp)*cell_y,atom_rel_pos[count][2]*cell_z])

  out_f.write("4\n\n")
  for atom in range(len(atom_abs_pos)):
    out_f.write("C " + str(atom_abs_pos[atom][0]) + " " + str(atom_abs_pos[atom][1]) + " " + str(atom_abs_pos[atom][2]) + "\n")
  
  return

if __name__ == "__main__":
  gen_pos(sys.argv[1])
