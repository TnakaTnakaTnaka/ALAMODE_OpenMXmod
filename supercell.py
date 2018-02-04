#!/usr/bin/env python
#
# supercell.py
#
# Simple script to generate super cell data file.
# OpenMX is supported.
#
# Copyright (c) 2018 Yuto Tanaka
#

""" 
--- How to use ---
If you create (l*m*n) supercell date file, please type the command.
(l,m,n should be integer.)

$ python supercell.py --file=file.dat --supercell=lmn

"""

import numpy as np
import math
import optparse

usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage=usage)
parser.add_option('--file', help="OpenMX input file (*.dat)")
parser.add_option('--supercell', help="l*m*n supercell (lmn)")

def get_prim_data(file_in):
    search_target1 = "Atoms.Number"
    search_target2 = "<Atoms.SpeciesAndCoordinates"
    search_target3 = "Atoms.SpeciesAndCoordinates.Unit"
    search_target4 = "<Atoms.UnitVectors"

    str_frac = ["Frac", "frac", "Frac"]

    fin = open(file_in, 'r')

    #set initial patameters
    lavec_flag = 0
    lavec_row = 0
    lavec = np.zeros([3, 3])
    nat = 0
    coord_flag = 0
    coord_row = 0


    #read oroginal file and pull out some infomations
    for line in fin:
        ss = line.strip().split()
        #number of atoms
        if len(ss) > 0 and ss[0] == search_target1:
            nat = int(ss[1])
            x_ang = np.zeros([nat, 3])
            species = []
            spin = np.zeros([nat, 2], dtype = np.float64)
 
        #coordinates_unit
        if len(ss) > 0 and ss[0] == search_target3:
            coord_unit = ss[1]
        
        #coordinates
        if coord_flag == 1:
            species.append(ss[1])
            for i in range(2):
                spin[coord_row][i] = ss[i+5]
    
            for i in range(3):
                x_ang[coord_row][i] = float(ss[i+2]) 
       
            coord_row += 1
            if coord_row == nat:
                coord_flag = 0

        else:
            pass

        #latice vector
        if lavec_flag == 1:
            for i in range(3):
                lavec[lavec_row][i] = float(ss[i])
            lavec_row += 1
            if lavec_row == 3:
                lavec_flag = 0
 
        if len(ss) > 0 and ss[0] == search_target4:
            lavec_flag = 1

        if len(ss) > 0 and ss[0] == search_target2:
            coord_flag = 1

    fin.close()

    #convert to ang
    conv = (np.linalg.inv(lavec)).T
    inv_conv = np.linalg.inv(conv)
    if coord_unit in str_frac:
        for i in range(nat):
            x_ang[i] = np.dot(inv_conv, x_ang[i])

    #errors
    if nat == 0:
        print "Could not read dat file properly." 
        exit(1)

    return nat, x_ang, lavec, species, spin


def create_supercell(m, file_in, file_out):
    
    f_out = open(file_out, 'w')
    Bohr_to_angstrom = 0.5291772108    
    atom = 0
    
    nat, x_ang, lavec, species, spin = get_prim_data(file_in)
    nat_sc = nat * m[0] * m[1] * m[2]
    lavec_sc = (m * lavec.T).T
    conv = (np.linalg.inv(lavec_sc)).T
    x_frac = np.zeros([nat_sc, 3])

    for i in range(m[0]):
        for j in range(m[1]):
            for k in range(m[2]):
                tran_vec = np.zeros([3])
                m_sc = np.array([i, j, k], dtype = np.float64)
                for l in range(3):
                    tran_vec += (m_sc * lavec.T).T[l]
                for l in range(nat):
                    x_frac[atom + l] = x_ang[l] + tran_vec 
                    x_frac[atom + l] = np.dot(conv, x_frac[atom + l])
                atom += nat


    
    f_out.write("Atoms.Number       %d \n" % (nat_sc))
    f_out.write("Atoms.SpeciesAndCoordinates.Unit   frac \n")
    f_out.write("<Atoms.SpeciesAndCoordinates \n")

    for i in range(nat_sc):
        f_out.write("%3d  %s  %12.10f  %12.10f  %12.10f  %2.1f %2.1f \n" % 
                (i+1, species[int(i%nat)], 
                    x_frac[i][0], x_frac[i][1], x_frac[i][2], 
                    spin[int(i%nat)][0], spin[int(i%nat)][1]))

    f_out.write("Atoms.SpeciesAndCoordinates> \n")

    #wirte unit vector (ang)
    f_out.write("\n")
    f_out.write("Atoms.UnitVectors.Unit   Ang \n")
    f_out.write("<Atoms.UnitVectors \n")
    for i in range(3):
        f_out.write("%12.10f  %12.10f  %12.10f \n" % (lavec_sc[i][0], lavec_sc[i][1], lavec_sc[i][2]))

    f_out.write("Atoms.UnitVectors> \n")

    
    # write unit vector (bohr)
    f_out.write("\n")
    f_out.write("<alamode alm cell fields> \n")
    f_out.write("&cell \n")
    lavec_sc /= Bohr_to_angstrom
    f_out.write("  %12.10f \n" %(lavec_sc[0][0]))
    lavec_sc /= lavec_sc[0][0]
    for i in range(3):
        f_out.write("  %12.10f  %12.10f  %12.10f \n" % (lavec_sc[i][0], lavec_sc[i][1], lavec_sc[i][2]))

    f_out.write("/ \n")

    f_out.close()



def main():
    
    file_out = "coord.data"

    options, args = parser.parse_args()

    if options.file:
        file_in = options.file
    else:
        print "input file is not selected."
    
    if options.supercell:
        sc = options.supercell
    else:
        print "supercell is not defined."


    m = np.array([sc[0], sc[1], sc[2]], dtype = np.int64)
    create_supercell(m, file_in, file_out)
    print "coord.data is generated."

if __name__ == "__main__":
  main()
