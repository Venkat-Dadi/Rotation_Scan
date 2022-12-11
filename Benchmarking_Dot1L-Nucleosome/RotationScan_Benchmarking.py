#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 23:51:39 2022

@author: Venkatareddy
"""

'''
Usage:
$ python ProteinEuler_v3.1.py <parameters.ini>

eg. $python ProteinEuler_v3.py parameters.inin

'''

import numpy as np
import re
from scipy.spatial.transform import Rotation as R
from scipy.spatial.distance import cdist
from datetime import date, time, datetime
import datetime
import time as time
import pandas as pd
import sys
import os
import shutil
from tqdm import tqdm
import argparse
import configparser

exe_time_start = datetime.datetime.now()

YEAR        = datetime.date.today().year
MONTH       = datetime.date.today().month
DATE        = datetime.date.today().day
HOUR        = datetime.datetime.now().hour
MINUTE      = datetime.datetime.now().minute
SECONDS     = datetime.datetime.now().second

if len(sys.argv) == 1:
        print("\nProteinEuler")
        print('''\nPlease specify the input file from the command line!
Example: python ProteinEuler.py parameter_file.ini

Example parameter_file.ini:
    
[INPUT_FILES]
EnzymePDB = <file name>    # input enzyme pdb file (.pdb)reaction center at
				  # origin and molecule is aligned to desired
                                 # position (eg. molecule's centroid on y-axis ), 
                                 # if required.
SubstratePDB = <file name>  # input substrate pdb file (.pdb) reaction
                                     # center at origin and molecule is aligned 
                                     # to desired position (eg. molecule's 
                                     # centroid on y-axis ), if required.
                                                                
[SUFFIX]
EnzymeSuffix = <string>    # name of an enzyme used as suffix for output (eg. dot1L)
SubstrateSuffix = <string>  # name of the substrate used as suffix for file output (eg. NCP)

[ROTATION_PARA]
Resolution = <integer>   # Resolution of a rotation scan. High resolution scans 
                       # are computationally expensive
                       
# Trait-Bryan angles for roation in 3D                          
AlphaLower = <integer>    # lower value of an rotation angle range about x-axis
AlphaUpper = <integer>    # upper value of an rotation angle range about x-axis  
BetaLower = <integer>     # lower value of an rotation angle range about y-axis
BetaUpper = <integer>     # upper value of an rotation angle range about y-axis
GammaLower = <integer>    # lower value of an rotation angle range about z-axis
GammaUpper = <integer>    # upper value of an rotation angle range about z-axis
''')
        
       # print
        sys.exit()
        
pars = sys.argv[1]
config = configparser.ConfigParser()
config.read(pars)

EnzymePDB = config['INPUT_FILES']['EnzymePDB'].split('#')[0].strip()
SubstratePDB = config['INPUT_FILES']['SubstratePDB'].split('#')[0].strip()
EnzymeSuffix = config['SUFFIX']['EnzymeSuffix'].split('#')[0].strip()
SubstrateSuffix = config['SUFFIX']['SubstrateSuffix'].split('#')[0].strip()
Resolution = int(config['ROTATION_PARA']['Resolution'].split('#')[0].strip())
AlphaLower = int(config['ROTATION_PARA']['AlphaLower'].split('#')[0].strip())
AlphaUpper = int(config['ROTATION_PARA']['AlphaUpper'].split('#')[0].strip())
BetaLower = int(config['ROTATION_PARA']['BetaLower'].split('#')[0].strip())
BetaUpper = int(config['ROTATION_PARA']['BetaUpper'].split('#')[0].strip())
GammaLower = int(config['ROTATION_PARA']['GammaLower'].split('#')[0].strip())
GammaUpper = int(config['ROTATION_PARA']['GammaUpper'].split('#')[0].strip())



ENZ_file = open(EnzymePDB, 'r') # Ligand pdb file
ENZ_prefix = EnzymeSuffix # Ligand pdb file prefix

SUB_file = open(SubstratePDB, 'r')
SUB_prefix = SubstrateSuffix

'''
vdw radii for atoms taken from Probe program
 ref: Probe program; DOI: 10.1006/jmbi.1998.2400 
atom_vdw = {'H': -2.0, 'C': 1.7, 'N': 1.625, 'O': 1.480, 'P': 1.871, \
            'S': 1.782} # parm99
'''
atom_vdw = {'H': -2.0, 'C': 1.75, 'N': 1.55, 'O': 1.40, 'P': 1.80, 'S': 1.80}


''' *xyz_vdw* function extracts the atomic coordinates from a pdb file and 
assigns vdW radius to each atom based on ELEMENT record '''

def xyz_vdw(pdb):
    atoms_list = []
    pdb_xyzvdw = []
    for b in pdb.readlines():
        if b.startswith(('ATOM', 'HETATM')): # and b[13:16].strip(' ') != 'DUM':
            with open(new_job_path+'/'+str(pdb.name), 'w') as f:
                f.write(b+'\n')
            for key in atom_vdw:
                if b[77:78] == key:
                    xyz_vdw_e = np.fromstring(b[31:54]+' '+\
                    b[61:66].replace(b[61:66], str(atom_vdw[key])), sep = ' ')
                    pdb_xyzvdw.append(xyz_vdw_e)
    xyz_vdw_array = np.array(pdb_xyzvdw)                 
    return xyz_vdw_array     

''' *ProteinEuler* function rotats coordinates of one pdb (ENZ) with respect to
the other (SUB) and calculates the number of clashes between each rotation.
Input parameters are:
    ENZ_xyz = input pdb coordinates in the form of ndarray (nX3); mobile pdb
    SUB_xyz = input pdb coordinates in the form of ndarray (nX3); fixed pdb
    alpha = Tait-Bryan angle, rotation about 

'''
def ProteinEuler(ENZ_xyz, SUB_xyz, alpha, beta, gamma):
    ''' rotation is extrinsic i.e. axes attaced to the object do not move.
     global frame of reference '''
    r = R.from_euler('zyx', [gamma, beta, alpha], degrees=True ) # generates euler rotation zyx extrinsic 
    Rot_matrix = r.as_matrix() # Euler rotation as matrix
    ENZ_xyz_T = ENZ_xyz.T # transpose of ENZ_xyz (enz cooordinates)
    MM = np.matmul(Rot_matrix, ENZ_xyz_T) # multiplication of xyz coord with rotation matrix
    ENZ_xyz_MM = MM.T # transpose of rotated coord 
    pair_wise_dist = cdist(SUB_xyz, ENZ_xyz_MM, 'euclidean') # calculates pairwise distances between rotated ligand coord and receptor coord
    ##################################################################
    clash_matrix = np.less_equal(pair_wise_dist, vdw_sum_limit_04).astype(dtype=int) # element-wise comparison of two ndarrays
    ##################################################################
    num_clashes_ENZ = np.sum(clash_matrix, axis = 0) # num. of clashed with enzyme involved array
    num_clashes_SUB = np.sum(clash_matrix, axis = 1) # num. of clashed with substrate involved array
    total_num_clashes = np.sum(clash_matrix)  # total number of clashes 
    num_clash_atoms_ENZ_array = np.greater_equal(num_clashes_ENZ, 1).astype(dtype=int)
    num_clash_atoms_ENZ = np.sum(num_clash_atoms_ENZ_array)
    num_clash_atoms_SUB_array = np.greater_equal(num_clashes_SUB,1).astype(dtype=int)
    num_clash_atoms_SUB = np.sum(num_clash_atoms_SUB_array)
    total_num_clash_atoms = num_clash_atoms_ENZ + num_clash_atoms_SUB
    clash_data = [alpha, beta, gamma, num_clash_atoms_ENZ, num_clash_atoms_SUB, total_num_clash_atoms, total_num_clashes]
    #np.savez_compressed(new_npz_path+'/MT67_{}_{}_{}'.format(alpha, beta, gamma), ENZ_xyz_MM.round(3))
    return clash_data    

###############################################################################
###########          SET SCAN RANGE and RESOLUTION          ###################
###############################################################################
'''
Tait-Bryan angles used for rotating ENZ pdb.
alpha = angle of rotation about x-axis
beta = angle of rotation about y-axis
gamma = angle of rotation about z-axis
res = angle interval 
angles are in degrees
'''
Resolution = int(Resolution)
alpha_range = np.arange(AlphaLower, AlphaUpper+Resolution, Resolution)  # rotation about z axis in zyx extrinsic convention
beta_range = np.arange(BetaLower, BetaUpper+Resolution, Resolution)  # rotation about y axis in zyx extrinsic convention
gamma_range = np.arange(GammaLower, GammaUpper+Resolution, Resolution)  # rotation about x axis in zyx extrinsic convention

###############################################################################

startTime   = time.time()


new_job_dir = 'results_PE_{}_{}_res{}_{}{}{}_{}h{}m{}s'.format(ENZ_prefix, SUB_prefix, str(Resolution).zfill(3), YEAR, str(MONTH).zfill(2), str(DATE).zfill(2), HOUR, MINUTE, SECONDS)
#new_npz_dir = 'ENZ_coord_npz_{}{}{}_{}h{}m{}s'.format(YEAR, str(MONTH).zfill(2), str(DATE).zfill(2), HOUR, MINUTE, SECONDS)

pwd_path = os.getcwd()
new_job_path = os.path.join(pwd_path, new_job_dir)
#new_npz_path = os.path.join(new_job_path, new_npz_dir)
os.mkdir(new_job_path)
#os.mkdir(new_npz_path)
#output = open(new_job_path+'/results.data', 'w')

           
# Prepared enzyme cooridnates from pdb using *pdb_xyz* function
ENZ_xyz_vdw = xyz_vdw(ENZ_file)    
ENZ_xyz = ENZ_xyz_vdw[:, :3]  # enz coordinates
ENZ_vdw = ENZ_xyz_vdw[:, -1]  # enz vdW radii

# Prepared substrate cooridnates from pdb using *pdb_xyz* function
SUB_xyz_vdw = xyz_vdw(SUB_file)
SUB_xyz = SUB_xyz_vdw[:, :3]  # sub coordinates
SUB_vdw = SUB_xyz_vdw[:, -1]  # sub vdW radii
'''
Prepares two matrices with grid point elements of each other.
ENZ_vdw_m = a matrix with rows as ENZ atoms' vdw radii with number of rows 
            equal to length of SUB atoms.
            # of columns = # of atoms in ENZ.
            # of rows = # of atoms in SUB.
SUB_vdw_m = a matrix with columns as SUB atoms' vdw radii with number of 
            columns is equal to length of ENZ atoms. 
            # of columns = # of atoms in ENZ.
            # of rows = # of atoms in SUB.            
size of ENZ_vdw_m = size of SUB_vdw_m            
'''
ENZ_vdw_m, SUB_vdw_m = np.meshgrid(ENZ_vdw, SUB_vdw) 
vdw_sum = np.add(ENZ_vdw_m, SUB_vdw_m) # sum of two matrices 
clash_threshold = 0.4 # input clash threshold value
vdw_sum_limit_04 = vdw_sum - clash_threshold   # clash threshold matrix
#vdw_sum_limit_avr = vdw_sum*0.5

pbar = tqdm(total = len(alpha_range)*len(beta_range)*len(gamma_range))
clash_data = []
for a in alpha_range:
    for b in beta_range:
        for g in gamma_range:
            clash_data.append(ProteinEuler(ENZ_xyz, SUB_xyz, a, b, g))
            pbar.update(1)
pbar.close()            
            
# =============================================================================
# clash_data_array = np.array(clash_data)            
# header = ['alpha', 'beta', 'gamma', 'ENZ_CA', 'SUB_CA', 'TOT_CA', 'clashes'] 
# clash_data_array_df = pd.DataFrame(clash_data_array, columns=header)
# sorted_val_df = clash_data_array_df.sort_values(by = ['TOT_CA'], axis=0, ascending=True)
# output_filename = new_job_path+'/output_ProteinEuler_{}_{}_res{}_{}{}{}_{}h{}m{}s.csv'.format(\
#     ENZ_pdb, SUB_pdb, str(resolution).zfill(3), YEAR, str(MONTH).zfill(2), str(DATE).zfill(2), HOUR, MINUTE, SECONDS)
# sorted_val_df.to_csv(output_filename, header = True, index=False, sep=',')    
# =============================================================================

'''
ENZ_CA = num_clash_atoms_ENZ; # of clash atoms in enzyme
SUB_CA = num_clash_atoms_SUB; # of clash atoms in substrate 
TOT_CA = total_num_clash_atoms; total # of clashing atoms    
clashes = total_num_clashes; total # of clashes
'''            
endTime = time.time()
task_time = endTime-startTime

num_rotations = len(alpha_range)*len(beta_range)*len(gamma_range)
#out_put_data = 'Output_ProteinEuler_{}_{}_{}{}{}_{}h{}m{}s.csv'.format(ENZ_pdb, SUB_pdb, \
       #                         YEAR, MONTH, DATE, HOUR, MINUTE, SECONDS)



# =============================================================================
# output = open(new_job_path+'/results_PE_{}_{}_res{}_{}{}{}_{}h{}m{}s.data'.format(\
#     ENZ_pdb, SUB_pdb, str(resolution).zfill(3), YEAR, str(MONTH).zfill(2), \
#         str(DATE).zfill(2), HOUR, MINUTE, SECONDS), 'w')
# =============================================================================
exe_time_end = datetime.datetime.now()
results = new_job_path+'/'+new_job_dir+'.data'
with open(results, 'w') as output:
    output.write('# Job title: MT67_6nj9_NCP job_1' + '\n')
    output.write('# Input ENZYME file:'+ os.getcwd()+'/'+ENZ_file.name + '\n')
    output.write('# Input SUBSTRATE file:'+os.getcwd()+'/'+SUB_file.name + '\n')
    output.write('# Number of ATOMs in ENZYME file: '+ str(ENZ_xyz.shape[0]) + '\n')
    output.write('# Number of ATOMs in SUBSTRATE file: '+str(SUB_xyz.shape[0]) + '\n')
    output.write('# Rotation scan range:' + '\n')
    output.write('# \t'+'alpha_range:'+ ' '+str(AlphaLower)+','+str(AlphaUpper) + '\n')
    output.write('# \t'+'beta_range:'+ ' '+str(BetaLower)+','+str(BetaUpper) + '\n')
    output.write('# \t'+'gamma_range:'+ ' '+str(GammaLower)+','+str(GammaUpper) + '\n')
    notation = '# Elemental rotations (extrinsic) in a sequence- zyx by gamma, beta and \
        alpha (Tait-Bryan angles), respectively. \
        \n# Extrinsic rotation- elementary rotations are carried out about \
            the fixed world coordinate system. '
    output.write(notation + '\n')
    output.write('# Clash criteria: >=0.4 vdw radii overlap' + '\n' )
    output.write('# Scan resolution: '+str(Resolution) + '\n')
    output.write('# Number of Rotations: '+str(num_rotations) + '\n')
    output.write('# Output file: '+str(output.name) + '\n')
    output.write('# Scan start Time: '+str(exe_time_start) + '\n')
    output.write('# Scan end time: '+str(exe_time_end) + '\n')
    output.write('# Execuation Time:' + str(datetime.timedelta(seconds = task_time)) + '\n') 
    output.write('# \n'*3)   
    output.write('# ENZ_CA: number of clashing atoms in ENZYME' + '\n')
    output.write('# SUB_CA: number of clashing atoms in SUBSTRATE' + '\n')
    output.write('# TOT_CA: total number of clashing atoms (ENZ_CA+SUB_CA)' + '\n')
    output.write('# clashes: total number of clashes' + '\n')
    output.write('# \n'*3)   
    output.write('# alpha, beta, gamma, ENZ_CA, SUB_CA, TOT_CA, clashes'  + '\n')  

    for cd in clash_data:
        #output.write(str(cd).strip('\[|\]')  + '\n')
        output.write('{:>5}\t{:>5}\t{:>5}\t{:<7}\t{:<7}\t{:<7}\t{:<7}'.format(\
                    cd[0],cd[1], cd[2], cd[3],cd[4], cd[5], cd[6])+'\n')

output.close()

