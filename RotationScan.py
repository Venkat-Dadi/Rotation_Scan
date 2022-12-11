#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 21:56:25 2022

@author: venkatareddy
"""

import numpy as np
from scipy.spatial.transform import Rotation as R
from scipy.spatial.distance import cdist
import datetime
import time as time
import pandas as pd
import sys
import os
from tqdm import tqdm
import configparser

YEAR = datetime.date.today().year
MONTH = datetime.date.today().month
DATE = datetime.date.today().day
HOUR = datetime.datetime.now().hour
MINUTE = datetime.datetime.now().minute
SECONDS = datetime.datetime.now().second

###############################################################################

def Get_Reacting_Atom_Coord(inputPDB, Reacting_Atoms):
    RAs = Reacting_Atoms.split(' ')
    if len(RAs) == 4:
        RA_chain = RAs[0]
        RA_resName = RAs[1]
        RA_resNum = RAs[2]
        RA_atomName = RAs[3]
    else:
        print('''Please provide Reaction_Atom details for '{}'':
        # chain_id residue_number residue_name atom_name'''.format(inputPDB))
        sys.exit()
    with open(inputPDB, 'r') as input_file:
        chain_lines = []
        for line in input_file.readlines():
            if line.startswith(('ATOM', 'HETATM')):
                if line[21:22].strip().casefold() == RA_chain.casefold():
                    chain_lines.append(line.strip('\n'))
        if len(chain_lines) == 0:
            sys.exit('Chain id \'{}\' is not found in \
                     PDB \'{}\''.format(RA_chain, inputPDB))
        else:
            pass
        
        resName_lines = []      
        for cl in chain_lines:
            if cl[17:20].strip().casefold() == RA_resName.casefold():
                resName_lines.append(cl.strip('\n'))
        if len(resName_lines) == 0:
            sys.exit('Residue name \'{}\' not found in chain \'{}\' of \
                     PDB \'{}\''.format(RA_resName, RA_chain, inputPDB))
        else:
            pass
        
        resNum_lines = []
        for rnaml in resName_lines:
            if rnaml[22:26].strip() == RA_resNum:
                resNum_lines.append(rnaml.strip('\n'))
        if len(resNum_lines) == 0:
            sys.exit('Residue number \'{}\' corresponding to\
 \'\{}:{}\' not found in PDB \'{}\''.format(RA_resNum, RA_chain, RA_resName, \
                                        inputPDB))
        else:
            pass
        
        atomName_line = []
        for rnuml in resNum_lines:
            if rnuml[12:16].strip().casefold() == RA_atomName.casefold():
                atomName_line.append(rnuml.strip('\n'))
        if len(atomName_line) == 0:
            sys.exit('Atom name \'{}\' not found in \
 \'\{}:{}{}\' of PDB \'{}\''.format(RA_atomName, RA_chain, RA_resName, \
                                    RA_resNum, inputPDB))
        if len(atomName_line) == 1:
            temp_arr0 = np.empty((0,3), float)          
            temp_xyz = np.fromstring(atomName_line[0][30:54], sep = ' ')
            temp_arr1 = np.append(temp_arr0, temp_xyz)
            RA_xyz = temp_arr1
            if len(RA_xyz) != 3:
                print(RA_xyz) 
                sys.exit('Something wrong with \'{}\' Reacting atom \
\'\{}:{}{}@{}\' coordinates'.format(inputPDB, RA_chain, RA_resName, RA_resNum,\
                                 RA_atomName))
            else:
                return RA_xyz
        else:
            alines = ''
            for al1 in atomName_line:
                alines += al1+'\n '
            sys.exit('More than one entry for \'\{}:{}{}@{}\' exists \
in PDB \'{}\'!!! Check your PDB!!!\n {}'.format(RA_chain, RA_resName, \
                                RA_resNum, RA_atomName, inputPDB, alines))

###############################################################################

def compute_rxncenCoord(A1_coord, A2_coord, vdw_sum_RAs):


    
    '''
    Computes the reaction center coordinates.
    Parameters
    ----------
    A1_coord : numpy array of shape (1,3)
        atomic coordinates of atom1 of enzyme's reacting atom
        (e.g., CE atom of SAM).
    A2_coord : numpy array of shape (1,3)
        atomic coordinates of atom2 that is connected to atom1 of 
        enzyme's reacting atom (e.g., SD atom of SAM).
    vdw_dist12 : float, optional
        sum of vdw radii of atom1 and substrate reacting atom
        e.g., CE atom of SAM and NZ atom of LYS.
        vdw radii are taken form Probe software(DOI: 10.1006/jmbi.1998.2400).
        The default is 3.3 (sum of vdw radii of CE and NA atoms)

    Returns
    -------
    numpy array of shape (1,3)
        location of reaction center in 3D coordinate system.

    '''

    rxn_vec = A1_coord - A2_coord 
    ''' vector that passes through bond between 
    atom1 and atom2; nothing but a reaction vector'''
    rxn_vec_norm = np.linalg.norm(rxn_vec)  # length of the above vector
    '''vdW radii C = 1.75A, N= 1.55A (Probe/Reduce)
    distance between CE and NZ is chosen to the the
    vdW distance between CE and NZ i.e 1.75+1.55'''
    unit_rxn_vec = rxn_vec/rxn_vec_norm
    rxn_cen_coord = A1_coord + vdw_sum_RAs*unit_rxn_vec
    return rxn_cen_coord        

###############################################################################

def Split_pdb_to_atomLines_xyz_vdw(inputPDB):
    

    '''
    Parameters
    ----------
    inputPDB : text file
    PDB file in .pdb format.

    Returns
    -------
    pdb_atomLines : TUPLE
        Contains PDB atom lines from column 1-30 (0 index of TUPLE) and
        column 55-80 (1 index of TUPLE)..
    xyz_vdw_array : numpy array
        pdb coordinate values in nx3 numpy array; n = # of atoms in the pdb.
    '''
    '''
    Extracts atomic coordinates from PDB file to an array and appends 
    vdw radius value to the each coordinate based on the ELEMENT record
    found in the atom line.
    
    vdw radii for atoms taken from Probe program
     ref: Probe program; DOI: 10.1006/jmbi.1998.2400 
    atom_vdw = {'H': -2.0, 'C': 1.7, 'N': 1.625, 'O': 1.480, 'P': 1.871, \
                'S': 1.782} # parm99 
    Hydrogen atoms are not considered in the computation of clashes.
    Due to that reason H vdw radius is assigned to be '-2'  
    '''
 
    atom_vdw = {'H': -2.0, 'C': 1.75, 'N': 1.55, 'O': 1.40, 'P': 1.80, 'S': 1.80}

    pdb_atomLines_0to31 = []
    pdb_atomLines_55to80 = []
    pdb_xyzvdw = []
    PDBfile = open(inputPDB, 'r')
    for pdbline in PDBfile.readlines():
        if pdbline.startswith(('ATOM', 'HETATM')):
            for key in atom_vdw:
                if pdbline[77:78] == key:
                    xyz_vdw_e = np.fromstring(pdbline[31:54]+' '+\
                    pdbline[61:66].replace(pdbline[61:66], str(atom_vdw[key])), sep = ' ')
                    pdb_xyzvdw.append(xyz_vdw_e)
            pdb_atomLines_0to31.append(pdbline[:30])
            pdb_atomLines_55to80.append(pdbline[54:81])
    xyz_vdw_array = np.array(pdb_xyzvdw)   
    pdb_atomLines = zip(pdb_atomLines_0to31, pdb_atomLines_55to80)
    PDBfile.close()              
    return pdb_atomLines, xyz_vdw_array

###############################################################################

def Move_SUB_rxnAtom2orgin(SUBpdb_coord_array, SUB_rxnAtom_array):
    '''
    Moves substrate coordinates in such way that the substrate's 
    reacting atom (NZ in case of protein LYS substrate) is placed at
    the origin.

    Parameters
    ----------
    SUBpdb_coord_array : numpy array of shape (N, 3)
        atomic coordinates for a substrate.
    SUB_rxnAtom_array : numpy array of shape (1,3)
        atomic coordinate values for a substrate's reacting atom.

    Returns
    -------
    Sub_atOrign_coord_array : numpy array of shape (N, 3)
        Transformed atomic coordinate array.

    '''
    Sub_atOrign_coord_array = SUBpdb_coord_array - SUB_rxnAtom_array
    return Sub_atOrign_coord_array

###############################################################################  

def Move_ENZ_rxnCen2origin(ENZpdb_coord_array, rxn_cen_coord):
    ENZpdb_coord_atOri = ENZpdb_coord_array - rxn_cen_coord 
    return ENZpdb_coord_atOri  

############################################################################### 

def Combine_atomLines_xyz_to_pdb(atomLines, xyz_coord):
    '''
    Combines PDB file atoms lines (column 1-30, 55-80) with 
    a input coordinate values in a numpy array format
    N x 3; N = # of atoms in a input pdb file.
    Number of atoms and the number of rows in the coordinate numpy array
    must be same.
    Parameters
    ----------
    atomLines_tuple : TUPLE
        Contains PDB atom lines from column 1-30 (0 index of TUPLE) and
        column 55-80 (1 index of TUPLE).
    xyz_coord : numpy array
        pdb coordinate values in nx3 numpy array; n = # of atoms in the pdb.

    Returns
    -------
    PDB file.

    '''
    out_atomLines_0to31 = []
    out_atomLines_55to80 = []
    for b1 in atomLines:
        out_atomLines_0to31.append(b1[0])
        out_atomLines_55to80.append(b1[1].strip('\n')) 
    
    if len(out_atomLines_0to31) != len(xyz_coord):
        sys.exit('Number of atoms do not match number of atom coordinates')
    
    else:
        coord = []    
        for c in range(len(xyz_coord)):
            str_x = np.array2string(xyz_coord[c,0])
            str_y = np.array2string(xyz_coord[c,1]) 
            str_z = np.array2string(xyz_coord[c,2]) 
            coord.append('{:8.3f}{:8.3f}{:8.3f}'.format(float(str_x), float(str_y), float(str_z)))
        
        combined_pdbLines = []
        for c1 in zip(out_atomLines_0to31, coord, out_atomLines_55to80):
            merged_atomLine = c1[0] + c1[1] + c1[2]+'\n'
            combined_pdbLines.append(merged_atomLine)
        return combined_pdbLines

###############################################################################

def writePDB(pdbLines_list, path, outputPrefix,  alpha = 0, beta = 0, gamma = 0):
    with open(path+'/'+outputPrefix+'_{}_{}_{}.pdb'.format(alpha, beta, gamma), 'w') as outPDB:
        for line in pdbLines_list:
            outPDB.write(line)   
            
###############################################################################

def rot_mat_align(vec1, vec2):
    '''
    adapted from https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    
    Generates a rotation matrix that aligns two vectors onto one another.
    Parameters
    ----------
    vec1 : numpy array of shape (3,3)
         vector.
    vec2 : numpy array of shape (1,3)
        vector.

    Returns
    -------
    rotation_matrix : numpy array of shape (3,3)
        returns a rotation matrix that align vec1 and vec2.

    '''
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
        ])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix    

###############################################################################

def Pairwise_vdw_SumLimit(ENZ_vdw, SUB_vdw, VDW_Clash_Threshold = 0.4):
    '''
    Creates an 2D array of shape (N, M) from two 1D arrays (A, B) with N and M atoms, respectively.
    Each element of the 2D array is a sum of vdw radii of i atom of A and j atom of B.
    N: number of atoms in Enzyme pdb
    M: number of atoms in Substrate pdb
    
    Parameters
    ----------
    ENZ_vdw : 1d array of size N
        Elements of the array are the vdw radius of atoms from Enzyme PDB.
    SUB_vdw : 1d array of size M
        Elements of the array are the vdw radius of atoms from Substrate PDB.
    VDW_Clash_Threshold : float, optional
        Value for clash tolerance in Angstronms. This value is subtracted from
        the sum of vdw radii of atom pairs (Ai, Bj). 
        The default is 0.4 (reference: Probe, DOI: 10.1006/jmbi.1998.2400)

    Returns
    -------
    vdw_sum_limit : 2D array
        sum of vdw radii of Ai and Bj minus VDW_Clash_Threshold.
    '''
    ENZ_vdw_m, SUB_vdw_m = np.meshgrid(ENZ_vdw, SUB_vdw) 
    vdw_sum = np.add(ENZ_vdw_m, SUB_vdw_m) # sum of two matrices 
    vdw_sum_limit = vdw_sum - VDW_Clash_Threshold   # clash threshold matrix
    #vdw_sum_limit_avr = vdw_sum*0.5
    return vdw_sum_limit

###############################################################################

def RotationScan(MOBILE_xyz, STATIC_xyz, alpha, beta, gamma, Pairwise_vdw_SumLimit):
    '''
    Rotates the MOBILE structure relative to STATIC structure using
    set of elemental rotations about X, Y and Z-axes. For each rotation, 
    it calculates the number of clashing atom/clashes between MOBILE and STATIC 
    structures. 
    Rotation is extrinsic i.e. axes attaced to the object do not move.
    Global frame of reference
    Parameters
    ----------
    MOBILE_xyz : numpy array of shape (N, 3)
        Coordinates of the MOBILE structure.
    STATIC_xyz : numpy array of shape (N, 3)
        Coordinates of the STATIC structure.
    alpha : float or int
        elemental rotation angle about X-axis.
    beta : float or int
        elemental rotation angle about Y-axis.
    gamma : float or int
        elemental rotation angle about Z-axis.
    Pairwise_vdw_SumLimit : numpy array of shape (N, M)
    Returns
    -------
    rotated_MOBILE_xyz : numpy array of shape (N, 3)
        Transformed MOBILE coordinates after rotation.
    clash_data_one : list
        Number of Clashes/clasing atoms as a function of alpha, beta and gamma.
        [alpha, beta, gamma, num_clash_atoms_MOBILE, num_clash_atoms_STATIC, total_num_clash_atoms, total_num_clashes]
        num_clash_atoms_MOBILE: Number of clasing atoms from MOBILE structure
        num_clash_atoms_STATIC: Number of clasing atoms from STATIC structure
    '''
    r = R.from_euler('zyx', [gamma, beta, alpha], degrees=True ) # generates euler rotation zyx extrinsic 
    Rot_matrix = r.as_matrix() # Euler rotation as matrix
    rotated_MOBILE_xyz = np.matmul(Rot_matrix, MOBILE_xyz.T).T.round(3) # multiplication of xyz coord with rotation matrix
    pair_wise_dist = cdist(STATIC_xyz, rotated_MOBILE_xyz, 'euclidean') # calculates pairwise distances between rotated ligand coord and receptor coord
    ##################################################################
    clash_matrix = np.less_equal(pair_wise_dist, Pairwise_vdw_SumLimit).astype(dtype=int) # element-wise comparison of two ndarrays
    ##################################################################
    num_clashes_mobile = np.sum(clash_matrix, axis = 0) # num. of clashed with enzyme involved array
    num_clashes_static = np.sum(clash_matrix, axis = 1) # num. of clashed with substrate involved array
    total_num_clashes = np.sum(clash_matrix)  # total number of clashes 
    num_clash_atoms_mobile_array = np.greater_equal(num_clashes_mobile, 1).astype(dtype=int)
    num_clash_atoms_mobile = np.sum(num_clash_atoms_mobile_array)
    num_clash_atoms_static_array = np.greater_equal(num_clashes_static,1).astype(dtype=int)
    num_clash_atoms_static = np.sum(num_clash_atoms_static_array)
    total_num_clash_atoms = num_clash_atoms_mobile + num_clash_atoms_static
    clash_data_one = [alpha, beta, gamma, num_clash_atoms_mobile, num_clash_atoms_static, total_num_clash_atoms, total_num_clashes]
    return rotated_MOBILE_xyz, clash_data_one   

###############################################################################

def Run(parametersFile):
    '''
    Runs the Protein Euler Script on input PDB files.

    Parameters
    ----------
    parametersFile : Configureation file
        Contains input data files and parameters to run the scrip.
        e.g., parameter.ini .

    Returns
    -------
    Run log file with clash data (*.data) and only clash data (*.csv).
    Writes PDB files of MOBILE pdb at every rotation.

    '''

    Execution_Start_Time = datetime.datetime.now()
    startTime = time.time()
    
    config = configparser.ConfigParser()
    config.read(parametersFile)
    
    EnzymePDB = config['INPUT_FILES']['EnzymePDB'].split('#')[0].strip()
    SubstratePDB = config['INPUT_FILES']['SubstratePDB'].split('#')[0].strip()
    JOB_Name = config['PARAMETERS']['JOB_Name'].split('#')[0].strip()
    EnzymePrefix = config['PARAMETERS']['EnzymePrefix'].split('#')[0].strip()
    SubstratePrefix = config['PARAMETERS']['SubstratePrefix'].split('#')[0].strip()
    Enzyme_Reaction_Atom = config['PARAMETERS']['Enzyme_Reaction_Atom'].split('#')[0].strip()
    Substrate_Reaction_Atom = config['PARAMETERS']['Substrate_Reaction_Atom'].split('#')[0].strip()
    Enzyme_Scissile_Atom2 = config['PARAMETERS']['Enzyme_Scissile_Atom2'].split('#')[0].strip()
    Resolution = float(config['PARAMETERS']['Resolution'].split('#')[0].strip())
    AlphaLower = float(config['PARAMETERS']['AlphaLower'].split('#')[0].strip())
    AlphaUpper = float(config['PARAMETERS']['AlphaUpper'].split('#')[0].strip())
    BetaLower = float(config['PARAMETERS']['BetaLower'].split('#')[0].strip())
    BetaUpper = float(config['PARAMETERS']['BetaUpper'].split('#')[0].strip())
    GammaLower = float(config['PARAMETERS']['GammaLower'].split('#')[0].strip())
    GammaUpper = float(config['PARAMETERS']['GammaUpper'].split('#')[0].strip())
    VDW_Radii_Sum_RAs = float(config['PARAMETERS']['VDW_Radii_Sum_RAs'].split('#')[0].strip())
    VDW_Clash_Threshold = float(config['PARAMETERS']['VDW_Clash_Threshold'].split('#')[0].strip())
    WriteOut_transfromed_PDBs = config['PARAMETERS']['WriteOut_transformed_PDBs'].split('#')[0].strip().upper()
    
    mandate_parameters = [EnzymePDB, SubstratePDB, EnzymePrefix, SubstratePrefix, 
                  Enzyme_Reaction_Atom, Enzyme_Scissile_Atom2, Substrate_Reaction_Atom,
                  Resolution, AlphaLower, AlphaUpper, BetaLower, BetaUpper, GammaLower,
                  GammaUpper, VDW_Radii_Sum_RAs, VDW_Clash_Threshold, WriteOut_transfromed_PDBs]
    if '' in mandate_parameters:
        sys.exit('One or more parameters are missing form \'{}\' file'.format(parametersFile))
    else: 
        print('####################################################')
        print('####   Rotation Scan using ProteinEuler script  ####')
        print('####################################################')
        print('''INPUT AND RUN PARAMETERS:
    EnzymePDB:                       {}
    SubstratePDB:                    {}
    JOB_Name:                        {}
    EnzymePrefix:                    {}
    SubstratePrefix:                 {}
    Enzyme_Reaction_Atom:            {}
    Enzyme_Scissile_Atom2:           {}
    Substrate_Reaction_Atom:         {}
    Resolution:                      {}\u00b0
    AlphaLower:                      {}\u00b0
    AlphaUpper:                      {}\u00b0
    BetaLower:                       {}\u00b0
    BetaUpper:                       {}\u00b0
    GammaLower:                      {}\u00b0
    GammaUpper:                      {}\u00b0
    VDW_Radii_Sum_RAs:               {} \u212B
    VDW_Clash_Threshold:             {} \u212B
    WriteOut_transformed_PDBs:       {}
    '''.format(EnzymePDB.split('/')[-1], SubstratePDB.split('/')[-1], JOB_Name, EnzymePrefix, SubstratePrefix,\
    Enzyme_Reaction_Atom, Enzyme_Scissile_Atom2, Substrate_Reaction_Atom, Resolution,\
    AlphaLower, AlphaUpper, BetaLower, BetaUpper, GammaLower, GammaUpper, \
    VDW_Radii_Sum_RAs, VDW_Clash_Threshold, WriteOut_transfromed_PDBs))
    
    ENZ_Reacting_Atom_xyz = Get_Reacting_Atom_Coord(EnzymePDB, Enzyme_Reaction_Atom)
    ENZ_Scissile_Atom_xyz = Get_Reacting_Atom_Coord(EnzymePDB, Enzyme_Scissile_Atom2)
    SUB_Reacting_Atom_xyz = Get_Reacting_Atom_Coord(SubstratePDB, Substrate_Reaction_Atom)

    Rxn_Cen_xyz = compute_rxncenCoord(ENZ_Reacting_Atom_xyz, ENZ_Scissile_Atom_xyz, VDW_Radii_Sum_RAs)
    
    ENZ_AtomLines, ENZ_xyz_vdw = Split_pdb_to_atomLines_xyz_vdw(EnzymePDB)
    SUB_AtomLines, SUB_xyz_vdw = Split_pdb_to_atomLines_xyz_vdw(SubstratePDB)
    
    ENZ_AtomLines = list(ENZ_AtomLines)
    SUB_AtomLines = list(SUB_AtomLines)
    
    ENZ_xyz_atOrigin = Move_ENZ_rxnCen2origin(ENZ_xyz_vdw[:, :3], Rxn_Cen_xyz)
    SUB_xyz_atRxnCen = Move_SUB_rxnAtom2orgin(SUB_xyz_vdw[:, :3], SUB_Reacting_Atom_xyz)
    
   
    ENZ_atOrigin_centroid = (ENZ_xyz_atOrigin.sum(axis = 0)/len(ENZ_xyz_atOrigin))
    SUB_atOrigin_centroid = (SUB_xyz_atRxnCen.sum(axis = 0)/len(SUB_xyz_atRxnCen))
    
    positive_Y = np.array([0.000, abs(np.linalg.norm(ENZ_atOrigin_centroid)), 0.000])
    negative_Y = np.array([0.000, -np.linalg.norm(SUB_atOrigin_centroid), 0.000])
    
    rotMat_for_Postive_Y = rot_mat_align(ENZ_atOrigin_centroid, positive_Y)
    rotMat_for_Negative_Y = rot_mat_align(SUB_atOrigin_centroid, negative_Y)
    
    ENZ_onPositive_Y = np.matmul(rotMat_for_Postive_Y, ENZ_xyz_atOrigin.T).T
    SUB_onNegative_Y = np.matmul(rotMat_for_Negative_Y, SUB_xyz_atRxnCen.T).T
    
    ENZ_xyz = ENZ_onPositive_Y
    ENZ_vdw = ENZ_xyz_vdw[:, -1]
    SUB_xyz = SUB_onNegative_Y
    SUB_vdw = SUB_xyz_vdw[:, -1]
    
    pairwise_vdw_sum_limit = Pairwise_vdw_SumLimit(ENZ_vdw, SUB_vdw, VDW_Clash_Threshold)
        
    if len(ENZ_xyz) <= len(SUB_xyz):
        MOBILE_xyz = ENZ_xyz
        MOBILE_AtomLines = ENZ_AtomLines
        MOBILE_PDB_Prefix = EnzymePrefix
        STATIC_PDB_Prefix = SubstratePrefix
        STATIC_xyz = SUB_xyz
        print('''Enzyme PDB has less number of atoms compared to Substrate PDB.
    Hence Enzyme PDB is chosen as MOBILE and Substrate PDB as STATIC for Rotation Scan''')
    else:
        MOBILE_xyz = SUB_xyz
        MOBILE_AtomLines = SUB_AtomLines
        MOBILE_PDB_Prefix = SubstratePrefix
        STATIC_PDB_Prefix = EnzymePrefix
        STATIC_xyz = ENZ_xyz
        print('''Substrate PDB has less number of atoms compared to Enzyme PDB.
Hence Substrate PDB is chosen as MOBILE and Enzyme PDB as STATIC for Rotation Scan.''')
    

    Resolution = int(Resolution)
    alpha_range = np.arange(AlphaLower, AlphaUpper+Resolution, Resolution)  # rotation about z axis in zyx extrinsic convention
    beta_range = np.arange(BetaLower, BetaUpper+Resolution, Resolution)  # rotation about y axis in zyx extrinsic convention
    gamma_range = np.arange(GammaLower, GammaUpper+Resolution, Resolution)  # rotation about x axis in zyx extrinsic convention
    
    pwd_path = os.getcwd()
    Results_dir = 'Results_PE_{}_{}_{}_res{}_{}{}{}_{}h{}m{}s'.format(JOB_Name, EnzymePrefix, SubstratePrefix, str(Resolution).zfill(3), YEAR, str(MONTH).zfill(2), str(DATE).zfill(2), HOUR, MINUTE, SECONDS)
    Results_dir_path = os.path.join(pwd_path, Results_dir)
    transformedPDB_output_dir = 'Rotated_PDBs_{}_{}{}{}_{}h{}m{}s'.format(MOBILE_PDB_Prefix, YEAR, str(MONTH).zfill(2), str(DATE).zfill(2), HOUR, MINUTE, SECONDS)
    transformedPDB_output_dir_path = os.path.join(Results_dir_path, transformedPDB_output_dir)
    os.mkdir(Results_dir_path)
    if WriteOut_transfromed_PDBs == 'TRUE':
        os.mkdir(transformedPDB_output_dir_path)
    
    ENZ_AtomLines_xyz_onY = Combine_atomLines_xyz_to_pdb(ENZ_AtomLines, ENZ_onPositive_Y.round(3))
    SUB_AtomLines_xyz_onY = Combine_atomLines_xyz_to_pdb(SUB_AtomLines, SUB_onNegative_Y.round(3))
    
    writePDB(ENZ_AtomLines_xyz_onY, Results_dir_path, 'input_' + EnzymePrefix)
    writePDB(SUB_AtomLines_xyz_onY, Results_dir_path, 'input_' + SubstratePrefix)
    results_file_data = Results_dir_path+'/'+Results_dir+'.data'
    results_file_csv = Results_dir_path+'/'+Results_dir+'.csv'
    
    Num_Rotations = len(alpha_range)*len(beta_range)*len(gamma_range)   
    print('Total Number of Rotations: {}'.format(Num_Rotations))
    print('Running Rotation Scan ...')
    clash_data = []
    pbar = tqdm(total = len(alpha_range)*len(beta_range)*len(gamma_range))
    for alpha in alpha_range:
        for beta in beta_range:
            for gamma in gamma_range:
                rotated_mobile_xyz, clash_data_one = RotationScan(MOBILE_xyz, STATIC_xyz, alpha, beta, gamma, pairwise_vdw_sum_limit) 
                clash_data.append(clash_data_one)
                if WriteOut_transfromed_PDBs == 'TRUE':
                    combined_rotated_MOBILE_xyz_AtomLines = Combine_atomLines_xyz_to_pdb(MOBILE_AtomLines, rotated_mobile_xyz)
                    writePDB(combined_rotated_MOBILE_xyz_AtomLines, transformedPDB_output_dir_path, MOBILE_PDB_Prefix, int(alpha), int(beta), int(gamma))
                else:
                    pass
                pbar.update(1)
    pbar.close() 
    
    clash_data_array = np.array(clash_data)
    columns = ['alpha', 'beta', 'gamma', 'MOB_CA', 'STAT_CA', 'TOT_CA', 'TOT_CLASHES']
    df = pd.DataFrame(clash_data_array, columns = columns)
    df.to_csv(results_file_csv, header = True, index = False)
    
    Execution_End_Time = datetime.datetime.now()
    endTime = time.time()
    taskTime  = endTime - startTime
    
    print('''\nWriting data to \'{}.data\'
       \'{}.csv\''''.format(Results_dir, Results_dir))
    with open(results_file_data, 'w') as output:
        output.write('''####################################################
####   Rotation Scan using ProteinEuler script  ####
####################################################
# INPUT AND RUN PARAMETERS:
#        EnzymePDB:                       {}
#        SubstratePDB:                    {}
#        JOB_Name:                        {}
#        EnzymePrefix:                    {}
#        SubstratePrefix:                 {}
#        Enzyme_Reaction_Atom:            {}
#        Enzyme_Scissile_Atom2:           {}
#        Substrate_Reaction_Atom:         {}
#        Resolution:                      {}\u00b0
#        AlphaLower:                      {}\u00b0
#        AlphaUpper:                      {}\u00b0
#        BetaLower:                       {}\u00b0
#        BetaUpper:                       {}\u00b0
#        GammaLower:                      {}\u00b0
#        GammaUpper:                      {}\u00b0
#        VDW_Radii_Sum_RAs:               {} \u212B
#        VDW_Clash_Threshold:             {} \u212B
# PDB DETAILS:        
#        Number of ATOMS (EnzymePDB)      {}
#        Number of ATOMS (SubstratePDB)   {}
#        
# ROTATION DETAILS:
#        ROTATION SCAN RANGE:
#                     alpha:              ({}\u00b0, {}\u00b0)
#                      beta:              ({}\u00b0, {}\u00b0)
#                     gamma:              ({}\u00b0, {}\u00b0)
#        Scan Resolution:                 {}\u00b0
#        ROTATION TYPE:                   Extrinsic elemental rotiaons in a sequence about z, y, and x
#                                         axes by gamma, beta and alpha, respectively (Tait_Bryan convention).
#                                         Extrinsic rotation- elementary rotations are carried out about the fixed world coordinate system. 
#        Number of Composite Rotations:   {};  Each composite rotation contains 3 elemental rotations.
#        MOBILE PDB:                      {}
#        STATIC PDB:                      {}   
#
# OUTPUT DETAILS:
#        Output Directory:                {}
#            Total output file (.data):   {}.data
#            Clash data (.csv):           {}.csv
#        Rotation Scan Start Time:        {}
#        Rotation Scan End Time:          {}
#        Execution Time:                  {} sec
#        WriteOut_transformed_PDBs:       {}
#
#
#        MOB_CA:                          Number of Clashing atoms from the \'{}\' PDB
#        STAT_CA:                         Number of Clashing atoms from the \'{}\' PDB
#        TOT_CA:                          Total number of clashing atoms ( = MOB_CA + STAT_CA)
#        TOT_CLASHES:                     Total number of clashes
#
#
# alpha, beta, gamma, MOB_CA, STAT_CA, TOT_CA, TOT_CLASHES
'''.format(EnzymePDB.split('/')[-1], SubstratePDB.split('/')[-1], JOB_Name, EnzymePrefix, SubstratePrefix,\
        Enzyme_Reaction_Atom, Enzyme_Scissile_Atom2, Substrate_Reaction_Atom, Resolution,\
        AlphaLower, AlphaUpper, BetaLower, BetaUpper, GammaLower, GammaUpper, \
        VDW_Radii_Sum_RAs, VDW_Clash_Threshold, len(ENZ_xyz), len(SUB_xyz),\
        AlphaLower, AlphaUpper, BetaLower, BetaUpper, GammaLower, GammaUpper,\
        Resolution, Num_Rotations, MOBILE_PDB_Prefix, STATIC_PDB_Prefix, Results_dir_path,\
            Results_dir, Results_dir,\
        Execution_Start_Time, Execution_End_Time, round(taskTime, 3), WriteOut_transfromed_PDBs,\
        MOBILE_PDB_Prefix, STATIC_PDB_Prefix))
            
        for cd in clash_data:
            output.write('{:>6}\t{:>6}\t{:>6}\t{:<7}\t{:<7}\t{:<7}\t{:<7}'.format(\
                        cd[0],cd[1], cd[2], cd[3],cd[4], cd[5], cd[6])+'\n')

###############################################################################

if len(sys.argv) != 2:
        print('####################################################')
        print('####   Rotation Scan using ProteinEuler script  ####')
        print('####################################################')
        print('''\nWARNING: Input parameters file is missing.
              
Usage: python ProteinEuler.py parameters.ini

Example parameters.ini:   
    
[INPUT_FILES]
EnzymePDB =         <pdb file>             # input enzyme pdb file (.pdb)reaction center at
				                           # origin and molecule is aligned to desired 
                                           # position (eg. molecule's centroid on y-axis ), 
                                           # if required.
SubstratePDB =     <pdb file>              # input substrate pdb file (.pdb) reaction 
                                           # center at origin and molecule is aligned 
                                           # to desired position (eg. molecule's 
                                           # centroid on y-axis ), if required.
                                     
[PARAMETERS]
JOB_Name =      <string>                   # Name of the job
EnzymePrefix = <string>                    # name of an enzyme used as suffix for output (eg. dot1L).
SubstratePrefix = <string>                 # name of the substrate used as suffix for file output (eg. NCP).
Substrate_Reaction_Atom = <string>         # chain_id residue_number residue_name atom_name
Enzyme_Reaction_Atom = <string>            # chain_id residue_number residue_name atom_name   
Enzyme_Scissile_Atom2 = <string>           # chain_id residue_number residue_name atom_name; 
                                           # it is the atom bonded to enzyme reacting atom             
Resolution = int                           # Resolution of a rotation scan. High resolution scans are computationally expensive                       

# Trait-Bryan angles for roation in 3D                          
AlphaLower = int                           # lower value of an rotation angle range about x-axis
AlphaUpper = int                           # upper value of an rotation angle range about x-axis  
BetaLower = int                            # lower value of an rotation angle range about y-axis
BetaUpper = int                            # upper value of an rotation angle range about y-axis
GammaLower =  int                          # lower value of an rotation angle range about z-axis
GammaUpper = int                           # upper value of an rotation angle range about z-axis 
VDW_Radii_Sum_RAs = float                  # is the sum of van der Waal's radii of two reacting atoms. 
                                           # e.g., in case of methyltransferase, it the sum of vdw radii of CE of SAM and NE of Lys of substrate.
VDW_Clash_Threshold = float                # Atom overlap in Angstroms for assigning the clash between atoms.
WriteOut_transformed_PDBs = boolean        # TRUE or FALSE; wether to write pdb files from the roation scan or not. 
''')   
        sys.exit()
        
else:
    Run(sys.argv[1])  
