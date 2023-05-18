# **Rotation Scan**

To check the accessibility of the *reaction center* of an enzyme to the its substrate by modeling the enzyme-substrate *reaction complex* for protein methylatransferases.

**Reaction complex:** a complex formed between an enzyme and its substrate during an enzymatic reaction. \
**Reaction center:** location of the methyl acceptor atom (NZ atom of lysine, in case of protein lysine methyltransferases) the substrate during the methyl group transfer. \

During protein lysine methylation, a *reaction complex* must be formed between the methyltransferase and its protein substrate (true for any enzymatic reaction).  The accessibility of the enzyme's *reaction center* to the methyl acceptor substrate depends on the milieus of the *reaction center* and the acceptor atom of the substrate. The accessibility can be calculated by modeling the enzyme-substrate *reaction complex*. In such complex, during the methyl transfer, the methyl acceptor atom of the substrate (NZ atom of lysine in case of lysine methyltransferases) lies at the *reaction center* in contact with methyl carbon atom (CE) of S-adenosyl-L-methionine. The substrate must approach the *reaction cetner* such that the binding pose between the enzyme and its substrate has no steric hidrance. The pose with minimum number of steric clashes can be calculated by constructing the enzyme-substrate *reaction complex* by imposing the criteria of the *reaction center* and rotating either the enzyme or substrate about the *reaction center* while the other is fixed immobile. Rotations were carried out with the help of Euler's rotation theorem, that covers all the possible orientation that a substrate can approach the *reaction center* of an enzyme. All these rotations constitutes the *Rotation Scan*. For each rotation, the number of atoms invovled in the clashes are calculated according to the clash criteria used by Probe ([DOI: 10.1006/jmbi.1998.2400](https://www.sciencedirect.com/science/article/pii/S0022283698924007?via%3Dihub)). The pose with minimum number of atomic clashes is considered as the possible binding mode of the substrate to the enzyme. This is implemented using python script: `RotationScan.py`.

Each pose of a *Rotation Scan* is represented by three elemental roations (α, β and γ) that constitute a single rotation. The total number of clashing atoms (TOT_CA or TCA) are plotted as a funciton of these elemental roations. For simplicity, two elemental rotation combinations were plotted while the third corresponds to the minimum TCA value.

**Caveats:**
* This model does not consider the flexibility either in the substrate or the enzyme.
* This model can only be used to show numberically that the inaccessibility of the active stie to the bulky substrates due to the steric hindrance.

## **System requirements:**
* Operating system: Linux, Windows.
* 4 - 8 GB RAM

`RotationScan.py` is a python script written in python3. 
The code is tested on Ubuntu 20.04 operating system with Intel i5 or i7 processor and 16 GB RAM. The usage of memory depends on the number of atoms in the model system. \
**Tested sytems and memory usage:**

Number of atoms <br> in a model system | Memory usage
--------------------------------|----------------------
14764 | 1.12 GB
18308 | 2.5 GB 


**Dependencies:**
* python3.x
* pip
* numpy
* pandas
* scipy
* matplotlib
* tqdm

**Installation of dependencies:** \
`python3.x` comes preinstalled on Ubuntu 20.04 or else can be installed before installing `pip` and other modules as follows: 

* Ubuntu and Debian: open a terminal on Ubuntu by pressing `Ctrl + Alt + T` and type the follwoing in the prompt.
```
sudo apt update
sudo apt install python3
sudo apt install python3-pip
pip install numpy
pip install pandas
pip install scipy
pip install tqdm
pip install matplotlib
```
Installation of the above dependencies takes a few minutes.

**Note:** Other Linux based operating systems and Windows with python3, pip and other python modules (numpy, pandas, scipy, tqdm and matplotlib) installed can use `RotationScan.py` as python script is OS-agnostic. However, the script is not tested on those operating systems.

## **How rotation scan works:**
The script is developed and explained in the context of **protein lysine methytransferases** and their **protein substrates**.
### **Step 1: Preparation of enzyme-substrate complex model system for rotation scan**
* The reacting atoms of enzyme and substrate are identified from the `pararmeters_RotationScan.ini` file. \
`Substrate_Reaction_Atom`: reacting atom of the substrate. Methyl acceptor atom (nitrogen, `NZ`) in case of lysine substrate. \
`Enzyme_Reaction_Atom`: reacting atom of the enzyme cofactor. Methyl carbon atom (`CE`) in case of cofactor/cosubstrate SAM.
`Enzyme_Scissile_Atom2`: atom covalently linked to the `Enzyme_Reaction_Atom`. Sulfur atom (`SD`) in case of cofactor/cosubstrate SAM.
* The coordinates of the reacting atoms are extracted.
* The reaction axis vector is computed from the coordinates of `Enzyme_Reaction_Atom` and `Enzyme_Scissile_Atom2`. The reaction axis is the line that passes through `Enzyme_Reaction_Atom` and `Enzyme_Scissile_Atom2`. 
* A point on the reaction axis vector is defined on the *reaction center*. The *reaction center* is the point where the substrate acceptor atom (`Substrate_Reaction_Atom`) lies during methylation and it lies about 3.3Å (`VDW_Radii_RAs`) away from the `Enzyme_Reaction_Atom`. This is the sum of van der Waals radii of `Substrate_Reaction_Atom` (`NZ` of lysine substrate) and `Enzyme_Reaction_Atom` (`CE` atom of SAM).
* The coordinates of the enzyme molecule is moved such that the *reaction center* is at the origin (0,0,0) of the three dimensional cartesian coordinate system.
* Similarly, the coordinates of the substrate molecules are moved to the origin. This makes the substrate reacing atom lies at the *reaction center* of the enzyme. This resembles a enzyme-substrate complex.
* To start with the initial orientation of enzyme and substrate molecules, enzyme is transformed such that its centroid is placed on the positive y-axis. Whereas, the substrate's centroid is place on the negative y-axis. This gives rise to the initial complex model on which the rotation scan is executed.

### **Step 2: Execution of rotation scan**
* The molecule with less number of atoms is identified and regarded as the mobile molecule which under goes rotation and the other is static.
* From the `pararmeters_RotationScan.ini` file, the range of elemental rotation angles (`AlphaLower`, `AlphaUpper`, `BetaLower`, `BetaUpper`, `GammaLower`, `GammaUpper`) are obtained. The `Resolution` determines the interval of elemental rotation angle range.
* The permutation of three angles (α, β and γ) for whole range of elemental rotation angles is prepared and for each set the rotation matrix is computed.
* The rotation matrix was used to transform the coordinates of the mobile molecule. For each rotation of mobile molecule, the number of clashes and the clasing atoms between the mobile and static molecules are computed. This is repeated for all rotations.
* The number of clashing atoms and clashes were output as a function of rotation (α, β, γ).



## **How to run rotation scan:**
Type the follwing in the terminal prompt.
```
python RotationScan.py parameters_RotationScan.ini

```

The input to the `RotationScan.py` is supplied in `parameters_RotationScan.ini` file. The content of the file as follows:

**parameters_RotationScan.ini file:**
```
[INPUT_FILES]
EnzymePDB =         <enzyme pdb>              # input enzyme pdb file (*.pdb)
				          
SubstratePDB =     <substrate pdb>               # input substrate pdb file (*.pdb) 
                                                                              
[PARAMETERS]
JOB_Name =      <job name>                           # Name of the job
EnzymePrefix = <enzyme prefix>                       # name of an enzyme used as suffix for output (eg. dot1L).
SubstratePrefix = <substrate prefix>                      # name of the substrate used as suffix for file output (eg. NCP).
Substrate_Reaction_Atom = <pdb chain id> <residue name> <residue number> <atom name> # eg. In "E LYS 79 NZ", chain id = E, residue name = LYS, residue number = 79 atom name = NZ
                                                      # This atom involve in enzymatic reaction.

Enzyme_Reaction_Atom = <pdb chain id> <cofactor name> <cofactor number> <atom name> # eg. In "K SAM 500 CE ", chain id = K, residue name = SAM, residue number = 500 atom name = CE       
						# This atom of cofactor involve in enzymatic reaction.
Enzyme_Scissile_Atom2 = <pdb chain id> <cofactor name> <cofactor number> <atom name> # eg. In "K SAM 500 SD", chain id = K, residue name = SAM, 
									             # residue number = 500 atom name = SD       
                                                                                     # it is the atom bonded to enzyme's cofactor's reacting atom             
Resolution =  <integer>                           # Resolution of a rotation scan. High resolution scans 
					  # are computationally expensive. As the resolution increases, the number of elemental rotations increase exponentially.                       

# Trait-Bryan angles for roation in 3D                       
# Each rotation in a rotation scan constitutes 3 elemental rotations alpha, beta and gamma about x, y, and z-axes, respectively.
# A scan range of -180 to + 180 for α, -90 to +90 for β and -180 to +180 for γ covers all the space 
# that an object under roation could explore in three-dimensions.   
AlphaLower = <integer>                          # lower value of an elemental rotation angle (alpha) about x-axis
AlphaUpper = <integer>                            # upper value of an elemental rotation angle (alpha) about x-axis 
BetaLower = <integer>                             # lower value of an elemental rotation angle (beta) about y-axis
BetaUpper = <integer>                              # upper value of an elemental rotation angle (beta) about y-axis
GammaLower = <integer>                           # lower value of an elemental rotation angle (gamma) about z-axis
GammaUpper = <integer>                            # upper value of an elemental rotation angle (gamma) about z-axis

VDW_Radii_Sum_RAs = <float>                    # is the sum of van der Waal's radii of two reacting atoms. 
                                           # e.g., in case of methyltransferase, it the sum of vdw radii of CE of SAM and NE of Lys of substrate.
                                           # the vdw radii for atoms are taken from  taken from Probe program
                                           # ref: Probe program; DOI: 10.1006/jmbi.1998.2400 
                                           # 'atom': <vdw> = {'H': -2.0, 'C': 1.7, 'N': 1.625, 'O': 1.480, 'P': 1.871, 'S': 1.782}
VDW_Clash_Threshold = 0.4                  # Atom overlap in Angstroms for assigning the clash between atoms. Default is 0.4 according to Probe program
WriteOut_transformed_PDBs = false           # TRUE or FALSE; whether to write pdb files from the roation scan or not.
```
## **Output of rotation scan run:**
`RotationScan.py` run generates a results directory with a name `Results_PE_<enzyme prefix>_<substrate prefix>_res<int>_<run date>_<run time>`. The results directory contains:
* **input enzyme pdb file** (`input_<enzyme prefix>_0_0_0.pdb`) before rotation scan start.
* **input substrate pdb file** (`input_<substrate prefix>_0_0_0.pdb`) before rotation scan start.
* `Results_PE_<enzyme prefix>_<substrate prefix>_res<int>_<run date>_<run time>.data` file with input parameters used in the rotation scan followd by three elemental rotation angles and their corresponding number of clashes and clasing atoms.
* `Results_PE_<enzyme prefix>_<substrate prefix>_res<int>_<run date>_<run time>.csv` file that contains elemental rotation angles and their corresponding number of clashes and clasing atoms: `alpha,beta,gamma,MOB_CA,STAT_CA,TOT_CA,TOT_CLASHES`. \
`alpha`: rotation about x-axis. \
`beta`: rotation about y-axis. \
`gamma`: rotation about z-axis. \
`MOB_CA`: number of clashing atoms in the mobile molecule (molecule with less number of atoms is rotated). \
`STAT_CA`: number of clashing atoms in the static molecule (molecule with more number of atoms is kept static). This reduces the cost of computation. \
`TOT_CA`: Total number of clasing atoms (TCA) i.e., `MOB_CA + STAT_CA`.
`TOT_CLASHES`: Total number of clashes between Mobile and static or enzyme and substrate.

## **Plotting rotation scan results:**
Type the following in the terminal prompt.
```
python plotting_RotationScan_results.py Results_PE_<enzyme prefix>_<substrate prefix>_res<int>_<run date>_<run time>.data
```
It generates two files:
* A `Results_PE_<enzyme prefix>_<substrate prefix>_res<int>_<run date>_<run time>.pdf` file heat map showing the total clasing atoms (TCA) as a function of elemental rotation angles (α, β and γ). Since the rotation scan generates 4-dimensional data (α, β, γ, TCA), three heat maps are prepared at each elemental rotation angle with minimum TCA value.
i.e., TCA(β, γ, α= at minimum TCA value), TCA(α, γ, β= at minimum TCA value) and TCA(β, α, γ= at minimum TCA value).
* A `Results_PE_<enzyme prefix>_<substrate prefix>_res<int>_<run date>_<run time>_sorted.csv` file with TOT_CA sorted in the ascending order.

## **Normal run time:**
A reaction complex model system with mobile molecule (6195 atoms) and stationary molecule (18209 atoms) for 26011 rotations takes around 8 hrs on computer with Intel(R) Core(TM) i5-2400 CPU @ 3.10GHz processor , i.e., on average 1.15 sec/rotation, and used 2.5 GB memory.


## **Benchmarking rotation scan:**

The script `RotationScan.py` tested on the known protein methylatransferase complex: Nucleosome-Dot1L (PDB: 6NJ9). 

The PDB: 6NJ9 contains ubiquitinated nucleosomes with two copies of Dot1L bound.
The substrate H3K79 is present in its active state, that can access the Dot1L reaction center. \

For benchmarking purpose, the ubiquitin and an extra copy of Dot1L were removed.
The H3Nle79 was mutated to lysine (H3K79), the location of NZ atom of H3K79 was considered as the *reaction center* and the *reaction center* was translated to the origin (0,0,0). Three random orientations of Dot1L were generated and these orientations were taken as initial models for *Rotation Scan*.

### **Benchmark runs with random orientations of DOT1L:**

```
cd Benchmarking_Dot1L-Nucleosome
python RotationScan_Benchmarking.py parameters_Dot1L_randomOrientation1.ini
python RotationScan_Benchmarking.py parameters_Dot1L_randomOrientation2.ini
python RotationScan_Benchmarking.py parameters_Dot1L_randomOrientation3.ini
```

## **Demo data:**
To demo `RotationScan.py`, the benchmark data set was used. For demo purpose, the *reaction center* was defined based on the coordinates of `CE` and `SD` atoms of the SAM instead of considering the `NZ` atom of H3K79 of nucleosome and the `NZ` atom was moved to the *reaction center*. The resolution of rotation scan was set to 45̊  to shorten the run time. 

### **How to run demo data**
```
cd demo_data
python ../RotationScan.py parameters_demo.ini
cd Results_PE__dot1l_nuc_res045_20230518_21h47m7s
python ../../plotting_RotationScan_results.py Results_PE__dot1l_nuc_res045_20230518_21h47m7s.data
```
**Note:** Demo data run time 201.2 sec


## **Rotation scan: Rv2067c-Nucleosome reaction complex model**

Accessiblity of Rv2067c's *reaction center* to nuclosomal H3K79 was tested using the enzyme-substrate *reaction complex* model between Rv2067c and nucleosome. Nucleosome coordinates were taken from PDB: 6NJ9. The model was constructed as described above.

### **How to run rotation scan on Rv2067c-Nucleosome complex model:**
#### **New run, results and plotting:**

```
cd Rv2067c-Nucleosome
./new_run1.sh
cd Results_PE__Rv2067c_Nuc_res010_20221211_2h8m56s
python ../../plotting_RotationScan_results.py Results_PE__Rv2067c_Nuc_res010_20221211_2h8m56s.data
```

### **Old run results (reported in the publication: https://www.biorxiv.org/content/10.1101/2023.02.24.529973v1)**

```
cd Rv2067c-Nucleosome/old_run

```
**Remark:** Old and new runs produce identical results as the input parameters are same. The script used in new run was more automated interms of processing of input pdb files and run parameters supplied with `parameters_RotationScan.ini` file.



