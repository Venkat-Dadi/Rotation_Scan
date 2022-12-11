# **Rotation Scan**

To check the accessibility of the *reaction center* of an enzyme to the its substrate by modeling the enzyme-substrate *reaction complex* for protein methylatransferases.

**Reaction complex:** a complex formed between an enzyme and its substrate during an enzymatic reaction. \
**Reaction center:** location of the methyl acceptor atom (NZ atom of lysine, in case of protein lysine methyltransferases) the substrate during the methyl group transfer. \

During protein lysine methylation, a *reaction complex* must be formed between the methyltransferase and its protein substrate (true for any enzymatic reaction).  The accessibility of the enzyme's *reaction center* to the methyl acceptor substrate depends on the milieus of the *reaction center* and the acceptor atom of the substrate. The accessibility can be calculated by modeling the enzyme-substrate *reaction complex*. In such complex, during the methyl transfer, the methyl acceptor atom of the substrate (NZ atom of lysine in case of lysine methyltransferases) lies at the *reaction center* in contact with methyl carbon atom (CE) of S-adenosyl-L-methionine. The substrate must approach the *reaction cetner* such that the binding pose between the enzyme and its substrate has no steric hidrance. The pose with minimum number of steric clashes can be calculated by constructing the enzyme-substrate *reaction complex* by imposing the criteria of the *reaction center* and rotating either the enzyme or substrate about the *reaction center* while the other is fixed immobile. Rotations were carried out with the help of Euler's rotation theorem, that covers all the possible orientation that a substrate can approach the *reaction center* of an enzyme. All these rotations constitutes the *Rotation Scan*. For each rotation, the number of atoms invovled in the clashes are calculated according to the clash criteria used by Probe ([DOI: 10.1006/jmbi.1998.2400](https://www.sciencedirect.com/science/article/pii/S0022283698924007?via%3Dihub)). The pose with minimum number of atomic clashes is considered as the possible binding mode of the substrate to the enzyme. This is implemented using python script: `RotationScan.py`.

Each pose of a *Rotation Scan* is represented by three elemental roations (α, β and γ) that constitute a single rotation. The total number of clashing atoms (TOT_CA or TCA) are plotted as a funciton of these elemental roations. For simplicity, two elemental rotation combinations were plotted while the third corresponds to the minimum TCA value.

**Caveats:**
* This model does not consider the flexibility either in the substrate or the enzyme.
* This model can only be used to show numberically that the inaccessibility of the active stie to the bulky substrates due to the steric hindrance.

## **Linux: install required packages**

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh 
conda activate base
pip install numpy
pip install pandas
pip install scipy
pip install tqdm
pip install matplotlib
```
**Note:** Windows OS with Miniconda and other python packages (numpy, pandas, scipy, tqdm) installed can use `RotationScan.py` script. \

## **Usage**
```
# Run rotation scan
python RotationScan.py parameters_RotationScan.ini

# Plotting rotation scan results
python plotting_RotationScan_results.py file_name.data
```

## **Benchmarking**

The script `RotationScan.py` tested on the known protein methylatransferase complex: Nucleosome-Dot1L (PDB: 6NJ9). 

The PDB: 6NJ9 contains ubiquitinated nucleosomes with two copies of Dot1L bound.
The substrate H3K79 is present in its active state, that can access the Dot1L reaction center. \

For benchmarking purpose, the ubiquitin and an extra copy of Dot1L were removed.
The H3Nle79 was mutated to lysine (H3K79), the location of NZ atom of H3K79 was considered as the *reaction center* and the *reaction center* was translated to the origin (0,0,0). Three random orientations of Dot1L were generated and these orientations were taken as initial models for *Rotation Scan*.

### **Benchmark runs with random orientations of Dot1L**

```
cd Benchmarking_Dot1L-Nucleosome
python RotationScan_Benchmarking.py parameters_Dot1L_randomOrientation1.ini
python RotationScan_Benchmarking.py parameters_Dot1L_randomOrientation2.ini
python RotationScan_Benchmarking.py parameters_Dot1L_randomOrientation3.ini
```
## **Rv2067c-Nucleosome reaction complex model**

Accessiblity of Rv2067c's *reaction center* to nuclosomal H3K79 was tested using the enzyme-substrate *reaction complex* model between Rv2067c and nucleosome. Nucleosome coordinates were taken from PDB: 6NJ9. The model was constructed as described above.

### **Rv2067c-Nucleosome Rotation Scan**
#### **Old run (results reported in the publication)**

```
cd Rv2067c-Nucleosome/old_run

```
#### **New run**
```
# new run 
cd Rv2067c-Nucleosome
./new_run1.sh

# new run results
cd Rv2067c-Nucleosome/Results_PE__Rv2067c_Nuc_res010_20221211_2h8m56s

# plotting rotation scan results
python ../../plotting_RotationScan_results.py Results_PE__Rv2067c_Nuc_res010_20221211_2h8m56s.data
```

