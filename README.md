# Self-assembly of amyloid-$\beta$ peptides with chitosan using WEPPROM and MARTINI forcefields

This repository source files, and protocols to perform self-assembly simulations similar to those presented in ["Mechanistic insights...," _Phys. Chem. Chem. Phys._ (2023)](placeholder_link).

## File Contents

* `ff/`: Forcefield files in GROMACS formats.
    * `forcefield.itp` : Definitions of coarse-grained beads used in this study, and their non-bonded interactions.  
    The special-bead types used in WEPPROM are patched with standard MARTINI types used in 
    * `abeta.itp` : Topology of the 7-residue chentral hydrophobic core region of amyloid-$\beta$, K<sub>16</sub>LVFFAE<sub>22</sub>, parameterized with the WEPPROM forcefield.  
    See [Sahoo et. al. "Pathways of amyloid-beta absorption and aggregation in a membranous environment" _Phys. Chem. Chem. Phys._ (2019)](https://par.nsf.gov/biblio/10174912) for details.
    * `posre_abeta.itp` : Restrains the position of residue F<sub>19</sub>'s BB bead of amyloid-$\beta$.
    * `NGLU/NGLU*.itp` : TBA
    * `water.itp` : MARTINI polarizable water, where bonds between the central bead and dummy charges are treated with the LINCS constraint algorithm. 
    * `water_em.itp` : MARTINI polarizable water, where bonds between the central bead and dummy charges are treated with harmonic bonds.  
    Only used during energy minimization in case of instabilities. 

* `generate_chitosan_itps/`: Scripts to generate itps of ensembles of chitosan molecules for a given number of chains, number of monomers per chain, and pH-dependent degree of cationicity.

* `protocols/`: MD parameter files for energy minimization (`em.mdp`), equilibration (`eq_1fs.mdp`, `eq.mdp`), production MD (`md.mdp`) and for ion insertion with genion (`ions.mdp`)

## Self-asssembly simulations of amyloid-$\beta$ with chitosan 

Simulation setup involves the following modules. First, create a new working directory for the simulation set up. We will call it `my_system/`
> `mkdir my_system/`

### I. **Prepare forcefield files for chitosan.**

1. Calculate the fraction of chitosan's 2'-amine groups that will be protonated/cationic at the selected pH using the Hendersen-Hasselbalch equation.

    $$ pH = pKa + log\frac{ concentration\ of\ conjugate\ base }{ concentration\ of\ conjugate\ acid } $$

    $$ R-NH_3^+ \longleftrightarrow R-NH_2 + H^+ $$

    For chitosan, $R-NH_3^+$ is the conjugate acid; $R-NH_2$ is the conjugate base, and the pKa for this reaction is 6.5.

    $$ pH = 6.5 + log\frac{[R-NH_2]}{[R-NH_3^+]} $$
    $$ \frac{[R-NH_2]}{[R-NH_3^+]} = 10^{pH - 6.5}$$
    $$ if\ \alpha = degree\ of\ cationicity/protonation $$ 
    $$ \frac{1-\alpha}{\alpha} = 10^{pH - 6.5}$$
    $$ {\alpha} = \frac{1}{1+10^{pH - 6.5}}$$

    For pH 7.5, $ {\alpha} = \frac{1}{1+10^{7.5 - 6.5}} = 0.0\bar{9} = 0.1 = 10\%$

2. **Distribute the cationic charge randomly across all the chitosan chains.** 

    Follow `generate_chitosan_itps/README.md` to create a set of itp files. In `generate_chitosan_itps/gen_itp.py`, change  `n_chain`, `n_monomer` and `prot_percentage` as  required.

    `generate_chitosan_itps/nglu/` contains the output files.

    Copy this directory to `my_system`
    > `cp -r generate_chitosan_itps/nglu/ my_system/`

### II. Initializing a random configuration of solutes.

1. Navigate to the `build_confs/` directory

2. Place 20 chitosan chains at random positions in a 18 nm cubic box.  
    **Input**: c1.gro, contains coordinates for 1 30-mer chitosan chain.  
    **Output**: c20.gro, coordinates of 100 peptides.  
    Using random seed 10 for reproducibility.

    > `gmx insert-molecules -ci c1.gro -nmol 20 -o c20.gro -seed 10 -box 18 18 18`

2. Place 100 peptides at random positions in a 18 nm cubic box.  
    **Input**: c20.gro (coordinates of 20 randomly placed chitosan chains ) a1.gro (coordinates of 1 unstructured amyloid beta peptide)
    **Output**: c20_a100.gro, coordinates of 20 randomly placed chitosan chains and 100 amyloid-beta peptides.  
    (Using random seed 10 for reproducibility.)
    > `~/programs/bin/gmx_mpi insert-molecules -f c20.gro -ci a1.gro -nmol 100 -o c20_a100.gro -box 18 18 18 -seed 10`

3. Solvate the box with MARTINI polarizable water.  
    Number of solvent molecules capped at 41243.
    **Input**: chitosan coordinates, `c20_a100.gro`, and water coordinates `polarize-water.gro.`
    **Output**: `c20_a100_solv.gro` : Solvated box containing 100 peptides.
    > `gmx solvate -cp c20_a100.gro -cs ../polarize-water.gro -maxsol 41243 -o c20_a100_solv.gro`

4. Create a new directory for further processing and simulation.
    
    > `cp c20_100_solv.gro ../my_system/`

### III. Create an initial topology file (.top) for the solvated system

Here is an example file. Modify paths, molecule names, and numbers of molecules as needed, and save as `system.top` in `my_system/`

```
; Include the main forcefield file
#include "../ff/ChitoABff.itp"

; Include itps for chitosan molecules
#include "./nglu/nglu_0000.itp"
#include "./nglu_0001.itp"

; Include itps for water
#ifdef FLEXIBLE
#include "../ff/water_em.itp"
#else
#include "../ff/water.itp"

; Include itps for monovalent ions
#include "../ff/ions.itp"

[ system ] 
Chitosan-AB self-assembly

[ molecules ]
; molecule_name    number_of_molecules
; For n chitosan (NGLU) molecules, list molecule names as NGLU<nglu_id>,
; where 0 <= nglu_id < n, and nglu has 4 digits with padding zeros.
; Example for 2 chitosan chains
NGLU0000            1                   
NGLU0001            1
Protein             100
PW                  41243

```

### IV. Neutralize excess charges with monovalent ions

1. Navigate to `my_system/`

1. Create an dummy tpr for ion addition.
    > `~/programs/bin/gmx_mpi grompp -f ../protocols/ions.mdp -c c20_a100_solv.gro -r c20_a100_solv.gro -p system.top -o ions.tpr`

3. Insert ions. Select group number corresponding to polarizable water. In this case, select Group 13: PW.
    > `~/programs/bin/gmx_mpi genion -s ions.tpr -p system.top -neutral -seed 10 -o c20_a100_ions.gro`

### V. Energy minimization
Create a `em ` directory, `em.tpr` file.
> `mkdir em; ~/programs/bin/gmx_mpi grompp -f ../protocols/em.mdp -c c20_a100_ions.gro -r c20_a100_ions.gro -p system.top -o em/em`

Perform energy minimization.
> `~/programs/bin/gmx_mpi mdrun -s em/em.tpr -deffnm em/em -v`

### VI. Equilibration

We will set up 2 replicas, each initialized at different velocities.

1. Make sub-directories for each iteration.
>> `mkdir iter1 iter2`

2. Create an equilibration tpr for each iteration.

We use Parrinello-Rahman barostat for equilibration, which raises a warning in GROMACS. The warning will be suppressed with the `-maxwarn 1` option in the `gmx grompp` call.

> `mkdir iter1/eq/; ~/programs/bin/gmx_mpi grompp -f ../protocols/eq_1fs.mdp -c em/em.gro -r em/em.gro -p system.top -o iter1/eq/eq -maxwarn 1`
> `mkdir iter2/eq/; ~/programs/bin/gmx_mpi grompp -f ../protocols/eq_1fs.mdp -c em/em.gro -r em/em.gro -p system.top -o iter2/eq/eq -maxwarn 1`

3. Perform equilibration.
Run on an HPCC for faster performance. If using non-mpi installation of gromacs, try the `-nt` option for parallelization with multithreading.
> `~/programs/bin/gmx_mpi mdrun -s iter1/eq/eq.tpr -deffnm iter1/eq/eq -v`
>
> `~/programs/bin/gmx_mpi mdrun -s iter2/eq/eq.tpr -deffnm iter2/eq/eq -v`

### VII. Production MD
Perform 700 ns of unrestrained NPT simulation

1. Create tprs:
> `mkdir iter1/md/; ~/programs/bin/gmx_mpi grompp -f ../../protocols/md.mdp -c iter1/eq/eq.gro -r iter1/eq/eq.gro -t iter1/eq/eq.cpt -p system.top -o iter1/md/md -n index.ndx`
> `mkdir iter2/md/; ~/programs/bin/gmx_mpi grompp -f ../../protocols/md.mdp -c iter2/eq/eq.gro -r iter2/eq/eq.gro -t iter2/eq/eq.cpt -p system.top -o iter2/md/md -n index.ndx`

2. Perform production MD. Running on HPCC is strongly recommended. For testing, the user may use the `-nsteps` option to stop the simulation after a certain number of steps, although output options in md.mdp should be edited to output files at smaller intervals.
