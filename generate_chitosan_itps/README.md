Generates a series of chitosan itp files with randomzied protonation sites

Usage:

- Generates random numbers for protonation sites

python gen_proto.py

- Uses generated random number to modify template system.protein.itp file into randomized protonation itp file

python gen_itp.py

In this demo, number of chitosan chains was set to 40 (n_nglu = 40, in gen_proto.py) and the number of chitosan monomers in each chain was set to 30 (n_mono = 30, in gen_proto.py). As a result, a temporary file named proto.txt was generated with 30x40 numbers (0 represents unprotonated, 1 represents protonated).

The gen_itp.py reads from the generated proto.txt file for randomized protonated site, and reads from system.protein.itp for 30-mer chitosan itp template. It applied modfications of B3 bead protonation state (charge and bead type changes) to the template itp file (system.protein.itp) and generates 40 random chains of chitosan with different protonation states. The overall protonation rate is constrained, but for each individual chain the protonation rate may vary. 
