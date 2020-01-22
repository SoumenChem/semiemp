# SEMIEMP

This code can be used as an interface to NWChem. It can be used to perform real-time dynamics using INDO/S Hamiltonian. For sample input see the manual. 

To set up the code-

(i) Move "semiemp" folder to "nwchem/src/".

(ii) Copy task.F, task_input.F, task_energy, task_gradient.F files from "nwchem/src/semiemp/setup" folder to "nwchem/src/task/" folder.

(iii) Copy make_nwchem_config from "nwchem/src/semiemp/setup" folder to  "nwchem/src/config" folder.

(iv) Copy util_module_avail.F from "nwchem/src/semiemp/setup" folder to "nwchem/src/util" folder.

(v) Compile nwchem. (It might be necessary to run "make realclean" before compiling the code.




Please cite our original publications when you use Semiemp-

Soumen Ghosh, Jason C. Asher, Laura Gagliardi, Christopher J. Cramer, Niranjan Govind. A Semiempirical Effective Hamiltonian Based Approach for Analyzing Excited State Wave Functions and Computing Excited State Absorption Spectra Using Real-Time Dynamics. J. Chem. Phys. 2019, 150, 104103 .

Soumen Ghosh, Amity Andersen, Laura Gagliardi, Christopher J. Cramer and Niranjan Govind. Modeling Optical Spectra of Large Organic Systems using Real-Time Propagation of Semiempirical Effective Hamiltonians. J. Chem. Theory Comput. 2017, 13, 4410−4420.


