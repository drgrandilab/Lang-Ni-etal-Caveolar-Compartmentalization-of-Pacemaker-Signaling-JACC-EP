# Mouse 3D SANC Model

Please cite: 

Lang, D.，Ni, H.，et al. Caveolar Compartmentalization of Pacemaker Signaling Ensures Stable Sinoatrial Rhythmicity Which Is Disrupted in Heart Failure. JACC: Clinical Electrophysiology https://doi.org/10.1016/j.jacep.2026.01.003 (2026) doi:10.1016/j.jacep.2026.01.003.


Code Author: Haibo Ni <haibo.ni02@gmail.com>

*** last update Dec 25 2025

### Requires
* Intel Compiler: Intel® oneAPI DPC++/C++ Compiler
* Operating System: Tested with ubuntu 25.04
* Python scripts: Pandas, numpy, scipy, matplotlib  

### To run the code, run 
```./SAN_Pace```
or submit through slurm system
```sbatch run.slurm```


### Sample output 
in ```Sample_ouput``` folder
```bash
tar -xvzf Sample_output.tar.gz
```

### Plot by 

run ```bash
python check_line_data.py```


<img width="778" height="705" alt="image" src="https://github.com/user-attachments/assets/1c0282ac-bacb-468f-8c97-ce85c231568e" />
