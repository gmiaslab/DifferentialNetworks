# Applying Differential Network Analysis toLongitudinal Gene Expression in Response to Perturbations
This repository implements the research paper by Xue et al., where the data, code, and supplementary materials can all be found.

### Requirements and Dependencies
The code has been executed in Python 3.8 (and there may be incompatibilities with different Python versions). The Wolfram Language was used to generate the network plots, for better visualization of - to obtain the plots shown in this paper Mathematica (version 12.0.0+) must be installed. Other dependencies include: 
* Matplotlib >= 3.3.2
* Seaborn >= 0.11.0
* pandas = 1.1.3
* NetworkX >= 2.5

The saliva network computation can be executed on a PC with a 2.50 GHz processor and 6.00 GB RAM, taking no more than 15 minutes. The B cell network computation needs at least 30 GB RAM. 40 GB RAM is recommended for smooth execution. If any script takes longer than an hour on a PC, consider running it on high performance computing (HPC) clusters. The file, Bcell_slurm.sb, is an example SLURM queue submission script for the HPC clusters at Michigan State University's High Performance Computing Center. This file can be referenced and customized according to your computational resources (if using SLURM on your cluster). 

### Directories and Files
* SLV_data & Bcell_data: the input folders containing the experimental data. The first Column is the GENECODE identifier. The remaining colums are normalized gene expression at different time points. For the saliva data, the file "TimesHourly.csv" is for the time point (24 hour sampling hourly). For B cells, the columns listed in a format T##, U##, RT##, and RU##, where T means treated, U meams untreated, R means repeated experiment. The ## are the times corresponding to 0, 1, 2, 4, 7 and 15 hours after the start of the experiment. 
* Results: the supplementary results for the paper. The SLV_results and Bcell_results each have complete results for the corresponding datasets in the paper  
* SLV_results / Bcell_results: the resulting directories created after the execution of the `SLV_script.py` and the `Bcell_script.py`

### Reproducing the paper results
To produce the results as in the paper, run the script `SLV_script.py` and `Bcell_script.py`, where each generates one subdirectory ending with "_results" (i.e., SLV_results and Bcell_results). The outputs in these two folders should match the contents in the Manuscript_online_file with the same names.

For code testing on a PC (where the hardware specifcation is limited), adapt as follows:
1. In the `network_plot_by_mathematica.py` module (in line 26), set HPCC = False.
2. In either the `Bcell_script.py` or `SLV_script.py` script, uncomment line 23.
3. In the `dir_conf.py`, set the significant_comms (in line 37 or 40) to match the range of line 2. This helps to avoid the `list index out of range` error. 

For more details, please refer to the comments within the code. The results from the above modifications will differ from those in the paper, but the format will remain the same. 

## Funding
This study was funded by the Translational Research Institute for Space Health through National Aeronautics and Space Administration (NASA) Cooperative Agreement NNX16AO69A (Project Number T0412, PI: Prof. George Mias). 

## Contact Information
* Code Contributor/First Author: Shuyue Xue (xueshuy1@msu.edu)
* Project Principal Investigators (PIs): Professor George Mias (gmias@msu.edu) and Professor Carlo Piermarocchi (piermaro@msu.edu)