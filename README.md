# TCR-MHC-CD4_kinetic_modeling
Mathematical modeling for the study "Cooperative binding of TCR and CD4 to pMHC enhances TCR sensitivity"

The code provided encompasses modeling of: 
1. time dependent 2D adhesion frequency
2. zero-force kinetics from thermal fluctuation assay
3. force-dependent kinetics from BFP.

Software requirements are:

    1. Matlab 2020a 
    2. Julia 1.6

The code was tested on Windows 10, in the Matlab 2020a IDE and the VSCode v1.63.2 environment.

Installation instructions:
    1. clone this repo
    2. run analysis scripts in "/TCR-MHC-CD4 modeling/" subfolder. E.g to run adhesion frequency assay modeling, run the script "/TCR-MHC-CD4 modeling/TCR_CD4_MHC_adhesion_freq_fit.jl"
    3. any missing required depedencies will be prompted to be installed after running each script

Note: a recurring and simple abbreviation system used throughout the software for convenience:

    TCR -> T
    MHC -> M
    CD4 -> C

such that, for example, 'TM' and 'TMC' indicate TCR-MHC bimolecular and TCR-MHC-CD4 trimolecular bonds, respectively

Datasets used in the analysis are provided in the respective subfolders
