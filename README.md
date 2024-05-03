# HMSCPipeLine
A pipeline for running HMSC

This readme file is a work in progress

# Instuctions
## Data format
The code is written assuming that your data are stored in csv files, this is preferable for long term data storage as these files will remain readable. Alternatively, data can be stored as a RData file, with an object for each input matrix.
Matrix formats are as follows:

__ADD TABLE__

## HPC VS R
There are two methods for fitting the model, one outputs an RDS file which can be fitted using the Tensorflow based HPC fitting.

# File structure
There are multiple ways to structure how the outputs are saved. For the HPC approach, I currently make use of single model runs withing stored as follows:
- Hmsc outputs
  - <Model_Name>
    - Models
      - Fitted
      - INITS
      - Raw_HPC
      - Temp
      - Unfitted
    - Results
      - Preds

For R based fitting fewer output files are required and if you wish to set up multiple models as listed in the University of Helsinkis HMSC Course pipeline a completelly different file structure can be used.

# HPC
##General format
The Fit_HMSC_Model and Fit_HMSC_Model_KFold are set up so that they can be submitted using SLRUM as follows:

sbatch --export=ALL,SAMP=<samples>,THIN=<thin>,AREA=<Area>,NAME=<Model_name> Fit_HMSC_Model.sh

-Samples: The number of samples
-Thin: Model thinning
-Area: If you are formatting you files so that there is a folder for each area that contains the Hmsc outputs folder use this. Else feel free to remove it from you local version of the batch script
-Model name: The name of the model. This is the name of the colder which countains the models Models and Results subfolders.

NOTE: Samples and thin are only used for file matching, you must have ran and uploaded the output of S2_Model_Fitting_HPC_Version to the HPC storeage for it to run correctly 

Fuzzy matching works for job submition, personally I recomned Fit_*el.sh and Fit_*old.sh rather then typing out the full names

## File notes
The HPC scripts assume that your data on the HPC storage is in the same folder structure as described above, however, the unfitted folder is not required and the RAW_HPC folder is named sampled.
Sampled HPC files are merged on your local machine by S2b_HPC_output_merger, there is no need to retain the unmerged outputs locally.

Batch scripts should be stored in the folder that contains are the subfolders for each Area, the way they point to files expects this. There might still be some residual file locating errors in your local setup, please double check all pointers.
##R Scripts
S4a_Partition_Creation and S4b_HPC_Post_Processing are required to be stored on the HPC in the same file location as the batch scripts.

