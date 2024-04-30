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
