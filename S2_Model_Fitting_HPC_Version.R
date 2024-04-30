remove(list=ls())
set.seed(369)
require(Hmsc)
require(jsonify)

### Set up directories #### 

#If you are using RStudio this will set the working directory to exactly where the file is 
setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path)))
model_description = "Example_x1_TRphy_Site"
localDir = sprintf("./HMSC/Hmsc Outputs/%s",model_description)
ModelDir = file.path(localDir, "Models")

### Read in the unfitted models ####
load(file = file.path(ModelDir, "/Unfitted/unfitted_models.RData"))

samples_list = c(50, 100, 250, 250, 500, 500, 500, 750)
thin_list = c(10, 10, 10, 20, 20, 50, 60, 50)
nChains = 4
nParallel = 4
Lst = 1
verbose = 1

while(Lst <= length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  transient = ceiling(0.5*samples*thin)
  
  filename = file.path(ModelDir,sprintf("INITS/HPC_INIT_samples_%.4d_thin_%.2d_chains_%.1d.rds",samples,thin,nChains))
  m = sampleMcmc(models[[1]], samples = samples, thin=thin,
                 #adaptNf=rep(transient,models[[1]]$nr), 
                 transient = transient,
                 nChains = nChains,
                 verbose = verbose, 
                 engine = "HPC")
  
  saveRDS(to_json(m), filename)
  Lst = Lst + 1
}