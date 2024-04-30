#### Hmsc analyses on ####
#General cleaning of the workspace
remove(list=ls())
gc()

### Set up directories #### 

#If you are using RStudio this will set the working directory to exactly where the file is 
setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path)))

localDir = "./Hmsc Outputs"
dataDir = "./Data"

#If the data are saved as an RData file
load(file.path(dataDir,"Example4_Data.Rdata"))

#Else read in each matrix as a csv, be careful of full stops in the row and
#column names.

#Normally you would have to subset the varibales more then this, but because the
#example data is very clean those steps are skipped here and the study design
#matrix is defined
studyDesign = data.frame(Site = as.factor(Pi$Site))
xycoord = data.frame(x = Pi$x, y = Pi$y)
rownames(xycoord) = Pi$Site

rl.coord = HmscRandomLevel(sData = xycoord)
rm(Pi)

model = Hmsc(Y = Y, XData = X, XFormula = ~ x1, 
             TrData = Tr, phyloTree =  phy, TrFormula = ~Tr1, 
             studyDesign = studyDesign ,ranLevels = list(Site = rl.coord),
             distr = "normal")

model_description = "Example_x1_TRphy_Site"

#These lists are left over from the orginal pipeline which had multiple models
#defined at the same time
models = list(model)
names(models) = c("AB x1 Traits Tree Site")

#Check if the model with this description has been created before, if not create
#both the model and results directories and all sub folders
ModelDir = file.path(localDir, sprintf("%s",model_description))
save.dir = file.path(ModelDir, "Models")
resluts.dir = file.path(ModelDir, "Results")

if(!dir.exists(save.dir)){
  dir.create(ModelDir)
  dir.create(save.dir)
  for(x in c("Fitted","Raw_HPC","Unfitted","INITS","Temp")){
    dir.create(file.path(save.dir,x))
  }
  if(!dir.exists(resluts.dir)){
    dir.create(resluts.dir)
  }
} else {
  cat(sprintf("\n%1$s\nWARNING\n%1$s\nThis model has been ran before, unfitted model file being over written\n",strrep("-",10)))
}
save(models, file = file.path(save.dir, "/Unfitted/unfitted_models.RData"))