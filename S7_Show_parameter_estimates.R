remove(list=ls())
require(Hmsc)
require(colorspace)
require(corrplot)
require(writexl)

set.seed(369)

### Set up directories #### Because I run this on two difference computers this
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

model_description = "Example_x1_TRphy_Site"
localDir = sprintf("./Hmsc Outputs/%s",model_description)
ModelDir = file.path(localDir, "Models")
UnfittedDir = file.path(ModelDir, "Unfitted")
ResultDir = file.path(localDir, "Results")

support.level.beta = 0.95
support.level.gamma =  0.95
support.level.omega =  0.9
var.part.order.explained = NULL #Default: in variance partitioning of explained variance, species are shown in the order they are in the model
var.part.order.raw = NULL #Default: in variance partitioning of raw variance, species are shown in the order they are in the model
show.sp.names.beta = NULL #Default: species names shown in beta plot if there are at most 30 species and no phylogeny
plotTree = NULL #Default: tree is plotted in Beta plot if the model includes it
omega.order = NULL #Default: species shown in the order they are in the model
show.sp.names.omega = NULL #Default: species names shown in beta plot if there are at most 30 species

samples_list = c(50, 100, 250, 250, 500, 500, 500, 500, 500, 500)
thin_list = c(10, 10, 10, 20, 20, 30, 40, 50, 60, 75)
nst = length(thin_list)
nChains = 4

text.file = file.path(ResultDir,paste0(model_description,"parameter_estimates.txt"))
cat(c("This file contains additional information regarding parameter estimates.\n\n"),file=text.file)

load(file = file.path(UnfittedDir, "unfitted_models.RData"))

#Note: 20240223
#I should move this to a separate file, save the output and remove the individual chain files
for (Lst in nst:1) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  chains = 4
  #Note that I use different file names for the R fitted and HPC fitted models
  #just to keep track
  
  #filename = file.path(ModelDir,sprintf("Fitted/HPC_samples_%.4d_thin_%.2d_chains_%.1d.Rdata",samples, thin, nChains))
  filename = file.path(ModelDir,sprintf("Fitted/FittedR_samples_%.4d_thin_%.2d_chains_%.1d.Rdata", samples, thin, nChains))
  if(file.exists(filename)){
    cli_alert_success("File {filename} exists")
    break} else{
      cli_alert_danger("File {filename} exists")
      cli_alert_info(cli_par("Run computing model fit script for:\n Thin: {thin} \t Samples: {samples} \t Chains: {nChains}"))
      cli_rule()
    }
}

if(file.exists(filename)){
  load(filename)
  #If you are using R fitted models you don't need to run the following two lines as the model is saved differently.
  #m = fitted_model$posteriors
  #rm(fitted_model)
  cat(sprintf("Model Description\nSamples: %.4d\tThin: %.3d\tChains: %.1d\n\n", samples, thin, nChains),file=text.file,sep="",append=TRUE)
  if(is.null(var.part.order.explained)){
    var.part.order.explained = list()
    var.part.order.explained = 0
  }
  if(is.null(var.part.order.raw)){
    var.part.order.raw = list()
    var.part.order.raw = 0
  }
  if(is.null(omega.order)){
    omega.order = list()
    omega.order = 0
  }
  
  pdf(file= file.path(ResultDir,paste0(model_description,"parameter_estimates.pdf")))
  
  if(m$XFormula=="~."){
    covariates = colnames(m$X)[-1]
  } else {
    covariates = attr(terms(m$XFormula),"term.labels")
  }
  if(m$nr+length(covariates)>1 & m$ns>1){
    preds = computePredictedValues(m)
    VP = computeVariancePartitioning(m)
    vals = VP$vals
    mycols = rainbow(nrow(VP$vals))
    MF = evaluateModelFit(hM=m, predY=preds)
    R2 = NULL
    if(!is.null(MF$TjurR2)){
      TjurR2 = MF$TjurR2
      vals = rbind(vals,TjurR2)
      R2=TjurR2
    }
    if(!is.null(MF$R2)){
      R2=MF$R2
      vals = rbind(vals,R2)
    }
    if(!is.null(MF$SR2)){
      R2=MF$SR2
      vals = rbind(vals,R2)
    }
    filename = file.path(ResultDir, paste0(model_description,"parameter_estimates_VP_.csv"))
    write.csv(vals,file=filename)
    if(!is.null(VP$R2T$Beta)){
      filename = file.path(ResultDir,paste0(model_description,"parameter_estimates_VP_R2T_Beta.csv"))
      write.csv(VP$R2T$Beta,file=filename)
    }
    if(!is.null(VP$R2T$Y)){
      filename = file.path(ResultDir, paste0(model_description,"parameter_estimates_VP_R2T_Y.csv"))
      write.csv(VP$R2T$Y,file=filename)
    }
    if(all(var.part.order.explained==0)){
      c.var.part.order.explained = 1:m$ns
    } else {
      if(all(var.part.order.explained=="decreasing")){
        c.var.part.order.explained = order(R2, decreasing = TRUE)
      }
      else {
        c.var.part.order.explained  = var.part.order.explained
      }
    }
    VP$vals = VP$vals[,c.var.part.order.explained]
    cat(c("\n","var.part.order.explained","\n","\n"),file=text.file,sep="",append=TRUE)
    cat(m$spNames[c.var.part.order.explained],file=text.file,sep="\n",append=TRUE)
    plotVariancePartitioning(hM=m, VP=VP, main = paste0("Proportion of explained variance, ",model_description), cex.main=0.8, cols = mycols, args.leg=list(bg="white",cex=0.7))
    if(all(var.part.order.raw==0)){
      c.var.part.order.raw = 1:m$ns
    } else {
      if(all(var.part.order.raw=="decreasing")){
        c.var.part.order.raw = order(R2, decreasing = TRUE)
      }
      else {
        c.var.part.order.raw  = var.part.order.raw
      }
    }
    if(!is.null(R2)){
      VPr = VP
      for(k in 1:m$ns){
        VPr$vals[,k] = R2[k]*VPr$vals[,k]
      }
      VPr$vals = VPr$vals[,c.var.part.order.raw]
      cat(c("\n","var.part.order.raw","\n","\n"),file=text.file,sep="",append=TRUE)
      cat(m$spNames[c.var.part.order.raw],file=text.file,sep="\n",append=TRUE)
      plotVariancePartitioning(hM=m, VP=VPr,main=paste0("Proportion of raw variance, ",model_description),cex.main=0.8, cols = mycols, args.leg=list(bg="white",cex=0.7),ylim=c(0,1))
    }
  }
  
  if(m$nc>1){
    postBeta = getPostEstimate(m, parName="Beta")
    filename = file.path(ResultDir, paste0(model_description,"parameter_estimates_Beta_.xlsx"))
    me = as.data.frame(t(postBeta$mean))
    me = cbind(m$spNames,me)
    colnames(me) = c("Species",m$covNames)
    po = as.data.frame(t(postBeta$support))
    po = cbind(m$spNames,po)
    colnames(po) = c("Species",m$covNames)
    ne = as.data.frame(t(postBeta$supportNeg))
    ne = cbind(m$spNames,ne)
    colnames(ne) = c("Species",m$covNames)
    vals = list("Posterior mean"=me,"Pr(x>0)"=po,"Pr(x<0)"=ne)
    writexl::write_xlsx(vals,path = filename)
    if(is.null(show.sp.names.beta)){
      c.show.sp.names = (is.null(m$phyloTree) && m$ns<=30) 
    } else {
      c.show.sp.names = show.sp.names.beta
    }
    c.plotTree = !is.null(m$phyloTree)
    if(!is.null(plotTree)){
      c.plotTree = c.plotTree & plotTree
    }
    plotBeta(m, post=postBeta, supportLevel = support.level.beta, param="Sign",
             plotTree = c.plotTree,
             covNamesNumbers = c(TRUE,FALSE),
             spNamesNumbers=c(c.show.sp.names,FALSE),
             cex=c(0.6,0.6,0.8))
    mymain = paste0("BetaPlot, ",model_description)
    if(!is.null(m$phyloTree)){
      mpost = convertToCodaObject(m)
      rhovals = unlist(poolMcmcChains(mpost$Rho))
      mymain = paste0(mymain,", E[rho] = ",round(mean(rhovals),2),", Pr[rho>0] = ",round(mean(rhovals>0),2))
    }
    title(main=mymain, line=2.5, cex.main=0.8)
  }
  
  
  if(m$nt>1 & m$nc>1){
    postGamma = getPostEstimate(m, parName="Gamma")
    plotGamma(m, post=postGamma, supportLevel = support.level.gamma, param="Sign",
              covNamesNumbers = c(TRUE,FALSE),
              trNamesNumbers=c(m$nt<21,FALSE),
              cex=c(0.6,0.6,0.8))
    title(main=paste0("GammaPlot ",model_description), line=2.5,cex.main=0.8)
  }
  
  #Compute Species Associations for the models

  if(m$nr>0 & m$ns>1){
    OmegaCor = computeAssociations(m)
    for (r in 1:m$nr){
      toPlot = ((OmegaCor[[r]]$support>support.level.omega) + (OmegaCor[[r]]$support<(1-support.level.omega))>0)*sign(OmegaCor[[r]]$mean)
      if(is.null(show.sp.names.omega)){
        c.show.sp.names = (m$ns<=30) 
      } else {
        c.show.sp.names = show.sp.names.omega
      }
      if(!c.show.sp.names){
        colnames(toPlot)=rep("",m$ns)
        rownames(toPlot)=rep("",m$ns)
      }
      if(all(omega.order==0)){
        plotOrder = 1:m$ns
      } else {
        if(all(omega.order=="AOE")){
          plotOrder = corrMatOrder(OmegaCor[[r]]$mean,order="AOE")
        } else {
          plotOrder = omega.order
        }
      }
      cat(c("\n","omega.order","\n","\n"),file=text.file,sep="",append=TRUE)
      cat(m$spNames[plotOrder],file=text.file,sep="\n",append=TRUE)
      mymain = paste0("Associations, ",model_description, ": ",names(m$ranLevels)[[r]])
      if(m$ranLevels[[r]]$sDim>0){
        mpost = convertToCodaObject(m)
        alphavals = unlist(poolMcmcChains(mpost$Alpha[[r]][,1]))
        mymain = paste0(mymain,", E[alpha",r,"] = ",round(mean(alphavals),2),", Pr[alpha",r,">0] = ",round(mean(alphavals>0),2))
      }
      corrplot(toPlot[plotOrder,plotOrder], method = "color",
               type = "lower",
               col=colorRampPalette(c("blue","white","red"))(3),
               mar=c(0,0,1,0),
               main=mymain,cex.main=0.8)
      me = as.data.frame(OmegaCor[[r]]$mean)
      me = cbind(m$spNames,me)
      colnames(me)[1] = ""
      po = as.data.frame(OmegaCor[[r]]$support)
      po = cbind(m$spNames,po)
      colnames(po)[1] = ""
      ne = as.data.frame(1-OmegaCor[[r]]$support)
      ne = cbind(m$spNames,ne)
      colnames(ne)[1] = ""
      vals = list("Posterior mean"=me,"Pr(x>0)"=po,"Pr(x<0)"=ne)
      filename = file.path(ResultDir, paste0("parameter_estimates_Omega_",names(m$ranLevels)[[r]],".xlsx"))
      writexl::write_xlsx(vals,path = filename)
    }
  }
  while(!is.null(dev.list())){
    dev.off()
  }
}
