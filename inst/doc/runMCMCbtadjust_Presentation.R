## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(warnings = FALSE)

### first calculates conditional (logical) variables to activate - or not - child documents whether the environment is OK to build this vignette
condition_rstan<-TRUE
if (nchar(system.file(package='rstan'))==0) {condition_rstan<-FALSE}
condition_greta<-TRUE
if (nchar(system.file(package='greta'))==0) {condition_greta<-FALSE}
	  if (nchar(system.file(package='R6'))==0) {condition_greta<-FALSE}
	  if (nchar(system.file(package='tensorflow'))==0) {condition_greta<-FALSE}
	#removed because interfers with Rmarkdown.  test<-capture.output(greta::greta_sitrep(),type="message")
	# if (max(regexpr("greta is ready to use",test))<=0) {stop("greta should be ready to use, which is not the case - call \"greta_sitrep()\" for more details")}
  if (condition_greta) {temp<-try(greta::normal(0,1)); if (class(temp)[1]=="try-error") {condition_greta<-FALSE}}
condition_jags<-TRUE                             
	  if (nchar(system.file(package='rjags'))==0) {condition_jags<-FALSE}
	  if (nchar(system.file(package='runjags'))==0) {condition_jags<-FALSE}
	  
if (condition_jags) 
  {suppressWarnings(temp<-runjags::testjags(silent=TRUE))
	 if(!(temp$JAGS.available&temp$JAGS.found&temp$JAGS.major==4)) {condition_jags<-FALSE}}

condition_nimble<-TRUE   
	  if (nchar(system.file(package='nimble'))==0) {condition_nimble<-FALSE}



## ----Simulating data----------------------------------------------------------

set.seed(1)
y1000<-rnorm(n=1000,mean=600,sd=30)


## ----Nimble setup-------------------------------------------------------------

library(runMCMCbtadjust)
library(nimble)


## ----Nimble Data--------------------------------------------------------------

ModelData <-list(mass = y1000)
ModelConsts <- list(nobs = length(y1000))


## ----Nimble Code--------------------------------------------------------------
 ModelCode<-nimbleCode(
  {
    # Priors
    population.mean ~ dunif(0,5000)
    population.sd ~ dunif(0,100)
    
    # Normal distribution parameterized by precision = 1/variance in Nimble
    population.variance <- population.sd * population.sd
    precision <- 1 / population.variance
  
    # Likelihood
    for(i in 1:nobs){
      mass[i] ~ dnorm(population.mean, precision)
    }
  })


## ----Nimble Inits-------------------------------------------------------------
ModelInits <- function()
{list (population.mean = rnorm(1,600,90), population.sd = runif(1, 1, 30))}
  
Nchains <- 3

set.seed(1)
Inits<-lapply(1:Nchains,function(x){ModelInits()})

#specifying the names of parameters to analyse and save:
params <- c("population.mean", "population.sd") 

#devising the maximum level allowed for the Geweke diagnostic of convergence (cf. following)
npars<-length(params)
Gew.Max<-as.double(format(quantile(sapply(1:100000,function(x,N){max(abs(rnorm(N)))},npars),0.95),digits=3,scientific=FALSE))


## ----Nimble first run, cache=FALSE--------------------------------------------
out.mcmc.Coda.Geweke<-runMCMC_btadjust(code=ModelCode, constants = ModelConsts, data = ModelData, MCMC_language="Nimble",
    Nchains=1, params=params, inits=Inits[1],
    niter.min=1000, niter.max=300000,
    nburnin.min=100, nburnin.max=200000, 
    thin.min=1, thin.max=1000,
    conv.max=Gew.Max, neff.min=1000,
    control=list(neff.method="Coda", convtype="Geweke"))


## ----Nimble first run with print.diagnostics, cache=FALSE---------------------
out.mcmc.Coda.Geweke.with.print.diagnostics<-runMCMC_btadjust(code=ModelCode, constants = ModelConsts, data = ModelData, MCMC_language="Nimble",
    Nchains=1, params=params, inits=Inits[1],
    niter.min=1000, niter.max=300000,
    nburnin.min=100, nburnin.max=200000, 
    thin.min=1, thin.max=1000,
    conv.max=Gew.Max, neff.min=1000,
    control=list(neff.method="Coda", convtype="Geweke",print.diagnostics=TRUE))


## ----nature of the output, cache=FALSE----------------------------------------

length(out.mcmc.Coda.Geweke)

class(out.mcmc.Coda.Geweke)

head(out.mcmc.Coda.Geweke[[1]])



## ----attributes of the output, cache=FALSE------------------------------------

names(attributes(out.mcmc.Coda.Geweke))

names(attributes(out.mcmc.Coda.Geweke)$package.versions)

attributes(out.mcmc.Coda.Geweke)$final.params

attributes(out.mcmc.Coda.Geweke)$sessionInfo


## ----Nimble second run, cache=FALSE-------------------------------------------
out.mcmc<-runMCMC_btadjust(code=ModelCode, constants = ModelConsts, data = ModelData, MCMC_language="Nimble",
    Nchains=Nchains, params=params, inits=Inits,
    niter.min=1000, niter.max=300000,
    nburnin.min=100, nburnin.max=200000, 
    thin.min=1, thin.max=1000,
    conv.max=1.05, neff.min=1000)


## ----Nimble comparison first two runs,echo=FALSE------------------------------
compar_MCMC_times<-format(cbind(unlist(attributes(out.mcmc.Coda.Geweke)$final.params),unlist(attributes(out.mcmc)$final.params)),digits=4,scientific=FALSE)


## ----Nimble comparison first two runs invisible,echo=FALSE--------------------


knitr::kable(compar_MCMC_times[1:(dim(compar_MCMC_times)[1]-1),],col.names=c("Nimble.Coda.Geweke","Nimble.default"),align="r",caption=paste0("Comparison of the efficiency of first two NIMBLE models:"))


## ----Nimble third run, cache=FALSE--------------------------------------------

out.mcmc.Geweke<-runMCMC_btadjust(code=ModelCode, constants = ModelConsts, data = ModelData, MCMC_language="Nimble",
    Nchains=1, params=params, inits=Inits[1],
    niter.min=1000, niter.max=300000,
    nburnin.min=100, nburnin.max=200000, 
    thin.min=1, thin.max=1000,
    conv.max=Gew.Max, neff.min=1000,
    control=list(convtype="Geweke"))


## ----Nimble comparison first three runs, echo=FALSE---------------------------
compar_MCMC_times<-format(cbind(unlist(attributes(out.mcmc.Coda.Geweke)$final.params),unlist(attributes(out.mcmc.Geweke)$final.params),unlist(attributes(out.mcmc)$final.params)),digits=4,scientific=FALSE)


knitr::kable(compar_MCMC_times[1:(dim(compar_MCMC_times)[1]-1),],col.names=c("Nimble.Coda.Geweke","Nimble.Geweke","Nimble.default"),align="r", caption="Comparison of the efficiency of the three NIMBLE models:")


## ----Nimble comparison parameters of first three runs, warning=FALSE,echo=FALSE----
compar_MCMC_params<-format(cbind(c(ks.test(unlist(out.mcmc[,"population.mean"]),unlist(out.mcmc.Geweke[,"population.mean"]))$p.value
,ks.test(unlist(out.mcmc.Coda.Geweke[,"population.mean"]),unlist(out.mcmc.Geweke[,"population.mean"]))$p.value
,ks.test(unlist(out.mcmc[,"population.mean"]),unlist(out.mcmc.Coda.Geweke[,"population.mean"]))$p.value
),c(ks.test(unlist(out.mcmc[,"population.sd"]),unlist(out.mcmc.Geweke[,"population.sd"]))$p.value
,ks.test(unlist(out.mcmc.Coda.Geweke[,"population.sd"]),unlist(out.mcmc.Geweke[,"population.sd"]))$p.value
,ks.test(unlist(out.mcmc[,"population.sd"]),unlist(out.mcmc.Coda.Geweke[,"population.sd"]))$p.value
)),digits=4,scientific=FALSE)
rownames(compar_MCMC_params)<-c("default vs. Geweke","Coda.Geweke vs. Geweke","Default vs. Coda.Geweke")


## ----Nimble comparison parameters of first three runs invisible, echo=FALSE----

knitr::kable(compar_MCMC_params,col.names=c("mean","sd"),align="r",caption="P-values of paired Kolmogorov-Smirnov tests of output parameters (columns) between the first three NIMBLE models (rows):")


## ----Nimble summary stats, echo=FALSE-----------------------------------------

knitr::kable(format(summary(out.mcmc.Coda.Geweke)[[1]],digits=2,scientific=FALSE),align="r",caption="Summary of the statistical parameters of the Nimble Coda.Geweke model:")

knitr::kable(format(summary(out.mcmc.Geweke)[[1]],digits=2,scientific=FALSE),align="r",caption="Summary of the statistical parameters of the Nimble Geweke model:")

knitr::kable(format(summary(out.mcmc)[[1]],digits=2,scientific=FALSE),align="r",caption="Summary of the statistical parameters of the Nimble default model:")


## ----Jags Data----------------------------------------------------------------

ModelData.Jags <-list(mass = y1000, nobs = length(y1000))



## ----Jags ModelCode-----------------------------------------------------------

modeltotransfer<-"model {

		# Priors
			population.mean ~ dunif(0,5000)
			population.sd ~ dunif(0,100)

			# Normal distribution parameterized by precision = 1/variance in Jags
    	population.variance <- population.sd * population.sd
      precision <- 1 / population.variance

			# Likelihood
			for(i in 1:nobs){
			  mass[i] ~ dnorm(population.mean, precision)
			}
		}"



## ----Jags First model, cache=FALSE--------------------------------------------

set.seed(1)
out.mcmc.Jags<-runMCMC_btadjust(code=modeltotransfer,  data = ModelData.Jags, MCMC_language="Jags", 
    Nchains=Nchains, params=params, inits=Inits,
    niter.min=1000,niter.max=300000,
    nburnin.min=100,nburnin.max=200000,
    thin.min=1,thin.max=1000,
		conv.max=1.05,neff.min=1000)


## ----Jags First model inv, echo=FALSE-----------------------------------------

knitr::kable(format(summary(out.mcmc.Jags)[[1]],digits=2,scientific=FALSE),align="r",caption="Summary of the statistical parameters of the Jags model:")


## ----Nimble-Jags comparison parameters, echo=FALSE, warning=FALSE-------------
compar_MCMC_params<-format(cbind(c(ks.test(unlist(out.mcmc.Geweke[,"population.mean"]),unlist(out.mcmc.Jags[,"population.mean"]))$p.value
,ks.test(unlist(out.mcmc.Coda.Geweke[,"population.mean"]),unlist(out.mcmc.Jags[,"population.mean"]))$p.value
,ks.test(unlist(out.mcmc[,"population.mean"]),unlist(out.mcmc.Jags[,"population.mean"]))$p.value
),c(ks.test(unlist(out.mcmc.Geweke[,"population.sd"]),unlist(out.mcmc.Jags[,"population.sd"]))$p.value
,ks.test(unlist(out.mcmc.Coda.Geweke[,"population.sd"]),unlist(out.mcmc.Jags[,"population.sd"]))$p.value
,ks.test(unlist(out.mcmc[,"population.sd"]),unlist(out.mcmc.Jags[,"population.sd"]))$p.value
)),digits=4,scientific=FALSE)
rownames(compar_MCMC_params)<-c("Nimble.Geweke vs. Jags","Nimble.Coda.Geweke vs. Jags","Nimble.Default vs. Jags")
knitr::kable(compar_MCMC_params,col.names=c("mean","sd"),align="r", caption="P-values of paired Kolmogorov-Smirnov tests of output parameters (columns) of the Jags model with the three NIMBLE models (rows):")


## ----Nimble-Jags comparison efficiency, echo=FALSE----------------------------
compar_MCMC_times<-format(cbind(unlist(attributes(out.mcmc)$final.params),unlist(attributes(out.mcmc.Jags)$final.params)),digits=4,scientific=FALSE)
knitr::kable(compar_MCMC_times[1:(dim(compar_MCMC_times)[1]-1),],col.names=c("Nimble.default","Jags"),align="r", caption="Comparison of the efficiency of the default NIMBLE model and the Jags model:")


## ----Nimble-Jags second comparison efficiency, echo=FALSE---------------------
compar_MCMC_times<-format(cbind(c(attributes(out.mcmc)$final.diags$neff_synth[,"min"
], (attributes(out.mcmc)$final.params$CPUduration.MCMC.transient+attributes(out.mcmc)$final.params$CPUduration.MCMC.asymptotic)/attributes(out.mcmc)$final.diags$neff_synth[,"min"
]),c(attributes(out.mcmc.Jags)$final.diags$neff_synth[,"min"
], (attributes(out.mcmc.Jags)$final.params$CPUduration.MCMC.transient+attributes(out.mcmc.Jags)$final.params$CPUduration.MCMC.asymptotic)/attributes(out.mcmc.Jags)$final.diags$neff_synth[,"min"
])),digits=4,scientific=FALSE)
row.names(compar_MCMC_times)<-c("Min. Number Eff. values","MCMC CPU time per Effective Value")


## ----Nimble-Jags second comparison efficiency bis, echo=FALSE-----------------
knitr::kable(compar_MCMC_times,col.names=c("Nimble.default","Jags"),align="r",caption="Comparison of the number of effective values between the default NIMBLE model and the JAGS model:")


## ----greta model preparation--------------------------------------------------
#in my setting I need to load not only greta but R6 & tensorflow packages
library(greta)
library (R6)
library(tensorflow)

#first requirement of greta: declaring the data that will be analyzed with the function as_data
Y<-as_data(y1000)

#we then proceed by writing the model directly in R, starting with the priors of the parameters using greta functions for probability distributions - here uniform()
population.mean<-uniform(0,5000)
population.sd<-uniform(0,100)
    
#we then define the distribution of the data - here with the normal distribution - by default parametrized with a standard deviation in greta:
try({distribution(Y)<-normal(population.mean,population.sd) })

#we finally declare the greta model, which will be the object passed to runMCMC_btadjust 
m<-model(population.mean, population.sd)

### we finally have to prepare initial values with a specific greta function - initials:
ModelInits.Greta <- function()
    {initials(population.mean = rnorm(1,600,90), population.sd = runif(1, 1, 30))}

set.seed(1)
  Inits.Greta<-lapply(1:Nchains,function(x){ModelInits.Greta()})


## ----greta model running with runMCMC_btadjust, cache=FALSE-------------------
out.mcmc.greta<-runMCMC_btadjust(model=m, MCMC_language="Greta",
    Nchains=Nchains,params=params,inits=Inits.Greta,
		niter.min=1000,niter.max=300000,
    nburnin.min=100,nburnin.max=200000,
		thin.min=1,thin.max=1000,
		conv.max=1.05, neff.min=1000)



## ----greta model running with runMCMC_btadjust table, echo=FALSE--------------

knitr::kable(format(summary(out.mcmc.greta)[[1]],digits=2,scientific=FALSE),align="r",caption="Summary of the statistical parameters of the greta model:")


## ----Nimble-Jags-Greta comparison parameters, echo=FALSE, warning=FALSE-------
compar_MCMC_params<-format(cbind(c(ks.test(unlist(out.mcmc[,"population.mean"]),unlist(out.mcmc.greta[,"population.mean"]))$p.value, ks.test(unlist(out.mcmc.Jags[,"population.mean"]),unlist(out.mcmc.greta[,"population.mean"]))$p.value),c(ks.test(unlist(out.mcmc[,"population.sd"]),unlist(out.mcmc.greta[,"population.sd"]))$p.value
,ks.test(unlist(out.mcmc.Jags[,"population.sd"]),unlist(out.mcmc.greta[,"population.sd"]))$p.value
)),digits=4,scientific=FALSE)
rownames(compar_MCMC_params)<-c("Nimble vs. greta","Jags vs. greta")
knitr::kable(compar_MCMC_params,col.names=c("mean","sd"),align="r",caption="P-values of paired Kolmogorov-Smirnov tests of output parameters of the greta model with the default NIMBLE model and the JAGS model:")


## ----Nimble-Jags-Greta comparison efficiency, echo=FALSE----------------------
compar_MCMC_times<-format(cbind(unlist(attributes(out.mcmc)$final.params),unlist(attributes(out.mcmc.Jags)$final.params),unlist(attributes(out.mcmc.greta)$final.params)),digits=4,scientific=FALSE)
knitr::kable(compar_MCMC_times[1:(dim(compar_MCMC_times)[1]-1),],col.names=c("Nimble.default","Jags","Greta"),align="r",caption="Comparison of the efficiency of the default NIMBLE, the JAGS and the greta models:")


## ----second greta model running with runMCMC_btadjust, cache=FALSE------------
Nchains.Greta<-15
ModelInits.Greta <- function()
    {initials(population.mean = rnorm(1,600,90), population.sd = runif(1, 1, 30))}

set.seed(1)
Inits.Greta<-lapply(1:Nchains.Greta,function(x){ModelInits.Greta()})
  
  out.mcmc.greta.morechains<-runMCMC_btadjust(model=m, MCMC_language="Greta",
    Nchains=Nchains.Greta,params=params,inits=Inits.Greta,
		niter.min=1000,niter.max=300000,
    nburnin.min=100,nburnin.max=200000,
		thin.min=1,thin.max=1000,
		conv.max=1.05, neff.min=1000)
  


## ----second greta model running with runMCMC_btadjust - inv, echo=FALSE-------

knitr::kable(format(summary(out.mcmc.greta.morechains)[[1]],digits=2,scientific=FALSE),align="r",caption=paste0("Summary of the statistical parameters of the greta model with ",Nchains.Greta," chains:"))



## ----Nimble-Jags-second Greta comparison parameters, echo=FALSE, warning=FALSE----
compar_MCMC_params<-format(cbind(c(ks.test(unlist(out.mcmc[,"population.mean"]),unlist(out.mcmc.greta.morechains[,"population.mean"]))$p.value, ks.test(unlist(out.mcmc.Jags[,"population.mean"]),unlist(out.mcmc.greta.morechains[,"population.mean"]))$p.value),c(ks.test(unlist(out.mcmc[,"population.sd"]),unlist(out.mcmc.greta.morechains[,"population.sd"]))$p.value
,ks.test(unlist(out.mcmc.Jags[,"population.sd"]),unlist(out.mcmc.greta.morechains[,"population.sd"]))$p.value
)),digits=4,scientific=FALSE)
rownames(compar_MCMC_params)<-c("Nimble vs. greta.morechains","Jags vs. greta.morechains")
knitr::kable(compar_MCMC_params,col.names=c("mean","sd"),align="r",caption="P-values of paired Kolmogorov-Smirnov tests of output parameters of the greta model with 15 chains with the default NIMBLE model and the JAGS model:")


## ----Nimble-Jags-second Greta comparison efficiency, echo=FALSE---------------
compar_MCMC_times<-format(cbind(unlist(attributes(out.mcmc)$final.params),unlist(attributes(out.mcmc.Jags)$final.params),unlist(attributes(out.mcmc.greta.morechains)$final.params)),digits=4,scientific=FALSE)
knitr::kable(compar_MCMC_times[1:(dim(compar_MCMC_times)[1]-1),],col.names=c("Nimble.default","Jags","Greta.morechains"),align="r",caption="Comparison of the efficiency of the default NIMBLE, the JAGS and the greta.morechains models:")


