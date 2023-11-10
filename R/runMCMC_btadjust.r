
#' @title runMCMC_btadjust
#' @description  returns a mcmc.list object which is the output of a Markov Chain Monte Carlo obtained after adjusting burn-in & thinning parameters to meet pre-specified criteria in terms of convergence & effective sample size - i.e. sample size adjusted for autocorrelation - of the MCMC output
#'
#' @param MCMC_language character value: designates the \code{MCMC_language} used to write & fit the Bayesian model in R. Current choices are "Nimble" - the default-, "Greta" or "Jags". Note that in case it is "Nimble", package \code{nimble} should be loaded in your search list.
#' @param code R object: code for the model that will be used to build the MCMC when \code{MCMC_language} is "Nimble" or "Jags". If "Nimble", must be the name (in R) of the object which is the result of the function \code{nimbleCode}. If "Jags", should be either: (i) a character string which is the name of a txt file that contains the code of the model (as used in the function jags.model): should then end up by ".txt"; or (ii) a character string that contains the text of the Jags code.
#' @param data R list: a list that will contain the data when \code{MCMC_language} is "Nimble" or "Jags". If "Nimble", will be sent to the \code{data} argument of the \code{nimbleModel} function in \code{nimble} package, i.e. the data that have a random distribution in the model. If \code{MCMC_language} is "Greta", can be used just to document the summary of data in the output.
#' @param constants R list: a list that will contain the rest of the data (in addition to data) when \code{MCMC_language} is "Nimble". Will be sent to the \code{constants} argument of the \code{nimbleModel} function in \code{nimble} package, i.e. the data that do not have a random distribution in the model. If \code{MCMC_language} is "Greta", can be used just to document the summary of other data in the output.
#' @param model R object: should be the result of the \code{model} command of Greta when \code{MCMC_language} is "Greta".
#' @param Nchains integer value : the number of Markov chains to run in the MCMC.
#' @param inits an R list, with \code{Nchains} components. Each component is a named list that contains the initial values of the parameters for the MCMC.
#' In case \code{MCMC_language=="Greta"}, each component should be the result of the \code{initials} function in \code{greta} package.
#' @param params character vector: contains the names of the parameters to save at the end of the MCMC and to monitor for convergence and effective sample size;
#' inactive for convergence/effective sample size if \code{params.conv} is specified;
#' inactive for saving if \code{params.save} is specified.
#' @param params.conv character vector: contains the names of the parameters to monitor for convergence and effective sample size.
#' @param params.save character vector: contains the names of the parameters to be saved at the end of the MCMC.
#' @param niter.min integer value: the minimum number of iterations in each chain of the MCMC.
#' @param niter.max integer value: the maximum number of iterations in each chain of the MCMC. Will stop the MCMC once the number of iterations will reach this limit.
#' @param nburnin.min integer value: the minimum number of burn-in (=transitory) iterations in each chain of the MCMC.
#' @param nburnin.max integer value: the maximum number of burn-in (=transitory) iterations in each chain of the MCMC. Will stay at this burn-in value once this limit is reached.
#' @param thin.min integer value: the minimum value of the thin parameter of the MCMC.
#' @param thin.max integer value: the maximum value of the thin parameter of the MCMC. Will stay at this thin value once this limit is reached.
#' @param neff.min positive real number: minimum effective sample size - over parameters used to diagnose convergence & effective sample size- , as calculated with \code{neff.method} (specified in \code{Control}).
#' The algorithm will not stop if the minimum number of efficient values is not above this value (unless another limit - e.g. \code{niter.max} - is reached).
#' @param neff.med positive real number: median effective sample size - over parameters used to diagnose convergence & effective sample size- , as calculated with \code{neff.method} (specified in \code{Control}).
#' The algorithm will not stop if the median number of efficient values is not above this value (unless another limit - e.g. \code{niter.max} - is reached).
#' @param neff.mean positive real number: mean effective sample size - over parameters used to diagnose convergence & effective sample size-, as calculated with \code{neff.method} (specified in \code{Control}).
#' The algorithm will not stop if the mean number of efficient values is not above this value (unless another limit - e.g. \code{niter.max} - is reached).
#' @param conv.max positive real number: maximum - over parameters used to diagnose convergence & effective sample size - convergence diagnostic, as calculated with \code{convtype} method (specified in \code{Control}).
#' The algorithm will not stop if the maximum convergence diagnostic is not below this value (unless another limit - e.g. \code{niter.max} - is reached).
#' @param conv.med positive real number: median - over parameters used to diagnose convergence & effective sample size - convergence diagnostic, as calculated with \code{convtype} method(specified in \code{Control}).
#' The algorithm will not stop if the median convergence diagnostic is not below this value (unless another limit - e.g. \code{niter.max} - is reached).
#' @param conv.mean positive real number: mean - over parameters used to diagnose convergence & effective sample size - convergence diagnostic, as calculated with \code{convtype} method (specified in \code{Control}).
#' The algorithm will not stop if the mean convergence diagnostic is not below this value (unless another limit - e.g. \code{niter.max} - is reached).

#' @param control list of \code{runMCMC_btadjust} control parameters: with the following components:\cr
#'  \itemize{ \item \code{time.max}: positive number (units: seconds): maximum time of the process in seconds; the program will organize itself to stop before \code{0.95*time.max}. Default to NULL, corresponding to no time constraint.
#'   \item \code{check.convergence}: logical value: should the program check convergence at all? Default to TRUE. See Details.
#'   \item \code{check.convergence.firstrun}: logical value: should we check convergence after the first run? Default to NULL in which case will depend on \code{MCMC_language}: if "Greta", will be TRUE because warmup phase separated from the rest; otherwise will be FALSE.
#'   \item \code{recheck.convergence}: logical value: should the algorithm recheck convergence once convergence has been found in a previous run? Default to TRUE.
#'   \item \code{contype}: character value: specifies the type of convergence diagnostic used. Currently implemented: "Gelman" for original Gelman-Rubin diagnostic (only possible if \code{Nchains>=2}), "Gelman_new" for the version of the Gelman-Rubin diagnostic in the second version of "Bayesian Data Analysis" (Gelman, Carlin, Stern and Rubin)(only possible if \code{Nchains>=2}), "Geweke" for Geweke diagnostic (at present applied only in case \code{Nchains==1}) and "Heidleberger" for the reciprocal of Heidleberger-Welch first part of convergence diagnostic based on the Cramer-von Mises test statistic.
#'   \item \code{convtype.Gelman}: integer value: when \code{convtype=="Gelman"}, do we target the Point estimate diagnostic (value 1) or the Upper C.I. diagnostic (value 2). Default to 2.
#'   \item \code{convtype.Geweke}: real vector with two components between 0 and 1: (i) the fraction of samples to consider as the beginning of the chain (frac1 in geweke.diag); (ii) the fraction of samples to consider as the end of the chain (frac2 in \code{gewke.diag}). Default to c(0.1,0.5) as in \code{geweke.diag}.
#'   \item \code{convtype.alpha}: real value between 0 and 1: significance level used in case \code{convtype=="Gelman"} and \code{convtype.Gelman==2}, or \code{convtype=="Heidleberger"}
#'   \item \code{neff.method}: character value: method used to calculate the effective sample sizes. Current choice between "Stan" (the default) and "Coda". If "Stan", uses the function \code{monitor} in package \code{rstan}. If "Coda", uses the function \code{effectiveSize} in package \code{coda}.
#'   \item \code{Ncycles.target}: integer value: targeted number of MCMC runs. Default to 2.
#'   \item \code{props.conv}: numeric vector: in case of non convergence: quantiles of number of iterations removed to recheck convergence. Values should be between 0 and 1.
#'   \item \code{min.Nvalues}: integer value: minimum number of values to diagnose convergence or level of autocorrelation.
#'   \item \code{round.thinmult}: logical value: should the thin multiplier be rounded to the nearest integer so that past values are precisely positioned on the modified iteration sequence? Default to TRUE. Value of FALSE may not be rigorous or may not converge well.
#'   \item \code{min.thinmult}: numeric value: minimum value of thin multiplier: if diagnostics suggest to multiply by less than this, this is not done and the current situation of autocorrelation is considered OK.
#'   \item \code{seed}: integer number: seed for the pseudo-random number generator inside runMCMC_btadjust.
#'   \item \code{identifier.to.print}: character string: printed each time an MCMC update is ran to identify the model (esp. if multiple successive calls to \code{runMCMC_btadjust} are made).
#'   \item \code{safemultiplier.Nvals}: positive number: number bigger than 1 used to multiply the targeted number of efficient values in calculations of additional number of iterations.
#'   \item \code{print.diagnostics}: logical value: should diagnostics be printed each time they are calculated? Default to FALSE.
#'   \item \code{print.thinmult}: logical value: should the raw multiplier of thin be printed each time it is calculated? Default to TRUE.
#'   \item \code{innerprint}: logical value: should printings be done inside the function \code{monitor} of \code{rstan} in case \code{neff.method=="Stan"}? Default to FALSE.
#'   \item \code{remove.fixedchains}: logical value: should we remove Markov chains that do not vary (i.e. whose all parameters have zero variances)? Default to TRUE.
#'   \item \code{check.installation}: logical value: should the function check installation of packages and programs? Default to TRUE.
#' }
#' @param control.MCMC list of MCMC control parameters: with the following components - that depend on \code{MCMC_language}:
#'  \itemize{ \item \code{confModel.expression.toadd} (only for \code{MCMC_language=="Nimble"}): expression to add to \code{confModel} to specify samplers, remove nodes... \code{confModel} should be referred to by \code{confModel[[i]]}. See Details for an example.
#' \item \code{sampler} (only for \code{MCMC_language=="Greta"}): expression used to specify the sampler used.
#' \item \code{warmup} (only for \code{MCMC_language=="Greta"}): integer value used as warmup parameter in the mcmc.
#' \item \code{n.adapt} (only for \code{MCMC_language=="Jags"}): integer value: number of iterations used for adaptation (in function \code{jags.model} in \code{rjags} package).
#' \item \code{RNG.names} (only for \code{MCMC_language=="Jags"}): character vector: name of pseudo-random number generators for each chain. Each component of the vector should be among "base::Wichmann-Hill", "base::Marsaglia-Multicarry", "base::Super-Duper", "base::Mersenne-Twister". If less values than \code{Nchains} are provided, they are specified periodically.
#' \item \code{n_cores} (only for \code{MCMC_language=="Greta"}): integer or NULL: maximum number of cores to use by each sampler.
#' \item \code{showCompilerOutput} (only for \code{MCMC_language=="Nimble"}): logical value indicating whether details of C++ compilation should be printed. Default to TRUE.
#' \item \code{buildDerivs} (only for \code{MCMC_language=="Nimble"}): logical value indicating derivatives should be prepared when preparing Nimble model (will esp. allow to use HMC sampler). Default to FALSE.
#' }
#'
#'
#' @return a \code{mcmc.list} object with attributes with the following components:
#'       \itemize{ \item \code{call.params}: a list containing most of the important arguments of the \code{runMCMC_btadjust} call as well as a summary of dimensions/lengths and mean of components of \code{data} and \code{constants} arguments.
#'       \item \code{final.params}: a list with the parameters of the MCMC at the end of fitting:
#'             \itemize{ \item \code{burnin}: number of iterations of the transient (burn-in) period
#'             \item \code{thin}: number of iterations used for thinning the final output
#'             \item \code{niter.tot}: total number of iterations (of each MCMC chain)
#'             \item \code{duration}: total duration (elapsed time) of the fit (in seconds)
#'             \item \code{duration.MCMC.preparation}: duration (elapsed time) of MCMC preparation (in seconds)
#'             \item \code{duration.MCMC.transient}: duration (elapsed time) of the MCMC transient (burn-in) phase (in seconds)
#'             \item \code{duration.MCMC.asymptotic}: duration (elapsed time) of the MCMC asymptotic phase (in seconds)
#'             \item \code{duration.btadjust}: duration (elapsed time) outside MCMC preparation & fitting (in seconds)
#'             \item \code{CPUduration}: total CPU duration (user+system) of the fit (in seconds)
#'             \item \code{CPUduration.MCMC.preparation}: CPU duration (user+system) of MCMC preparation (in seconds)
#'             \item \code{CPUduration.MCMC.transient}: CPU duration (user+system) of the MCMC transient (burn-in) phase (in seconds)
#'             \item \code{CPUduration.MCMC.asymptotic}: CPU duration (user+system) of the MCMC asymptotic phase (in seconds)
#'             \item \code{CPUduration.btadjust}: CPU duration (user+system) outside MCMC preparation & fitting (in seconds)
#'             \item \code{time}: time (from Sys.time) at the end of model fitting
#'             }
#'        \item \code{final.diags}: a list with final diagnostics of the fit:
#'             \itemize{ \item \code{params}: parameters of the MCMC (burn-in, thin, niter...)
#'             \item \code{conv_synth}: synthetic output of convergence diagnostics
#'             \item \code{neff_synth}: synthetic output for calculations of effective sample sizes
#'             \item \code{conv}: raw convergence values for all the parameters being diagnosed
#'             \item \code{neff}: raw effective sample size values for all the parameters being diagnosed
#'             }
#'        \item \code{package.versions}: a named vector with the versions of the packages used during MCMC fitting
#'        \item \code{R.version}: a list giving the different details of the R setup during MCMC fitting
#'        \item \code{warnings}: a list of the warning messages issued during fitting; unsure it still works with this version
#'        \item \code{error}: a list with the error messages issued during fitting; unsure it still works with this version
#'        }
#'
#'@details
#' Recap: \cr
#' If \code{MCMC_language=="Nimble"}, the code, data and constants arguments should be specified according to the requirements of \code{nimble} package. \cr
#' If \code{MCMC_language=="Jags"}, the code and  data arguments need to be specified as required by \code{rjags} package.  \cr
#' If \code{MCMC_language=="Greta"}, the model argument must be specified and should be the result of the \code{model} command in \code{greta} package.  \cr
#'
#' Details on \code{check.convergence}: \cr
#' If FALSE, no check of convergence at all, after \code{nburnin.min} (& \code{recheck.convergence} is put to FALSE & \code{check.convergence.firstrun} is dominated by \code{check.convergence}). \cr
#' If TRUE, the convergence behavior is governed by \code{check.convergence.firstrun} & \code{recheck.convergence}.\cr
#'
#' Example for \code{confModel.expression.toadd} component of \code{control.MCMC}:\cr
#' 		\code{confModel.expression.toadd<-expression({ConfModel[[i]]$removeSamplers(c("alpha","dzetad","beta","exper_bias[2]","exper_bias[3]","exper_precision[2]","exper_precision[3]"))
#' 			ConfModel[[i]]$addSampler(target = c("alpha","dzetad","beta"),type = "RW_block")
#' 			ConfModel[[i]]$addSampler(target = c("exper_bias[2]","exper_bias[3]"),type = "RW_block")
#' 			ConfModel[[i]]$addSampler(target = c("exper_precision[2]","exper_precision[3]"),type = "RW_block")
#' 			})}
#'
#' Remark for \code{params, params.conv, params.save}:\cr
#' in cases of parameters that are vectors, matrices... the \code{params} vector can contain only the name of the vector or matrix... in which case all its components will be used. It can also contain the names of individual components.
#'
#'
#' @examples
#'  #\code{
#' # for examples with Nimble or Greta, see the Vignette.
#' # condition variable of whether installation is OK with Jags to avoid error durong package check
#' condition_jags<-TRUE
#' if (nchar(system.file(package='rjags'))==0) {condition_jags<-FALSE}
#' if (nchar(system.file(package='runjags'))==0) {condition_jags<-FALSE}
#' if (condition_jags)
#' {suppressWarnings(temp<-runjags::testjags(silent=TRUE))
#'  if(!(temp$JAGS.available&temp$JAGS.found&temp$JAGS.major==4)) {condition_jags<-FALSE}}
#'
#' if (condition_jags) {
#' #generating data
#' set.seed(1)
#' y1000<-rnorm(n=1000,mean=600,sd=30)
#' ModelData <-list(mass = y1000,nobs = length(y1000))
#'
#' #writing the Jags code as a character chain in R
#' modeltotransfer<-"model {
#'
#' # Priors
#' population.mean ~ dunif(0,5000)
#' population.sd ~ dunif(0,100)
#'
#' # Precision = 1/variance: Normal distribution parameterized by precision in Jags
#' population.variance <- population.sd * population.sd
#' precision <- 1 / population.variance
#'
#' # Likelihood
#' for(i in 1:nobs){
#'   mass[i] ~ dnorm(population.mean, precision)
#'  }
#'  }"
#'
#' #specifying the initial values
#' ModelInits <- function()
#' {list (population.mean = rnorm(1,600,90), population.sd = runif(1, 1, 30))}
#' params <- c("population.mean", "population.sd", "population.variance")
#' K<-3
#' set.seed(1)
#' Inits<-lapply(1:K,function(x){ModelInits()})
#'
#' # running runMCMC_btadjust with MCMC_language="Jags":
#' set.seed(1)
#' out.mcmc.Coda<-runMCMC_btadjust(MCMC_language="Jags", code=modeltotransfer,
#' data=ModelData,
#' Nchains=K, params=params, inits=Inits,
#' niter.min=1000, niter.max=300000,
#' nburnin.min=100, nburnin.max=200000,
#' thin.min=1, thin.max=1000,
#' neff.min=1000, conv.max=1.05,
#' control=list(print.diagnostics=TRUE, neff.method="Coda"))
#'
#' summary(out.mcmc.Coda)
#' }
#' #}
#' @importFrom stats median rnorm qnorm quantile update var window
#' @importFrom utils sessionInfo capture.output
#' @export
#'
runMCMC_btadjust<-function(code=NULL,data=NULL,constants=NULL,model=NULL,MCMC_language="Nimble",
						Nchains,inits=NULL,
						params=NULL,params.conv=NULL,params.save=NULL,
						niter.min=100,niter.max=Inf,nburnin.min=10,nburnin.max=Inf,thin.min=1,thin.max=Inf,
						neff.min=NULL,neff.med=NULL,neff.mean=NULL,
						conv.max=NULL,conv.med=NULL,conv.mean=NULL,
						control=list(time.max=NULL,
						check.convergence=TRUE,check.convergence.firstrun=NULL,recheck.convergence=TRUE,
						convtype="Gelman",convtype.Gelman=2,convtype.Geweke=c(0.1,0.5),convtype.alpha=0.05,neff.method="Stan",
						Ncycles.target=2,props.conv=c(0.25,0.5,0.75),min.Nvalues=300,
						min.thinmult=1.1,safemultiplier.Nvals=1.2,round.thinmult=TRUE,
						identifier.to.print="",print.diagnostics=FALSE,print.thinmult=TRUE,innerprint=FALSE,seed=1,remove.fixedchains=TRUE,check.installation=TRUE),
						control.MCMC=list(confModel.expression.toadd=NULL,sampler=expression(hmc()),warmup=1000,n.adapt=1000,RNG.names=c("base::Wichmann-Hill", "base::Marsaglia-Multicarry", "base::Super-Duper","base::Mersenne-Twister"),n_cores=NULL,showCompilerOutput=TRUE,buildDerivs=FALSE))

{

time.MCMC.Preparation0<-Sys.time()
CPUtime.MCMC.Preparation<-0
CPUtime.MCMC<-0
CPUtime.btadjust<-0

#temporarily added waiting for the resolution of mere call to 'nimble::': cf. https://groups.google.com/g/nimble-users/c/-pDY9dB9ovE
#if (MCMC_language=="Nimble") {library(nimble)}


 	current.CPU.time<-system.time({
	#starting counting time
	time.start<-Sys.time()

	########################################
	## 00- code of functions that will be used in the sequel:
	########################################

	list.as.array<-function(x)
	{#transforms a list whose components are all of the same dims (either: vectors, matrices or arrays) into an array whose last component is the level in the list
		array(unlist(x),dim=c(dim(as.array(x[[1]])),length(x)),dimnames=c(dimnames(as.array(x[[1]])),list(names(x))))
	}

	list.as.matrixbisdim1<-function(x)
	{#transforms a list whose components are all matrices with the same number of columns into a matrix with the same number of columns

		## checking all elements should be matrices
		test1<-vapply(x,function(x){length(dim(x))},c(2))
		if (sum(test1!=2) )stop("All elements of the list should be matrices")

		## checking all elements should have the same number of columns
		test2<-vapply(x,function(x){(dim(x))[2]},c(2))
		if (var(test2)!=0) stop("All elements should have the same number of columns")

		  temp<-array(unlist(x),dim=c(dim(as.array(x[[1]])),length(x)))

		 temp<-apply(temp,2,I)
		dimnames(temp)[[2]]<-dimnames(x[[1]])[[2]]

		temp
	}


	formatC_adapted<-function(x,digits=3)
	{#puts numbers in an adequate format for printing
		as.double(formatC(as.double(formatC(x,digits=digits,format="f")),format="g"))
	}

	conveff<-function(out,print.diagnostics=FALSE)
	{### assumes the model is in object "out"
	#out<-samplesList
	#out<-window(samplesList,thin=12)

		if (length(grep("logProb_",dimnames(out[[1]])[[2]]))>0)
		{
			out<-coda::as.mcmc.list(lapply(out,function(x){coda::as.mcmc(x[,-grep("logProb_",dimnames(x)[[2]])])}))
		}
		tempres_params<-cbind.data.frame(Nchains=length(out),thin=thin,niter.tot=niter.tot,Nvalues=dim(out[[1]])[1]*length(out),nu.burn=nburnin.min+sum(numIter.samplesList[1:(index.conv.local-1)]))
		tempres_params_print<-tempres_params
		row.names(tempres_params_print)<-"MCMC parameters"

		if (control$convtype=="Gelman") {
			tempg<-coda::gelman.diag(coda::as.mcmc.list(lapply(out,function(x){coda::as.mcmc(x)})),confidence=1-control$convtype.alpha,autoburnin=FALSE,multivariate=FALSE)[[1]]
			tempres_conv<-cbind.data.frame(max_gelm_Upper_C_I=max(tempg[,2])[1],med_gelm_Upper_C_I=median(tempg[,2])[1],mean_gelm_Upper_C_I=mean(tempg[,2]),max_gelm_Point_Est=max(tempg[,1])[1],name_max_gelm_Point_Est=names(tempg[,1])[tempg[,1]==max(tempg[,1])][1],med_gelm_Point_Est=median(tempg[,1])[1],mean_gelm_Point_Est=mean(tempg[,1]),prop_gelm_Point_Est_above_1p2=mean(tempg[,1]>1.2),prop_gelm_Point_Est_above_1p05=mean(tempg[,1]>1.05),prop_gelm_Point_Est_above_1p01=mean(tempg[,1]>1.01))
			tempres_conv_print<-cbind.data.frame(max=formatC_adapted(c(max(tempg[,2])[1],max(tempg[,1])[1])),median=formatC_adapted(c(median(tempg[,2])[1],median(tempg[,1])[1])),mean=formatC_adapted(c(mean(tempg[,2]),mean(tempg[,1]))),name_max=c(names(tempg[,2])[tempg[,2]==max(tempg[,2])][1],names(tempg[,1])[tempg[,1]==max(tempg[,1])][1]),prop_ab_1p2=formatC_adapted(c(mean(tempg[,2]>1.2),mean(tempg[,1]>1.2))),prop_ab_1p05=formatC_adapted(c(mean(tempg[,2]>1.05),mean(tempg[,1]>1.05))),prop_ab_1p01=formatC_adapted(c(mean(tempg[,2]>1.01),mean(tempg[,1]>1.01))))
			row.names(tempres_conv_print)<-c("Gelman_Upper_C_I","Gelman_Point_Est")
		}

		if (control$convtype=="Gelman_new") {
		  tempg<-as.data.frame(ggmcmc::ggs_Rhat(ggmcmc::ggs(coda::as.mcmc.list(lapply(out,function(x){coda::as.mcmc(x)}))),plot=FALSE))
		  row.names(tempg)<-as.character(tempg[,1])
		  tempg<-tempg[,2,drop=FALSE]
		  tempres_conv<-cbind.data.frame(max_gelm_Point_Est=max(tempg[,1])[1],name_max_gelm_Point_Est=row.names(tempg)[tempg[,1]==max(tempg[,1])][1],med_gelm_Point_Est=median(tempg[,1])[1],mean_gelm_Point_Est=mean(tempg[,1]),prop_gelm_Point_Est_above_1p2=mean(tempg[,1]>1.2),prop_gelm_Point_Est_above_1p05=mean(tempg[,1]>1.05),prop_gelm_Point_Est_above_1p01=mean(tempg[,1]>1.01))
		  tempres_conv_print<-cbind.data.frame(max=formatC_adapted(c(max(tempg[,1])[1])),median=formatC_adapted(c(median(tempg[,1])[1])),mean=formatC_adapted(c(mean(tempg[,1]))),name_max=c(row.names(tempg)[tempg[,1]==max(tempg[,1])][1]),prop_ab_1p2=formatC_adapted(c(mean(tempg[,1]>1.2))),prop_ab_1p05=formatC_adapted(c(mean(tempg[,1]>1.05))),prop_ab_1p01=formatC_adapted(c(mean(tempg[,1]>1.01))))
		  row.names(tempres_conv_print)<-c("Gelman_Point_Est")
		}


		if (control$convtype=="Geweke") {
			tempg<-abs(coda::geweke.diag(coda::as.mcmc.list(lapply(out,function(x){coda::as.mcmc(x)}))[[1]],frac1=control$convtype.Geweke[1],frac2=control$convtype.Geweke[2])$z)
			tempres_conv<-cbind.data.frame(max_gew=max(tempg),med_gew=median(tempg),mean_gew=mean(tempg),name_max_gew=names(tempg)[tempg==max(tempg)][1],prop_gew_above_p975=mean(tempg>qnorm(0.975)),prop_gew_above_p995=mean(tempg>qnorm(0.995)),prop_gew_above_p9995=mean(tempg>qnorm(0.9995)))
			tempres_conv_print<-cbind.data.frame(max=formatC_adapted(max(tempg)),median=formatC_adapted(median(tempg)),mean=formatC_adapted(mean(tempg)),name_max=names(tempg)[tempg==max(tempg)][1],prop_ab_p975=formatC_adapted(mean(tempg>qnorm(0.975))),prop_ab_p995=formatC_adapted(mean(tempg>qnorm(0.995))),prop_ab_p9995=formatC_adapted(mean(tempg>qnorm(0.9995))))
			row.names(tempres_conv_print)<-c("Geweke")
		}

		if (control$convtype=="Heidleberger") {
		  tempg<-1-(coda::heidel.diag(coda::as.mcmc.list(lapply(out,function(x){coda::as.mcmc(x)}))[[1]],pvalue=control$convtype.alpha)[,"stest"])
		  tempres_conv<-cbind.data.frame(max_gew=max(tempg),med_gew=median(tempg),mean_gew=mean(tempg),name_max_gew=names(tempg)[tempg==max(tempg)][1],prop_gew_above_p975=mean(tempg>qnorm(0.975)),prop_gew_above_p995=mean(tempg>qnorm(0.995)),prop_gew_above_p9995=mean(tempg>qnorm(0.9995)))
		  tempres_conv_print<-cbind.data.frame(max=formatC_adapted(max(tempg)),median=formatC_adapted(median(tempg)),mean=formatC_adapted(mean(tempg)),name_max=names(tempg)[tempg==max(tempg)][1],prop_ab_p975=formatC_adapted(mean(tempg>qnorm(0.975))),prop_ab_p995=formatC_adapted(mean(tempg>qnorm(0.995))),prop_ab_p9995=formatC_adapted(mean(tempg>qnorm(0.9995))))
		  row.names(tempres_conv_print)<-c("Heidleberger")
		}

		### calculs squares_Ses_ratios
		tempsu<-summary(out)[[1]]
		squared_SE_ratio<-((tempsu[,4])/(tempsu[,3]))^2
		#print(squared_SE_ratio)
		tempres_SEs<-cbind.data.frame(max_sqSEratio=max(squared_SE_ratio)[1],name_max_sqSEratio=names(squared_SE_ratio)[squared_SE_ratio==max(squared_SE_ratio)][1],med_sqSEratio=median(squared_SE_ratio)[1],prop_sqSEratio_above_2=mean(squared_SE_ratio>2),prop_sqSEratio_above_1p2=mean(squared_SE_ratio>1.2))

		if (control$neff.method=="Stan")
		{
			MonModelStan<- aperm(list.as.array(out),c(1,3,2))
			dimnames(MonModelStan)<-list(NULL,NULL,dimnames(out[[1]])[[2]])
			tempm<-rstan::monitor(MonModelStan,print=control$innerprint,warmup=0)[,"n_eff"]
		}

		if (control$neff.method=="Coda")
		{
			tempm<-coda::effectiveSize(out)
		}


		tempres_neff<-cbind.data.frame(min_neff=min(tempm)[1],name_min_neff=names(tempm)[tempm==min(tempm)][1],med_neff=median(tempm),mean_neff=mean(tempm),prop_neff_below_1000=mean(tempm<1000),prop_neff_below_5000=mean(tempm<5000),prop_neff_below_10000=mean(tempm<10000))
		tempres_neff_print<-cbind.data.frame(min=formatC_adapted(min(tempm)[1]),median=formatC_adapted(median(tempm)),mean=formatC_adapted(mean(tempm)),name_min=names(tempm)[tempm==min(tempm)][1],prop_bel_1000=formatC_adapted(mean(tempm<1000)),prop_bel_5000=formatC_adapted(mean(tempm<5000)),prop_bel_10000=formatC_adapted(mean(tempm<10000)))
		row.names(tempres_neff_print)<-"Neff            "
		#tempressimple<-tempres[,c("maxgelm_Point_Est","medelm1","min_neff","med_neff")]
		if (print.diagnostics) {
			print("###################################################################################")
			print("Current state of diagnostics:"); print(tempres_params_print)
			print("###################################################################################")
			print(tempres_conv_print)
			print("###################################################################################")
			print(tempres_neff_print)
			print("###################################################################################")
		}
		tempres<-cbind.data.frame(tempres_params,tempres_conv,tempres_neff,tempres_SEs)

		tempres
	}	#END conveff function


	conveff_final<-function(out,print.diagnostics=FALSE)
	{### assumes the model is in object "out"
	#out<-samplesList
	#out<-window(samplesList,thin=12)

		if (length(grep("logProb_",dimnames(out[[1]])[[2]]))>0)
		{
			out<-coda::as.mcmc.list(lapply(out,function(x){coda::as.mcmc(x[,-grep("logProb_",dimnames(x)[[2]])])}))
		}
		tempres_params<-cbind.data.frame(Nchains=length(out),thin=thin,niter.tot=niter.tot,Nvalues=dim(out[[1]])[1]*length(out),nu.burn=nburnin.min+sum(numIter.samplesList[1:(index.conv.local-1)]))
		tempres_params_print<-tempres_params
		row.names(tempres_params_print)<-"MCMC parameters"

		if (control$convtype=="Gelman") {
			tempg<-coda::gelman.diag(coda::as.mcmc.list(lapply(out,function(x){coda::as.mcmc(x)})),confidence=1-control$convtype.alpha,autoburnin=FALSE,multivariate=FALSE)[[1]]
			tempres_conv<-cbind.data.frame(max_gelm_Upper_C_I=max(tempg[,2])[1],med_gelm_Upper_C_I=median(tempg[,2])[1],mean_gelm_Upper_C_I=mean(tempg[,2]),max_gelm_Point_Est=max(tempg[,1])[1],name_max_gelm_Point_Est=names(tempg[,1])[tempg[,1]==max(tempg[,1])][1],med_gelm_Point_Est=median(tempg[,1])[1],mean_gelm_Point_Est=mean(tempg[,1]),prop_gelm_Point_Est_above_1p2=mean(tempg[,1]>1.2),prop_gelm_Point_Est_above_1p05=mean(tempg[,1]>1.05),prop_gelm_Point_Est_above_1p01=mean(tempg[,1]>1.01))
			tempres_conv_print<-cbind.data.frame(max=formatC_adapted(c(max(tempg[,2])[1],max(tempg[,1])[1])),median=formatC_adapted(c(median(tempg[,2])[1],median(tempg[,1])[1])),mean=formatC_adapted(c(mean(tempg[,2]),mean(tempg[,1]))),name_max=c(names(tempg[,2])[tempg[,2]==max(tempg[,2])][1],names(tempg[,1])[tempg[,1]==max(tempg[,1])][1]),prop_ab_1p2=formatC_adapted(c(mean(tempg[,2]>1.2),mean(tempg[,1]>1.2))),prop_ab_1p05=formatC_adapted(c(mean(tempg[,2]>1.05),mean(tempg[,1]>1.05))),prop_ab_1p01=formatC_adapted(c(mean(tempg[,2]>1.01),mean(tempg[,1]>1.01))))
			row.names(tempres_conv_print)<-c("Gelman_Upper_C_I","Gelman_Point_Est")
		}

		if (control$convtype=="Gelman_new") {
		  tempg<-as.data.frame(ggmcmc::ggs_Rhat(ggmcmc::ggs(coda::as.mcmc.list(lapply(out,function(x){coda::as.mcmc(x)}))),plot=FALSE))
		  row.names(tempg)<-as.character(tempg[,1])
		  tempg<-tempg[,2,drop=FALSE]
		  tempres_conv<-cbind.data.frame(max_gelm_Point_Est=max(tempg[,1])[1],name_max_gelm_Point_Est=row.names(tempg)[tempg[,1]==max(tempg[,1])][1],med_gelm_Point_Est=median(tempg[,1])[1],mean_gelm_Point_Est=mean(tempg[,1]),prop_gelm_Point_Est_above_1p2=mean(tempg[,1]>1.2),prop_gelm_Point_Est_above_1p05=mean(tempg[,1]>1.05),prop_gelm_Point_Est_above_1p01=mean(tempg[,1]>1.01))
		  tempres_conv_print<-cbind.data.frame(max=formatC_adapted(c(max(tempg[,1])[1])),median=formatC_adapted(c(median(tempg[,1])[1])),mean=formatC_adapted(c(mean(tempg[,1]))),name_max=c(row.names(tempg)[tempg[,1]==max(tempg[,1])][1]),prop_ab_1p2=formatC_adapted(c(mean(tempg[,1]>1.2))),prop_ab_1p05=formatC_adapted(c(mean(tempg[,1]>1.05))),prop_ab_1p01=formatC_adapted(c(mean(tempg[,1]>1.01))))
		  row.names(tempres_conv_print)<-c("Gelman_Point_Est")
		}

		if (control$convtype=="Geweke") {
			tempg<-abs(coda::geweke.diag(coda::as.mcmc.list(lapply(out,function(x){coda::as.mcmc(x)}))[[1]],frac1=control$convtype.Geweke[1],frac2=control$convtype.Geweke[2])$z)
			tempres_conv<-cbind.data.frame(max_gew=max(tempg),med_gew=median(tempg),mean_gew=mean(tempg),name_max_gew=names(tempg)[tempg==max(tempg)][1],prop_gew_above_p975=mean(tempg>qnorm(0.975)),prop_gew_above_p995=mean(tempg>qnorm(0.995)),prop_gew_above_p9995=mean(tempg>qnorm(0.9995)))
			tempres_conv_print<-cbind.data.frame(max=formatC_adapted(max(tempg)),median=formatC_adapted(median(tempg)),mean=formatC_adapted(mean(tempg)),name_max=names(tempg)[tempg==max(tempg)][1],prop_ab_p975=formatC_adapted(mean(tempg>qnorm(0.975))),prop_ab_p995=formatC_adapted(mean(tempg>qnorm(0.995))),prop_ab_p9995=formatC_adapted(mean(tempg>qnorm(0.9995))))
			row.names(tempres_conv_print)<-c("Geweke")
		}

		if (control$convtype=="Heidleberger") {
		  tempg<-1-(coda::heidel.diag(coda::as.mcmc.list(lapply(out,function(x){coda::as.mcmc(x)}))[[1]],pvalue=control$convtype.alpha)[,"stest"])
		  tempres_conv<-cbind.data.frame(max_gew=max(tempg),med_gew=median(tempg),mean_gew=mean(tempg),name_max_gew=names(tempg)[tempg==max(tempg)][1],prop_gew_above_p975=mean(tempg>qnorm(0.975)),prop_gew_above_p995=mean(tempg>qnorm(0.995)),prop_gew_above_p9995=mean(tempg>qnorm(0.9995)))
		  tempres_conv_print<-cbind.data.frame(max=formatC_adapted(max(tempg)),median=formatC_adapted(median(tempg)),mean=formatC_adapted(mean(tempg)),name_max=names(tempg)[tempg==max(tempg)][1],prop_ab_p975=formatC_adapted(mean(tempg>qnorm(0.975))),prop_ab_p995=formatC_adapted(mean(tempg>qnorm(0.995))),prop_ab_p9995=formatC_adapted(mean(tempg>qnorm(0.9995))))
		  row.names(tempres_conv_print)<-c("Heidleberger")
		}

		### calculs squares_Ses_ratios
		tempsu<-summary(out)[[1]]
		squared_SE_ratio<-((tempsu[,4])/(tempsu[,3]))^2
		#print(squared_SE_ratio)
		tempres_SEs<-cbind.data.frame(max_sqSEratio=max(squared_SE_ratio)[1],name_max_sqSEratio=names(squared_SE_ratio)[squared_SE_ratio==max(squared_SE_ratio)][1],med_sqSEratio=median(squared_SE_ratio)[1],prop_sqSEratio_above_2=mean(squared_SE_ratio>2),prop_sqSEratio_above_1p2=mean(squared_SE_ratio>1.2))

		if (control$neff.method=="Stan")
		{
			MonModelStan<- aperm(list.as.array(out),c(1,3,2))
			dimnames(MonModelStan)<-list(NULL,NULL,dimnames(out[[1]])[[2]])
			tempm<-rstan::monitor(MonModelStan,print=control$innerprint,warmup=0)[,"n_eff"]
		}

		if (control$neff.method=="Coda")
		{
			tempm<-coda::effectiveSize(out)
		}


		tempres_neff<-cbind.data.frame(min_neff=min(tempm)[1],name_min_neff=names(tempm)[tempm==min(tempm)][1],med_neff=median(tempm),mean_neff=mean(tempm),prop_neff_below_1000=mean(tempm<1000),prop_neff_below_5000=mean(tempm<5000),prop_neff_below_10000=mean(tempm<10000))
		tempres_neff_print<-cbind.data.frame(min=formatC_adapted(min(tempm)[1]),median=formatC_adapted(median(tempm)),mean=formatC_adapted(mean(tempm)),name_min=names(tempm)[tempm==min(tempm)][1],prop_bel_1000=formatC_adapted(mean(tempm<1000)),prop_bel_5000=formatC_adapted(mean(tempm<5000)),prop_bel_10000=formatC_adapted(mean(tempm<10000)))
		row.names(tempres_neff_print)<-"Neff            "
		#tempressimple<-tempres[,c("maxgelm_Point_Est","medelm1","min_neff","med_neff")]
		if (print.diagnostics) {
			print("###################################################################################")
			print("Current state of diagnostics:"); print(tempres_params_print)
			print("###################################################################################")
			print(tempres_conv_print)
			print("###################################################################################")
			print(tempres_neff_print)
			print("###################################################################################")
		}
		#removed because too much unstructured: ,global_synth=cbind.data.frame(tempres_params,tempres_conv,tempres_neff,tempres_SEs)
		tempres<-list(params=tempres_params_print,conv_synth=tempres_conv_print,neff_synth=tempres_neff_print,conv=tempg,neff=tempm)

		tempres
	}	#END conveff_final function


	checking.convergence<-function(diags,convtype,convtype.Gelman,conv.max,conv.med,conv.mean)
	{
		converged <-TRUE
		#print(conv.max)
		#print(diags$max_gelm_Upper_C_I)
		if (convtype=="Gelman")
		{
			if (convtype.Gelman==1)
			{
				if (!is.null(conv.max)) {if (conv.max<diags$max_gelm_Point_Est) {converged<-FALSE}}
				if (!is.null(conv.mean)) {if (conv.mean<diags$mean_gelm_Point_Est) {converged<-FALSE}}
				if (!is.null(conv.med)) {if (conv.med<diags$med_gelm_Point_Est) {converged<-FALSE}}
			}
			if (convtype.Gelman==2)
			{
				if (!is.null(conv.max)) {if (conv.max<diags$max_gelm_Upper_C_I) {converged<-FALSE}}
				if (!is.null(conv.mean)) {if (conv.mean<diags$mean_gelm_Upper_C_I) {converged<-FALSE}}
				if (!is.null(conv.med)) {if (conv.med<diags$med_gelm_Upper_C_I) {converged<-FALSE}}
			}
		}

		if (convtype=="Gelman_new")
		{

		    if (!is.null(conv.max)) {if (conv.max<diags$max_gelm_Point_Est) {converged<-FALSE}}
		    if (!is.null(conv.mean)) {if (conv.mean<diags$mean_gelm_Point_Est) {converged<-FALSE}}
		    if (!is.null(conv.med)) {if (conv.med<diags$med_gelm_Point_Est) {converged<-FALSE}}
		}

		if (convtype=="Geweke")
		{
			if (!is.null(conv.max)) {if (conv.max<diags$max_gew) {converged<-FALSE}}
			if (!is.null(conv.mean)) {if (conv.mean<diags$mean_gew) {converged<-FALSE}}
			if (!is.null(conv.med)) {if (conv.med<diags$med_gew) {converged<-FALSE}}
		}

		if (convtype=="Heidleberger")
		{
		  if (!is.null(conv.max)) {if (conv.max<diags$max_gew) {converged<-FALSE}}
		  if (!is.null(conv.mean)) {if (conv.mean<diags$mean_gew) {converged<-FALSE}}
		  if (!is.null(conv.med)) {if (conv.med<diags$med_gew) {converged<-FALSE}}
		}
		converged
	}	#END checking.convergence function

	## function to check if criteria of neff are reached
	checking.neffs.reached<-function(diags,neff.min,neff.med,neff.mean)
	{
			neffs.reached<-TRUE
			if (!is.null(neff.min)) {if(diags$min_neff<neff.min) {neffs.reached<-FALSE} }
			if (!is.null(neff.med)) {if(diags$med_neff<neff.med) {neffs.reached<-FALSE}}
			if (!is.null(neff.mean)) {if(diags$mean_neff<neff.mean) {neffs.reached<-FALSE}}

		neffs.reached
	}

	## function to calculate the multiplier of thin to do in the next step to render replicates nearly independent
	calculate.thinmult.target<-function(diags,neff.min,neff.med,neff.mean)
	{
		thinmult<-1
		#neffmax<-max(c(neff.min,neff.med,neff.mean))
		N<-diags$Nvalues
		if (!is.null(neff.min)) {thinmult<-max(thinmult,N/diags$min_neff*neff.min/neff.max) }
		if (!is.null(neff.med)) {thinmult<-max(thinmult,N/diags$med_neff*neff.med/neff.max)}
		if (!is.null(neff.mean)) {thinmult<-max(thinmult,N/diags$mean_neff*neff.mean/neff.max)}


		thinmult
	}

	##function to put all the neffs on the same scale and then taking the minimum
	scale.available.neffs<-function(diags,neff.min,neff.med,neff.mean)
	{
		#neffmax<-max(c(neff.min,neff.med,neff.mean))
		neffs<-min(c(diags$min_neff/neff.min*neff.max,diags$med_neff/neff.med*neff.max,diags$mean_neff/neff.mean*neff.max))
		neffs
	}

	## function to "summarize" data (or if Null: ModelData) and constants (or if Null, ModelCosnts)
	summarize.data<-function(data)
	{
	names.data<-names(data)
	res<-as.list(rep(1,2*length(names.data)))
	names(res)<-paste0(rep(c("dims.","mean."),times=length(names.data)),rep(names.data,each=2))
	for (i in seq_len(length(names.data)))
		{if (length(dim(data[[i]]))>0)
			{
			res[[2*(i-1)+1]]<-dim(data[[i]])
			}
		else {
			res[[2*(i-1)+1]]<-length(data[[i]])

			}
		res[[2*(i-1)+2]]<-mean(data[[i]],na.rm=TRUE)
		}
	res
	}
	########################################
	###0-1. Initialisation & checkings
	########################################

	neff.max<-max(neff.min,neff.med,neff.mean)

	### variable that will contain the info on previous convergences
	previously.converged<-FALSE

	### putting the control argument in the good format in case it is specified partially
	control0<-list(time.max=NULL,check.convergence=TRUE,check.convergence.firstrun=NULL,recheck.convergence=TRUE,convtype="Gelman",
							convtype.Gelman=2,convtype.Geweke=c(0.1,0.5),convtype.alpha=0.05,neff.method="Stan",Ncycles.target=2,props.conv=c(0.25,0.5,0.75),min.Nvalues=300,
							min.thinmult=1.1,safemultiplier.Nvals=1.2,round.thinmult=TRUE,
							identifier.to.print="",print.diagnostics=FALSE,print.thinmult=TRUE,innerprint=FALSE,seed=1,remove.fixedchains=TRUE,check.installation=TRUE)

	### checking all the names of arguments of control are in control0; otherwise stops because arguments are not intrepretable
	if (length(setdiff(names(control),names(control0)))>0)
		{stop("The names of the control argument do not match the default names")}

	### puts the components of control into control0:
	control0[names(control)]<-control

	### INACTIVATED: putting in control0 parameters that could have been passed directly to the function, not in the control argument:
	#if (length(intersect(names(c(as.list(environment()), list(...))),names(control0)))>0) {control0[intersect(names(c(as.list(environment()), list(...))),names(control0))]<-c(as.list(environment()), list(...))[intersect(names(c(as.list(environment()), list(...))),names(control0))]}

	### transferring control0 to control
	control<-control0


	### putting the control.MCMC argument in the good format in case it is specified partially: same sequence as for control/control0:
	control.MCMC0<-list(confModel.expression.toadd=NULL,sampler=expression(hmc()),warmup=1000,n.adapt=1000,RNG.names=c("base::Wichmann-Hill", "base::Marsaglia-Multicarry", "base::Super-Duper","base::Mersenne-Twister"),n_cores=NULL,showCompilerOutput=TRUE,buildDerivs=FALSE)
	if (length(setdiff(names(control.MCMC),names(control.MCMC0)))>0)
		{stop("The names of the control.MCMC argument do not match the default names")}
	if (length(intersect(names(control.MCMC),names(control.MCMC0)))<length(names(control.MCMC0)))
		{control.MCMC0[names(control.MCMC)]<-control.MCMC}
	### INACTIVATED: putting in control.MCMC parameters that could have been passed directly to the function, not in the control.MCMC argument:
#	if (length(intersect(names(c(as.list(environment()), list(...))),names(control.MCMC0)))>0) {control.MCMC0[intersect(names(c(as.list(environment()), list(...))),names(control.MCMC0))]<-c(as.list(environment()), list(...))[intersect(names(c(as.list(environment()), list(...))),names(control.MCMC0))]}
	control.MCMC<-control.MCMC0

	## checking the installation of adequate packages and programs
if (control$check.installation)
  {
		if (control$neff.method=="Stan" & nchar(system.file(package='rstan'))==0) {stop("Since parameter neff.method in control argument is \"Stan\", package \"rtsan\" should be installed, which is not the case")}
	  if (control$convtype=="Gelman_new" & nchar(system.file(package='ggmcmc'))==0) {stop("Since parameter convtype in control argument is \"Stan\", package \"ggmcmc\" should be installed, which is not the case")}
	  if (MCMC_language=="Greta"& nchar(system.file(package='greta'))==0) {stop("Since MCMC_language is \"Greta\", package \"greta\" should be installed, which is not the case")}
	  if (MCMC_language=="Greta"& nchar(system.file(package='R6'))==0) {stop("Since MCMC_language is \"Greta\", package \"R6\" should be installed, which is not the case")}
	  if (MCMC_language=="Greta"& nchar(system.file(package='tensorflow'))==0) {stop("Since MCMC_language is \"Greta\", package \"tensorflow\" should be installed, which is not the case")}

  # this control has been temporarily removed since it takes some time and seems not to work within RMarkdown => dangerous
  # much as in: https://github.com/rstudio/rmarkdown/issues/1150
  # request made on greta forum: https://forum.greta-stats.org/t/request-for-a-modification-of-function-greta-sitrep/345
    #if (MCMC_language=="Greta") {test<-capture.output(greta::greta_sitrep(),type="message")
	   #                           if (max(max(regexpr("greta is ready to use!",test))<=0))
	   #                             {stop("Since MCMC_language is \"Greta\", greta should be ready to use, which is not the case - call \"greta_sitrep()\" for more details")}
	    #                          }
	  if (MCMC_language=="Jags"& nchar(system.file(package='rjags'))==0) {stop("Since MCMC_language is \"Jags\", package \"rjags\" should be installed, which is not the case")}
	  if (MCMC_language=="Jags"& nchar(system.file(package='runjags'))==0) {stop("Since MCMC_language is \"Jags\", package \"runjags\" should be installed, which is not the case")}
	  if (MCMC_language=="Jags") {suppressWarnings(temp<-runjags::testjags(silent=TRUE))
	                              if(!(temp$JAGS.available&temp$JAGS.found&temp$JAGS.major==4)) {stop("Since MCMC_language is \"Jags\", program \"JAGS\" with version greater than 4.x.y should be installed, which is not the case")}
	                              }

	  if (MCMC_language=="Nimble"& nchar(system.file(package='nimble'))==0) {stop("Since MCMC_language is \"Nimble\", package \"nimble\" should be installed, which is not the case")}
}

	## modifications of parameters & checks that were brought through control or control.MCMC:
		## if control$check.convergence.firstrun is NULL: will depend on MCMC_language: if "Greta", will be TRUE because warmup phase separated from the rest; same for "Jags" in which the adaptation phase is also separated from the rest; otherwise will be FALSE.
		if (is.null(control$check.convergence.firstrun)) { if (MCMC_language=="Nimble") {control$check.convergence.firstrun<-FALSE} else {control$check.convergence.firstrun<-TRUE}}

		##if control$check.convergence is FALSE, this should be the same for control$recheck.convergence
		if (!control$check.convergence) {control$recheck.convergence<-FALSE}

		## checking that control$props.conv are >=0 (if ==0 will double the first diagnostic) and <=1:
		if (min(control$props.conv)<0) {stop("props.conv values in control should all be nonnegative")}
		if (max(control$props.conv)>1) {stop("props.conv values in control should all be less than 1")}
		## if control$props.conv does not contain 1: add 1 because the algorithm requires it to stop correctly
		if (max(control$props.conv)<1) {control$props.conv<-c(control$props.conv,1)}

		if (length(control$safemultiplier.Nvals)>1)
		{
			print("safemultiplier.Nvals in control of length >1 although should be of length 1; Only first value kept")
			control$safemultiplier.Nvals<-control$safemultiplier.Nvals[1]
		}

		if (!is.numeric(control$safemultiplier.Nvals)) {stop("safemultiplier.Nvals in control is not numerical although it should be")}
		if (control$safemultiplier.Nvals<1) {stop("safemultiplier.Nvals in control is less than 1.0 although it should not be")}

		if (is.null(control$time.max) & !is.finite(niter.max)) {stop("time.max in control is NULL & niter.max is not finite; not possible: can cause MCMC not to stop")}



		## if control$identifier.to.print is specified add ", " afetr it because will be coerced with other characters
		if (control$identifier.to.print!="") { control$identifier.to.print<-paste(control$identifier.to.print,", ",sep="")}

		## checking adequacy of control$convtype
		if (!is.element(control$convtype, c("Gelman","Gelman_new","Geweke","Heidleberger"))) {stop("Parameter convtype in control argument should be equal either to \"Gelman\", \"Gelman_new\", \"Geweke\" or \"Heidleberger\": please change this parameter")}


		## checking adequacy between convtype and number of MCMC chains
		if ((control$convtype=="Gelman"|control$convtype=="Gelman_new")&Nchains==1) {stop("Impossible to use the Gelman-Rubin convergence diagnostic with only one MCMC chain: please change Nchain or control$convtype")}
		if ((control$convtype=="Geweke")&Nchains>1) {stop("Impossible to use the Geweke convergence diagnostic with more than one MCMC chain: please change Nchain or control$convtype")}
	  if ((control$convtype=="Heidleberger")&Nchains>1) {stop("Impossible to use the Heidleberger convergence diagnostic with more than one MCMC chain: please change Nchain or control$convtype")}

	  ## checking adequacy of convtype.alpha
	  if (length(control$convtype.alpha)!=1) {stop("control$convtype.alpha should be of length 1")}
	  if ((control$convtype.alpha)<=0|(control$convtype.alpha)>=1) {stop("control$convtype.alpha should be between 0 and 1")}


		## checking adequacy of control$neff.method
		if (control$neff.method!="Stan"&control$neff.method!="Coda") {stop("Parameter neff.method in control argument should be equal either to \"Stan\" or \"Coda\": please change this parameter")}


		## checking adequacy of MCMC_language
		if (!is.element(MCMC_language,c("Nimble","Jags","Greta"))) {stop("MCMC_language should be either \"Nimble\", \"Jags\" or \"Greta\"; please respecify this parameter accordingly")}

		## checking adequacy of arguments in case MCMC_language=="Nimble"
		if (MCMC_language=="Nimble"&(is.null(code)|is.null(data)|is.null(constants)))
				{stop("Either code or data or constants argument is null although it is required by Nimble to define the model and fit MCMC; please provide this missing parameter")}

  	if (MCMC_language=="Nimble"&!is.element("package:nimble",search()))
    	{stop("It is required that you have \"package:nimble\" in your search list for runMCMC_btadajust to run nicely. Consider doing library(nimble) or require(nimble). You should already have loaded it to build the argument code")}

		## checking adequacy of arguments in case MCMC_language=="Greta"
		if (MCMC_language=="Greta"&(is.null(model)))
				{stop("Model argument is null although it is required by Greta to define the model and fit MCMC; please provide this missing parameter")}
		if (MCMC_language=="Greta"&((!is.null(inits))&!(is.element("initials",c(class(inits),class(inits[[1]]))))))
				{stop("Initial values should be buildt with the greta::initials function in case Greta is the language")}

		## checking adequacy of arguments in case MCMC_language=="Jags"
		if (MCMC_language=="Jags"&(is.null(code)|is.null(data)))
				{stop("Either code or data argument is null although it is required by Jags to define the model and fit MCMC; please provide this missing parameter")}

  	if (MCMC_language=="Jags"&(!is.character(code)))
  	{stop("Code should be a character containing either the link to the text file containing the code or the code itself")}

	  if (is.finite(niter.max))
	  {if ((niter.max/thin.max)>(100*neff.max)) {warning("There is a high risk of very oversized Rdata files since (niter.max/thin.max)>(100*neff.max); consider increasing thin.max")}
	    else if ((niter.max/thin.max)>(10*neff.max)) {warning("There is a high risk of oversized Rdata files since (niter.max/thin.max)>(10*neff.max); consider increasing thin.max")}
	    else if ((niter.max/thin.max)>(3*neff.max)) {warning("There is a risk of slightly oversized Rdata files since (niter.max/thin.max)>(3*neff.max); consider increasing thin.max")}
	  }

	## initial values of variables used on the algorithm

	set.seed(control$seed)

	## initialisation of samplesList
		## samplesList will be the list that will store the MCMC outputs; one the central objects of this package
	samplesList <- vector("list", Nchains)
	names(samplesList) <- paste0("chain", 1:Nchains)

	numIter.samplesList<-NULL ##will contain the number of iterations between successive values in samplesList; same dimension as the number of rows in samplesList[[1]]

	index.conv<-1 ## will contain the index (in rows) of the transient period in number of rows in samplesList[[1]] (as diagnosed by convergence diagnostics)
	thin<-thin.min ## will contain the active thin value
	nburnin<-nburnin.min
	niter<-max(niter.min,nburnin+ceiling(control$min.Nvalues/Nchains)*thin+10)
	thinmult<-1
	Ncycles<-1
	niter.tot<-0
	neffs.reached<-FALSE
	params<-union(params,union(params.conv,params.save))
	if (is.null(params)) stop("Program stopped: at least one of the following arguments should be specified: params, params.conv, params.save.")
	if (is.null(params.conv)) {params.conv<-params}
	if (is.null(params.save)) {params.save<-params}

	########################################
	###0-2. initial preparation of MCMC models: process differs depending on MCMC_language
	########################################
		###NB: nothing to be done in case MCMC_language=="Greta" since all is already prepared and given in argument model


	if (MCMC_language=="Nimble") {
		Model <- vector("list", Nchains)
		CModelMCMC <- vector("list", Nchains)
		CModel <- vector("list", Nchains)
		ConfModel <- vector("list", Nchains)
		ModelMCMC <- vector("list", Nchains)

		for (i in 1:Nchains)
		{Model[[i]] <- nimble::nimbleModel(code = code, name = 'Nimble', constants = constants, data = data, buildDerivs=control.MCMC$buildDerivs)
		CModel[[i]] <- nimble::compileNimble(Model[[i]], showCompilerOutput = control.MCMC$showCompilerOutput)
		ConfModel[[i]] <- nimble::configureMCMC(Model[[i]], thin = thin, monitors = params)
		if (!is.null(control.MCMC$confModel.expression.toadd))
		{eval(control.MCMC$confModel.expression.toadd)}
		ModelMCMC[[i]] <- nimble::buildMCMC(ConfModel[[i]])
		CModelMCMC[[i]] <- nimble::compileNimble(ModelMCMC[[i]], project = CModel[[i]], showCompilerOutput = control.MCMC$showCompilerOutput)
		}
	}
	## END: MCMC_language=="Nimble"

	if (MCMC_language=="Jags") {
		code0<-code
		if (is.character(code)) {
			if (substring(code,first=nchar(code)-3,last=nchar(code))!=".txt")
			{
			 code<-textConnection(code)
			}
		}
    if (is.null(inits))
    {inits<-lapply(1:Nchains,function(x){list(".RNG.name" = control.MCMC$RNG.names[(x-1)%%Nchains+1],".RNG.seed" = control$seed+x)})
    } else
    {
      inits<-lapply(1:Nchains,function(x){c(inits[[x]],list(".RNG.name" = control.MCMC$RNG.names[(x-1)%%Nchains+1],".RNG.seed" = control$seed+x))})
    }
		myModel<-rjags::jags.model(code, data=data, n.chains = Nchains, n.adapt = control.MCMC$n.adapt, inits=inits)

	}
	## END: MCMC_language=="Jags"
	})
	## END: current.CPU.time<-system.time({
	time.MCMC.Preparation<-(Sys.time()-time.MCMC.Preparation0)
	units(time.MCMC.Preparation)<-"secs"
	CPUtime.MCMC.Preparation<-CPUtime.MCMC.Preparation+current.CPU.time[1]+current.CPU.time[2]

	########################################
	###1.1- First run of MCMCs:
	########################################
	current.CPU.time<-system.time({
	time.MCMC0<-Sys.time()
	samplesList <- vector("list", Nchains)
	 names(samplesList) <- paste0("chain", 1:Nchains)

	 message(control$identifier.to.print,"Cycle ", Ncycles, "...")
	 print("###################################################################################")

	 if (MCMC_language=="Nimble") {
	   nburnin.min0<-nburnin.min
	   for (i in 1:Nchains)
			{message("      Running chain ", i, "...")
			 set.seed(control$seed+i)
			 Modeltemp <- {if (nimble::is.Cnf(CModelMCMC[[i]]))
			   CModelMCMC[[i]]$Robject$model$CobjectInterface
				else CModelMCMC[[i]]$model}
			 Modeltemp$setInits(inits[[i]])
			 CModelMCMC[[i]]$run(niter, nburnin = nburnin, thin = thin, progressBar =TRUE)
			 samplesList[[i]] <- as.matrix(CModelMCMC[[i]]$mvSamples)
			}
		samplesList <- coda::as.mcmc.list(lapply(samplesList, coda::as.mcmc))
	 }	 ## END: MCMC_language=="Nimble"

	if (MCMC_language=="Jags") {

	  #running unstored nburnin.min iterations:
	  update(myModel, n.iter=nburnin.min)

	  #changing nburnin.min in case MCMC_language=="Jags", after update, for consistency of calculations of number of iterations:
	  nburnin.min0<-nburnin.min
	  nburnin.min<-0

	  ## as of 2022/12/14: added n.adapt to niter because otherwise output is void. probably because in the adaptive phase. Removed on 2023/04/27 because reversed behaviour: too much values.
	  ## as of 2023/04/20: nburnin.min0 is removed from niter for cross-compatibility with Nimble's behaviour
		resultat<-rjags::coda.samples(myModel,var=params, n.iter = niter-nburnin.min0, thin = thin)
		matres<-list.as.matrixbisdim1(resultat)
		#format changed due to particular value for start with Jags: https://stackoverflow.com/questions/63507610/regarding-a-warning-message-in-jags-for-r
		samplesList<- coda::as.mcmc.list(lapply(resultat,function(x){y<-x; attributes(y)$mcpar<-c(1, attributes(y)$mcpar[2]- attributes(y)$mcpar[1]+1, attributes(y)$mcpar[3]);y}))
	}	## END: MCMC_language=="Jags"


	if (MCMC_language=="Greta") {
	  #changing nburnin.min in case MCMC_language=="Greta", for consistency of calculations of number of iterations:
	  nburnin.min0<-0
	  nburnin.min<-0
		## added to try to resolve problems with Greta
		niter<-(floor(niter/thin)+1*as.double((niter/thin-floor(niter/thin))>0))*thin
		if (!is.null(inits))
		{samplesList<-m<-greta::mcmc(model,sampler = eval(control.MCMC$sampler),initial_values=inits, warmup=control.MCMC$warmup, n_cores=control.MCMC$n_cores,chains=Nchains, n_samples=niter, thin=thin, pb_update=(floor(50/thin)+1*as.double((50/thin-floor(50/thin))>0))*thin)}
		else
		{samplesList<-m<-greta::mcmc(model,sampler = eval(control.MCMC$sampler),warmup=control.MCMC$warmup, n_cores=control.MCMC$n_cores,chains=Nchains, n_samples=niter, thin=thin, pb_update=(floor(50/thin)+1*as.double((50/thin-floor(50/thin))>0))*thin)}

		#to be checked whether niter is the correct argument here:
		## should be: number of MCMC samples to draw per chain (after any warm-up, but before thinning)
	}	## END: MCMC_language=="Greta"

	}) ## END: current.CPU.time<-system.time({
	time.MCMC<-(Sys.time()-time.MCMC0)
	units(time.MCMC)<-"secs"
	CPUtime.MCMC<-CPUtime.MCMC+current.CPU.time[1]+current.CPU.time[2]

	current.CPU.time<-system.time({
	numIter.samplesList<-c(numIter.samplesList,rep(thin,dim(samplesList[[1]])[1]))
	size.samplesList<-dim(samplesList[[1]])[1]
	niter.tot<-niter.tot+niter

	## indices.conv wil contain the indices in samplesList of the parameters in params.conv to select only those columns for diagnosing convergence & neff
	indices.conv<-NULL
	for (i in params.conv)
	{indices.conv<-c(indices.conv,which(i==dimnames(samplesList[[1]])[[2]]|regexpr(paste(i,"\\[",sep=""),dimnames(samplesList[[1]])[[2]])==1))
	}
	indices.conv<-sort(indices.conv)

	## indices.save wil contain the indices in samplesList of the parameters in params.save to select only those columns for saving. Only used in the end.
	indices.save<-NULL
	for (i in params.save)
	{indices.save<-c(indices.save,which(i==dimnames(samplesList[[1]])[[2]]|regexpr(paste(i,"\\[",sep=""),dimnames(samplesList[[1]])[[2]])==1))
	}
	indices.save<-sort(indices.save)


	########################################
	###1.2- removing chains in case some chains get stuck
	########################################

	#variable that will contain the updated number of (valid) MCMC chains
	Nchains.updated<-Nchains
	if (control$remove.fixedchains)
		{chains.to.remove<-(1:Nchains)[apply(sapply(samplesList,function(x){apply(x,2,var)}),2,sum)==0]
		Nchains.updated<-Nchains-length(chains.to.remove)
		if (Nchains.updated==0)
			{stop("All MCMC chains have fixed parameters. Consider revising the initial values or the model.")
			}
		if (control$convtype=="Gelman"&Nchains.updated==1) {stop("Impossible to use the Gelman-Rubin convergence diagnostic with only one MCMC chain (***after update of Nchains***): please change Nchain or control$convtype")}

		}



	########################################
	###1.3- first diagnostic
	########################################

	index.conv0<-index.conv
	index.props.conv<-1
	if (length(chains.to.remove)>0)
	{index.conv.local<-index.conv
	diags<-conveff(window(coda::as.mcmc.list(samplesList[-chains.to.remove])[,indices.conv],start=index.conv,end=size.samplesList),control$print.diagnostics)}
	else
	{index.conv.local<-index.conv
	diags<-conveff(window(samplesList[,indices.conv],start=index.conv,end=size.samplesList),control$print.diagnostics)}

	########################################
	###1.4- interpreting convergence of first run:
	########################################

	index.conv.new<-floor(quantile(index.conv0:size.samplesList,control$props.conv[index.props.conv]))
	##in case has not converged, re-checks convergence at index.conv.new if acceptable in terms of
		## (i) nburnin.max, (ii) control$props.conv[index.props.conv] & (iii) number of remaining values (must be >control$min.Nvalues)
	while (!checking.convergence(diags,control$convtype,control$convtype.Gelman,conv.max,conv.med,conv.mean)&((nburnin.min+sum(numIter.samplesList[1:(index.conv.new-1)]))<nburnin.max)&(control$props.conv[index.props.conv]<1)&((length(index.conv.new:size.samplesList)*Nchains.updated)>control$min.Nvalues))
	{
		index.conv.local<-index.conv.new
		nburnin<-nburnin.min+sum(numIter.samplesList[1:(index.conv.local-1)])
		if (length(chains.to.remove)>0)
			{diags<-conveff(window(coda::as.mcmc.list(samplesList[-chains.to.remove])[,indices.conv],start=index.conv.local,end=size.samplesList),control$print.diagnostics)}
		else
			{diags<-conveff(window(samplesList[,indices.conv],start=index.conv.local,end=size.samplesList),control$print.diagnostics)}

		index.props.conv<-index.props.conv+1
		index.conv.new<-floor(quantile(index.conv0:size.samplesList,control$props.conv[index.props.conv]))
	}
	##final convergence diagnostic:
	converged<-checking.convergence(diags,control$convtype,control$convtype.Gelman,conv.max,conv.med,conv.mean)
	previously.converged<-as.logical(max(previously.converged,converged))

	## changing converged state in cases where we should not have checked convergence at first run, not checked convergence at all or if iterations above nburnin.max
	if (!control$check.convergence.firstrun) {converged<-FALSE}
	if (!control$check.convergence) {converged<-TRUE}
	if (converged) {index.conv<-index.conv.local}

	### FG:  inactivated 2023/03/11: I don't understand it : if ((nburnin.min+sum(numIter.samplesList[1:(index.conv-1)]))>=nburnin.max) {converged<-TRUE}

	########################################
	###1.5 in case of convergence, calculation of thinmult for future update of thin
	########################################

	thinmult<-1
	thinmult0<-1
	if (converged)
	{
		if ((length(index.conv:size.samplesList)*Nchains.updated)>control$min.Nvalues) {
			## first checking if neffs reached
			neffs.reached<-checking.neffs.reached(diags,neff.min,neff.med,neff.mean)

			## second, if neffs not reached, calculates multiplier of thin to be targeted
			if (!neffs.reached)
			{
				thinmult<-calculate.thinmult.target(diags,neff.min,neff.med,neff.mean)
				if (control$print.thinmult)
				{
					print(paste("raw multiplier of thin: ", formatC_adapted(thinmult)))
					print("###################################################################################")
				}
				if (thinmult<control$min.thinmult){thinmult<-1}
				if (control$round.thinmult) {thinmult<-round(thinmult)}
			}

		}

	}


	##duplicating samplesList to sampledList.temp
	if (length(chains.to.remove)>0)
		{samplesList.temp<-coda::as.mcmc.list(samplesList[-chains.to.remove])}
	else
		{samplesList.temp<-samplesList}

	size.samplesList.temp<-size.samplesList
	index.conv.temp<-index.conv
	}) ## END: current.CPU.time<-system.time({
	CPUtime.btadjust<-CPUtime.btadjust+current.CPU.time[1]+current.CPU.time[2]



	########################################
	## 2. Rerunning the model sequentially:
	########################################
	#default value for previously.converged0
	previously.converged0<-previously.converged

	while((!converged |!neffs.reached)&niter>=thin)
	{
	current.CPU.time<-system.time({
	##store previously.converged at the beginning of the last cycle
	previously.converged0<-previously.converged

	########################################
	### 2.1: specifying the new thinning level and the new number of iterations
	########################################

		thin<-min(thin*thinmult,thin.max)

		#estimation of number of efficient values already available
		if (!converged|((length(index.conv:size.samplesList)*Nchains.updated)<=control$min.Nvalues))
		{available.neffs<-0}
		else {available.neffs<-scale.available.neffs(diags,neff.min,neff.med,neff.mean)}

		if (!is.null(control$time.max)) {
			current.time<-Sys.time()
			duration<-(current.time-time.start)
			units(duration)<-"secs"
			niter.previous<-niter
			niter<--1
			duration<-as.double(duration)

			##if not converged: then roughly plan to redo niter.tot iterations - under conditions
				## implies that at the end of next run niter.tot would be approximately doubled
				## except especially if control$time.max might be exceeded
			if (! converged) {
				niter<-ceiling(min(c(niter.tot,niter.max-niter.tot,(control$time.max-duration)*0.95*niter.tot/duration)))
				if (!previously.converged)
				{
				niter<-ceiling(min(c(niter.tot,niter.max-niter.tot,(control$time.max-duration)*0.95*niter.tot/duration)))
				} else {
				#niter<-ceiling(min(c(niter.previous,niter.max-niter.tot,(control$time.max-duration)*0.95*niter.tot/duration)))
				niter<-ceiling(min(c(round(niter.tot/Ncycles),niter.max-niter.tot,(control$time.max-duration)*0.95*niter.tot/duration)))
				}
				print("Case of niter update: Non convergence")
				print("###################################################################################")
			}

			## if converged but not enough values to estimate thinmult
				### then add few iterations to just exceed control$min.Nvalues while fulfilling nitex.max & control$time.max conditions
			if ( converged & ((length(index.conv:size.samplesList)*Nchains.updated)<=control$min.Nvalues))
			{
				niter<-ceiling(min(c((control$min.Nvalues-(length(index.conv:size.samplesList)*Nchains.updated)/Nchains.updated)+10,niter.max-niter.tot,(control$time.max-duration)*0.95*niter.tot/duration)))
				print("Case of niter update: Convergence but not enough values after convergence to safely update thin")
				print("###################################################################################")
			}

			## if converged and enough values to estimate thinmult
			if ( converged & ((length(index.conv:size.samplesList)*Nchains.updated)>control$min.Nvalues))
			{


				if ((Ncycles+1)/control$Ncycles.target<0.95)
				{
					## if (Ncycles+1)/control$Ncycles.target<0.95
					## then do not plan to reach the number of efficient values in the next Cycle while fulfilling nitex.max & control$time.max conditions
					niter<-ceiling((Ncycles+1)/control$Ncycles.target*min(c(ceiling((control$safemultiplier.Nvals*neff.max-available.neffs)/Nchains.updated)*thin,niter.max-niter.tot,(control$time.max-duration)*0.95*niter.tot/duration)))
					print("Case of niter update: Convergence but not in the last planned cycle")
					print("###################################################################################")
				}
				else {
					## if (Ncycles+1)/control$Ncycles.target>=0.95
					## then plan to reach the number of efficient values in the next Cycle while fulfilling nitex.max & control$time.max conditions
					niter<-min(c(ceiling((control$safemultiplier.Nvals*neff.max-available.neffs)/Nchains.updated)*thin,niter.max-niter.tot,(control$time.max-duration)*0.95*niter.tot/duration))
					print("Case of niter update: Convergence and trying to reach end of MCMC at the end of next cycle")
					print("###################################################################################")
				}
			}
			## corrections of values of niter in case these values are inadequate
			if (niter<=10) {niter<-ceiling(min(niter.tot,niter.tot/duration*(0.95*control$time.max-duration)))}
			#if (niter==Inf){niter<-niter.tot}
			if (niter==Inf){ceiling(min(c(niter.tot,niter.max-niter.tot,(control$time.max-duration)*0.95*niter.tot/duration)))}
			if (niter<=0)
			{
				print("Number of planned new iterations non-positive: end of MCMC cycles")
				print("###################################################################################");
				niter<--1
			}

		}	##END if (!is.null(control$time.max))
		else
		{#is.null(control$time.max); thus niter.max is specified
			niter.previous<-niter
			niter<--1
			##if not converged: then roughly plan to redo niter.tot iterations - under conditions - except if previously.converged in which case niter is essentially left unchanged
				## implies that at the end of next run niter.tot would be approximately doubled
				## except especially if control$time.max might be exceeded
			if (! converged)
			{if (!previously.converged)
				{
				niter<-ceiling(min(c(niter.tot,niter.max-niter.tot)))
				} else {
				#niter<-ceiling(min(c(niter.previous,niter.max-niter.tot)))
				niter<-ceiling(min(c(round(niter.tot/Ncycles),niter.max-niter.tot)))

				}
				print("Case of niter update: Non convergence")
				#print(paste0("new niter: ",niter,", previously.converged:",previously.converged))
				print("###################################################################################")
			}

			## if converged but not enough values to estimate thinmult
				### then add few iterations to just exceed control$min.Nvalues while fulfilling nitex.max & control$time.max conditions
			if ( converged & ((length(index.conv:size.samplesList)*Nchains.updated)<=control$min.Nvalues))
			{
				#print(paste("case 2",length(index.conv:size.samplesList),niter.max-niter.tot))
				niter<-ceiling(min(c((control$min.Nvalues-(length(index.conv:size.samplesList)*Nchains.updated)/Nchains.updated)+10,niter.max-niter.tot)))
				print("Case of niter update: Convergence but not enough values after convergence to safely update thin")
				print("###################################################################################")
			}

			if ( converged & ((length(index.conv:size.samplesList)*Nchains.updated)>control$min.Nvalues))
			{
				if (((Ncycles+1)/control$Ncycles.target)<0.95)
				{
					## if (Ncycles+1)/control$Ncycles.target<0.95
					## then do not plan to reach the number of efficient values in the next Cycle while fulfilling nitex.max & control$time.max conditions
					niter<-ceiling((Ncycles+1)/control$Ncycles.target*min(c(ceiling((control$safemultiplier.Nvals*neff.max-available.neffs)/Nchains.updated)*thin,niter.max-niter.tot)))
					print("Case of niter update: Convergence but not in the last planned cycle")
					print("###################################################################################")
				}
				else {
					## if (Ncycles+1)/control$Ncycles.target>=0.95
					## then plan to reach the number of efficient values in the next Cycle while fulfilling nitex.max & control$time.max conditions
						#print(paste("case 5bis",available.neffs,diags$Nvals,niter.max-niter.tot, thin))
					niter<-min(c(ceiling((control$safemultiplier.Nvals*neff.max-available.neffs)/Nchains.updated)*thin,niter.max-niter.tot))
					print("Case of niter update: Convergence and trying to reach end of MCMC at the end of next cycle")
					print("###################################################################################")
				}
			}

			## corrections of values of niter in case these values are inadequate
			if (niter<=10)
			{
				#print("case 6")
				niter<-ceiling(min(niter.tot,niter.max-niter.tot))
			}
			if (niter==Inf){niter<-niter.tot}
			if (niter<=0)
			{
				print("Number of planned new iterations non-positive: end of MCMC cycles")
				print("###################################################################################")
				niter<--1
			}
		}
	}) ## END: current.CPU.time<-system.time({
	CPUtime.btadjust<-CPUtime.btadjust+current.CPU.time[1]+current.CPU.time[2]

		########################################
		### 2.2: rerunning the model:
		########################################
		if (niter>=thin)
		{
			current.CPU.time<-system.time({
			time.MCMC0<-Sys.time()
			##slight modifications of niter & thin prior to re-running models
			niter<-ceiling(niter)
			thin<-ceiling(thin)
			## added to try to resolve problems with Greta:
				## cf. https://forum.greta-stats.org/t/problem-with-extra-samples-command-when-thin-pb-update-and-or-n-samples-are-not-multiples/338
			if (MCMC_language=="Greta") {niter<-ceiling((floor(niter/thin)+1*as.double((niter/thin-floor(niter/thin))>0))*thin)}

			##updating Ncycles
			Ncycles<-Ncycles+1
			message(control$identifier.to.print,"Cycle ", Ncycles, "...")
			print("###################################################################################")

			if (MCMC_language=="Nimble")
			{
				for (i in 1:Nchains)
				{
					message("      Running chain ", i, ", niter: ", niter, "...")
						#message("dimension au début: ", dim(as.matrix(CModelMCMC[[i]]$mvSamples)))
					set.seed(control$seed+i+Nchains)
					CModelMCMC[[i]]$run(niter, thin = thin, progressBar =TRUE,reset = FALSE)
					samplesList[[i]] <- as.matrix(CModelMCMC[[i]]$mvSamples)
						#message("dimension à la fin: ", dim(as.matrix(CModelMCMC[[i]]$mvSamples)))
				}
				samplesList <- coda::as.mcmc.list(lapply(samplesList, coda::as.mcmc))
			}	## END: MCMC_language=="Nimble"

			if (MCMC_language=="Jags")
			{
				resultatt<-rjags::coda.samples(myModel,var=params, n.iter = niter, thin = thin)
				matres <- rbind(matres,list.as.matrixbisdim1(resultatt))
				samplesList<- coda::as.mcmc.list(lapply(seq_len(length(samplesList)),function(x){coda::as.mcmc(rbind(as.matrix(samplesList[[x]]),as.matrix(resultatt[[x]])))}))
			}		## END: MCMC_language=="Jags"

			if (MCMC_language=="Greta")
			{
				samplesList<-m<-greta::extra_samples(m, n_samples=niter, thin=thin, n_cores=control.MCMC$n_cores, pb_update=(floor(50/thin)+1*as.double((50/thin-floor(50/thin))>0))*thin)

			}		## END: MCMC_language=="Greta"
			}) ## END: current.CPU.time<-system.time({
			##updating time.MCMC
			time.MCMC.duration<-(Sys.time()-time.MCMC0)
			units(time.MCMC.duration)<-"secs"
			time.MCMC<-time.MCMC+time.MCMC.duration
			CPUtime.MCMC<-CPUtime.MCMC+current.CPU.time[1]+current.CPU.time[2]

			current.CPU.time<-system.time({
			##updating counters around samplesList & niter
			numIter.samplesList<-c(numIter.samplesList,rep(thin,dim(samplesList[[1]])[1]-size.samplesList))
			size.samplesList<-dim(samplesList[[1]])[1]
			niter.tot<-niter.tot+niter


			########################################
			###2.3- Reshaping samplesList so that the spacing between values is roughly of thin: result in samplesList.temp
			########################################

			index.conv0<-1

			### reversed vector of - decreasing- indices to save
			temp<-length(numIter.samplesList[index.conv:size.samplesList])-round(cumsum(rev(numIter.samplesList[index.conv:size.samplesList]))/thin)+1
			### same but not rounded
			tempunr<-length(numIter.samplesList[index.conv:size.samplesList])-(cumsum(rev(numIter.samplesList[index.conv:size.samplesList]))/thin)+1
			## unique values of temp; will be changed afterwards
			indices.samplesList<-unique(temp)
			## same, will not be changed afterwards
			indices.samplesList0<- indices.samplesList

			for (i in indices.samplesList)
			{
				tempbis<-which(temp==i)
				## will give the corresponding order on the scale of indices of [index.conv:size.samplesList]
				indices.samplesList[indices.samplesList0==i]<-length(numIter.samplesList[index.conv:size.samplesList])-tempbis[which.min(abs(tempunr[tempbis]-i))[1]]+1
			}
			## transforming indices on the scale of indices of numIter.samplesList
			indices.samplesList<-(index.conv:size.samplesList)[sort(indices.samplesList)]

			## tranferring the part of samplesList corresponding to indices.samplesList to samplesList.temp
			samplesList.temp<-samplesList
			## associated index.conv is 1 by definition
			index.conv.temp<-1
			index.conv0.temp<-index.conv.temp
			for (i in 1:Nchains)
			{
				samplesList.temp[[i]]<-samplesList[[i]][indices.samplesList,]
			}

			samplesList.temp <- coda::as.mcmc.list(lapply(samplesList.temp, coda::as.mcmc))

			if (length(chains.to.remove)>0)
				{samplesList.temp<-coda::as.mcmc.list(samplesList.temp[-chains.to.remove])}
			else
				{samplesList.temp<-samplesList.temp}
			size.samplesList.temp<-dim(samplesList.temp[[1]])[1]


			########################################
			###2.4- diagnostics based on new model
			########################################
			index.conv.local<-indices.samplesList[index.conv.temp]
			diags<-conveff(window(samplesList.temp[,indices.conv],start=index.conv.temp,end=size.samplesList.temp),control$print.diagnostics)

			########################################
			###2.4.1- interpreting convergence of next runs (done only if not already converged or if we have to control$recheck.convergence):
			########################################

			if (!converged|control$recheck.convergence)
			{
				index.props.conv<-1
				index.conv.new<-floor(quantile(index.conv0.temp:size.samplesList.temp,control$props.conv[index.props.conv]))
				index.conv.new0<-indices.samplesList[index.conv.new]
				##in case has not converged, re-checks convergence at index.conv.new if acceptable in terms of
					## (i) nburnin.max, (ii) control$props.conv[index.props.conv] & (iii) number of remaining values (must be >control$min.Nvalues)
				while (!checking.convergence(diags,control$convtype,control$convtype.Gelman,conv.max,conv.med,conv.mean)&((nburnin.min+sum(numIter.samplesList[1:(index.conv.new0-1)]))<nburnin.max)&(control$props.conv[index.props.conv]<1)&((length(index.conv.new:size.samplesList.temp)*Nchains.updated)>control$min.Nvalues))
				{
					index.conv.temp<-index.conv.new
					nburnin<-nburnin.min+sum(numIter.samplesList[1:(index.conv-1)])
					index.conv.local<-indices.samplesList[index.conv.temp]
					diags<-conveff(window(samplesList.temp[,indices.conv],start=index.conv.temp,end=size.samplesList.temp),control$print.diagnostics)
					index.props.conv<-index.props.conv+1
					index.conv.new<-floor(quantile(index.conv0.temp:size.samplesList.temp,control$props.conv[index.props.conv]))
					index.conv.new0<-indices.samplesList[index.conv.new]
				}

				##final convergence diagnostic:
				converged<-checking.convergence(diags,control$convtype,control$convtype.Gelman,conv.max,conv.med,conv.mean)
				previously.converged<-as.logical(max(previously.converged,converged))
				if (converged) {index.conv<-indices.samplesList[index.conv.temp]}
				## changing converged state in cases where we should not have checked convergence at first run, not checked convergence at all or if iterations above nburnin.max
				###FG: inactivated 2023/03/11: I don't understand it : if ((nburnin.min+sum(numIter.samplesList[1:(index.conv-1)]))>=nburnin.max) {converged<-TRUE}
			}	## END if (!converged|control$recheck.convergence)

			########################################
			###2.4.2. in case of convergence, calculating thinmult for future udate of thin
			########################################

			thinmult<-1
			if (converged)
			{
				if ((length(index.conv.temp:size.samplesList.temp)*Nchains.updated)>control$min.Nvalues)
				{
					## first checking if neffs reached
					neffs.reached<-checking.neffs.reached(diags,neff.min,neff.med,neff.mean)
					## second, if neffs not reached, calculates multiplier of thin to be targeted
					if (!neffs.reached)
					{
						thinmult<-calculate.thinmult.target(diags,neff.min,neff.med,neff.mean)
						if (control$print.thinmult)
						{
							print(paste("raw multiplier of thin: ", formatC_adapted(thinmult)))
							print("###################################################################################")
						}
						if (thinmult<control$min.thinmult){thinmult<-1}
						if (control$round.thinmult) {thinmult<-round(thinmult)}
					}
				}

			}

			#print(paste("underlying indices: ",index.conv.temp,", ",size.samplesList.temp))
			}) ## END: current.CPU.time<-system.time({
			CPUtime.btadjust<-CPUtime.btadjust+current.CPU.time[1]+current.CPU.time[2]
			} ## END: if (niter>=thin)
	}	##END of while((!converged |!neffs.reached)&niter>=thin)


	########################################
	###2.5- reshaping samplesList in case converged & !previously.converged0: otherwise: much too big object
	########################################
	current.CPU.time<-system.time({
	if (converged & !previously.converged0) {

		########################################
		### 2.5.1: specifying the new thinning level
		########################################

			thin<-min(thin*thinmult,thin.max)

		########################################
		###  2.5.2- Reshaping samplesList so that the spacing between values is roughly of thin: result in samplesList.temp
		########################################

				index.conv0<-1

				### reversed vector of - decreasing- indices to save
				temp<-length(numIter.samplesList[index.conv:size.samplesList])-round(cumsum(rev(numIter.samplesList[index.conv:size.samplesList]))/thin)+1
				### same but not rounded
				tempunr<-length(numIter.samplesList[index.conv:size.samplesList])-(cumsum(rev(numIter.samplesList[index.conv:size.samplesList]))/thin)+1
				## unique values of temp; will be changed afterwards
				indices.samplesList<-unique(temp)
				## same, will not be changed afterwards
				indices.samplesList0<- indices.samplesList

				for (i in indices.samplesList)
				{
					tempbis<-which(temp==i)
					## will give the corresponding order on the scale of indices of [index.conv:size.samplesList]
					indices.samplesList[indices.samplesList0==i]<-length(numIter.samplesList[index.conv:size.samplesList])-tempbis[which.min(abs(tempunr[tempbis]-i))[1]]+1
				}
				## transforming indices on the scale of indices of numIter.samplesList
				indices.samplesList<-(index.conv:size.samplesList)[sort(indices.samplesList)]

				## tranferring the part of samplesList corresponding to indices.samplesList to samplesList.temp
				samplesList.temp<-samplesList
				## associated index.conv is 1 by definition
				index.conv.temp<-1
				index.conv0.temp<-index.conv.temp
				for (i in 1:Nchains)
				{
					samplesList.temp[[i]]<-samplesList[[i]][indices.samplesList,]
				}

				samplesList.temp <- coda::as.mcmc.list(lapply(samplesList.temp, coda::as.mcmc))

				if (length(chains.to.remove)>0)
					{samplesList.temp<-coda::as.mcmc.list(samplesList.temp[-chains.to.remove])}
				else
					{samplesList.temp<-samplesList.temp}
				size.samplesList.temp<-dim(samplesList.temp[[1]])[1]


	}

	if (MCMC_language=="Nimble") {
	########################################
	###3. uncompiling Nimble objects: https://groups.google.com/g/nimble-users/c/-eoYs__eg0o
	### Only useful outside Windows
	########################################
		### Only useful if MCMC_language=="Nimble"
	  if (Sys.info()["sysname"]!="Windows") {
			for (i in 1:Nchains)
			{try(nimble::clearCompiled(CModelMCMC[[i]]),silent=TRUE)
			try(nimble::clearCompiled(Model[[i]]),silent=TRUE)
			}
	  try(nimble::clearCompiled(Modeltemp),silent=TRUE)
	  try({rm(CModelMCMC)})
	  try({rm(Modeltemp)})
	  }
	}	## END: MCMC_language=="Nimble"

	########################################
	####4. changing names of MCMC parameters
	########################################

	### this is done so that names of parameters are coherent between MCMC.languages
	### there is especially a discrepancy between Nimble & Greta in the names of elements of matrices (adding or not a space after commas)

	for (i in seq_len(length(samplesList.temp)))
	{
		dimnames(samplesList.temp[[i]])[[2]]<-gsub(", ",",",dimnames(samplesList.temp[[i]])[[2]])
	}


	########################################
	####5. assigning the results
	########################################

	if (!converged) warning("The MCMC did not converge")
	if (!neffs.reached) warning("The expected effective sample size was not reached")


	if (converged)
	  {result<-window(samplesList.temp[,indices.save],start=index.conv.temp,end=size.samplesList.temp)} else
	  {result<-samplesList.temp[,indices.save]}
	}) ## END: current.CPU.time<-system.time({
	CPUtime.btadjust<-CPUtime.btadjust+current.CPU.time[1]+current.CPU.time[2]
	total.duration<-Sys.time()-time.start
	units(total.duration)<-"secs"
	original.atrributes<-list(call.params=list(summarized.data= {if(!is.null(data)) {summarize.data(data)} else {NULL} },
												summarized.consts= {if(!is.null(constants)) {summarize.data(constants)} else {NULL} },
												params=params,params.conv=params.conv,params.save=params.save,
												niter.min=niter.min,niter.max=niter.max,nburnin.min=nburnin.min,nburnin.max=nburnin.max,thin.min=thin.min,thin.max=thin.max,
												time.max=control$time.max,
												inits=inits,
												neff.min=neff.min,neff.med=neff.med,neff.mean=neff.mean,
												conv.max=conv.max,conv.med=conv.med,conv.mean=conv.mean,
												control=control,control.MCMC=control.MCMC,
												Nchains=Nchains,final.Nchains=Nchains.updated),
							final.params=list(converged=converged,burnin=nburnin.min0+sum(numIter.samplesList[1:ifelse(converged,index.conv-1,size.samplesList)]),thin=thin,niter.tot=niter.tot,
											duration=total.duration,duration.MCMC.preparation=time.MCMC.Preparation,
											duration.MCMC.transient=time.MCMC*(ifelse(MCMC_language=="Greta",1,0)*control.MCMC$warmup+ifelse(MCMC_language=="Jags",1,0)*control.MCMC$n.adapt+nburnin.min0+sum(numIter.samplesList[1:ifelse(converged,index.conv-1,size.samplesList)]))/(ifelse(MCMC_language=="Greta",1,0)*control.MCMC$warmup+ifelse(MCMC_language=="Jags",1,0)*control.MCMC$n.adapt+nburnin.min0+sum(numIter.samplesList)),
											duration.MCMC.asymptotic=time.MCMC*(1-(ifelse(MCMC_language=="Greta",1,0)*control.MCMC$warmup+ifelse(MCMC_language=="Jags",1,0)*control.MCMC$n.adapt+nburnin.min0+sum(numIter.samplesList[1:ifelse(converged,index.conv-1,size.samplesList)]))/(ifelse(MCMC_language=="Greta",1,0)*control.MCMC$warmup+ifelse(MCMC_language=="Jags",1,0)*control.MCMC$n.adapt+nburnin.min0+sum(numIter.samplesList))),
											duration.btadjust=total.duration-time.MCMC.Preparation-time.MCMC,
											CPUduration=unname(CPUtime.MCMC.Preparation+CPUtime.MCMC+CPUtime.btadjust),
											CPUduration.MCMC.preparation=unname(CPUtime.MCMC.Preparation),
											CPUduration.MCMC.transient=unname(CPUtime.MCMC*(ifelse(MCMC_language=="Greta",1,0)*control.MCMC$warmup+ifelse(MCMC_language=="Jags",1,0)*control.MCMC$n.adapt+nburnin.min0+sum(numIter.samplesList[1:ifelse(converged,index.conv-1,size.samplesList)]))/(ifelse(MCMC_language=="Greta",1,0)*control.MCMC$warmup+ifelse(MCMC_language=="Jags",1,0)*control.MCMC$n.adapt+nburnin.min0+sum(numIter.samplesList))),
											CPUduration.MCMC.asymptotic=unname(CPUtime.MCMC*(1-(ifelse(MCMC_language=="Greta",1,0)*control.MCMC$warmup+ifelse(MCMC_language=="Jags",1,0)*control.MCMC$n.adapt+nburnin.min0+sum(numIter.samplesList[1:ifelse(converged,index.conv-1,size.samplesList)]))/(ifelse(MCMC_language=="Greta",1,0)*control.MCMC$warmup+ifelse(MCMC_language=="Jags",1,0)*control.MCMC$n.adapt+nburnin.min0+sum(numIter.samplesList)))),
											CPUduration.btadjust=unname(CPUtime.btadjust),
											time_end=Sys.time()),
							#NB: burnin and niter.tot are in MCMC iteration units for one chain; thin is in MCMC iteration units
							final.diags=if (converged){index.conv.local<-index.conv.temp;conveff_final(window(samplesList.temp[,indices.conv],start=index.conv.temp,end=size.samplesList.temp))} else {index.conv.local<-index.conv.temp;conveff_final(window(samplesList.temp[,indices.conv]))},
							package.versions=sapply(sessionInfo()$otherPkgs,function(x){x$Version}),
							R.version=sessionInfo()$R.version,
							# removed because a priori useless outside runMCMC_btadjust							numIter.samplesList=numIter.samplesList,
							warnings=c(NULL))
	attributes(result)<-original.atrributes
	result<-coda::as.mcmc.list(result)
	## if the MCMC has not converged: risk that there are much too many values. We
	if (!converged)
		{current.result.size<-length(result)*dim(result[[1]])[1]
		thin.end<-floor(current.result.size/(2*min(c(neff.min,neff.med,neff.mean),na.rm=TRUE)))
		if (thin.end>1)
			{
				result<-window(result,thin=thin.end)
				new.atrributes<-list(call.params=list(summarized.data= {if(!is.null(data)) {summarize.data(data)} else {NULL} },
														summarized.consts= {if(!is.null(constants)) {summarize.data(constants)} else {NULL} },
														params=params,params.conv=params.conv,params.save=params.save,
														niter.min=niter.min,niter.max=niter.max,nburnin.min=nburnin.min,nburnin.max=nburnin.max,thin.min=thin.min,thin.max=thin.max,
														time.max=control$time.max,
														inits=inits,
														neff.min=neff.min,neff.med=neff.med,neff.mean=neff.mean,
														conv.max=conv.max,conv.med=conv.med,conv.mean=conv.mean,
														control=control,control.MCMC=control.MCMC,
														Nchains=Nchains,final.Nchains=Nchains.updated),
									final.params=list(converged=converged,burnin=nburnin.min0+sum(numIter.samplesList[1:ifelse(converged,index.conv-1,size.samplesList)]),thin=thin*thin.end,niter.tot=niter.tot,
													duration=total.duration,duration.MCMC.preparation=time.MCMC.Preparation,
													duration.MCMC.transient=time.MCMC*(ifelse(MCMC_language=="Greta",1,0)*control.MCMC$warmup+ifelse(MCMC_language=="Jags",1,0)*control.MCMC$n.adapt+nburnin.min0+sum(numIter.samplesList[1:ifelse(converged,index.conv-1,size.samplesList)]))/(ifelse(MCMC_language=="Greta",1,0)*control.MCMC$warmup+ifelse(MCMC_language=="Jags",1,0)*control.MCMC$n.adapt+nburnin.min0+sum(numIter.samplesList)),
													duration.MCMC.asymptotic=time.MCMC*(1-(ifelse(MCMC_language=="Greta",1,0)*control.MCMC$warmup+ifelse(MCMC_language=="Jags",1,0)*control.MCMC$n.adapt+nburnin.min0+sum(numIter.samplesList[1:ifelse(converged,index.conv-1,size.samplesList)]))/(ifelse(MCMC_language=="Greta",1,0)*control.MCMC$warmup+ifelse(MCMC_language=="Jags",1,0)*control.MCMC$n.adapt+nburnin.min0+sum(numIter.samplesList))),
													duration.btadjust=total.duration-time.MCMC.Preparation-time.MCMC,
													CPUduration=unname(CPUtime.MCMC.Preparation+CPUtime.MCMC+CPUtime.btadjust),
													CPUduration.MCMC.preparation=unname(CPUtime.MCMC.Preparation),
													CPUduration.MCMC.transient=unname(CPUtime.MCMC*(ifelse(MCMC_language=="Greta",1,0)*control.MCMC$warmup+ifelse(MCMC_language=="Jags",1,0)*control.MCMC$n.adapt+nburnin.min0+sum(numIter.samplesList[1:ifelse(converged,index.conv-1,size.samplesList)]))/(ifelse(MCMC_language=="Greta",1,0)*control.MCMC$warmup+ifelse(MCMC_language=="Jags",1,0)*control.MCMC$n.adapt+nburnin.min0+sum(numIter.samplesList))),
													CPUduration.MCMC.asymptotic=unname(CPUtime.MCMC*(1-(ifelse(MCMC_language=="Greta",1,0)*control.MCMC$warmup+ifelse(MCMC_language=="Jags",1,0)*control.MCMC$n.adapt+nburnin.min0+sum(numIter.samplesList[1:ifelse(converged,index.conv-1,size.samplesList)]))/(ifelse(MCMC_language=="Greta",1,0)*control.MCMC$warmup+ifelse(MCMC_language=="Jags",1,0)*control.MCMC$n.adapt+nburnin.min0+sum(numIter.samplesList)))),
													CPUduration.btadjust=unname(CPUtime.btadjust),
													time_end=Sys.time()),
									#NB: burnin and niter.tot are in MCMC iteration units for one chain; thin is in MCMC iteration units
									final.diags=conveff_final(result),
									package.versions=sapply(sessionInfo()$otherPkgs,function(x){x$Version}),
									R.version=sessionInfo()$R.version,
									# removed because a priori useless outside runMCMC_btadjust							numIter.samplesList=numIter.samplesList,
									warnings=c(NULL),
									original.attributes=original.atrributes)
				attributes(result)<-new.atrributes
			}
		}

	coda::as.mcmc.list(result)
}
