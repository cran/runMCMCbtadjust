% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runMCMC_btadjust.r
\name{runMCMC_btadjust}
\alias{runMCMC_btadjust}
\title{runMCMC_btadjust}
\usage{
runMCMC_btadjust(
  code = NULL,
  data = NULL,
  constants = NULL,
  model = NULL,
  MCMC_language = "Nimble",
  Nchains,
  inits = NULL,
  params = NULL,
  params.conv = NULL,
  params.save = NULL,
  niter.min = 100,
  niter.max = Inf,
  nburnin.min = 10,
  nburnin.max = Inf,
  thin.min = 1,
  thin.max = Inf,
  neff.min = NULL,
  neff.med = NULL,
  neff.mean = NULL,
  conv.max = NULL,
  conv.med = NULL,
  conv.mean = NULL,
  control = list(time.max = NULL, check.convergence = TRUE, check.convergence.firstrun =
    NULL, recheck.convergence = TRUE, convtype = NULL, convtype.Gelman = 2,
    convtype.Geweke = c(0.1, 0.5), convtype.alpha = 0.05, ip.nc = 0, neff.method =
    "Stan", Ncycles.target = 2, props.conv = c(0.25, 0.5, 0.75), min.Nvalues = NULL,
    min.thinmult = 1.1, force.niter.max = FALSE, force.time.max = FALSE,
    time.max.turns.off.niter.max = FALSE, safemultiplier.Nvals = 1.2, max.prop.decr.neff
    = 0.1, round.thinmult = TRUE, thinmult.in.resetMV.temporary = TRUE, 
    
    check.thinmult = 2, decrease.thinmult.multiplier = 0.8, decrease.thinmult.threshold =
    20, only.final.adapt.thin = FALSE, identifier.to.print = "", print.diagnostics =
    FALSE, conv.thorough.check = FALSE, print.thinmult = TRUE, innerprint = FALSE, seed =
    NULL, remove.fixedchains = TRUE, check.installation = TRUE, save.data = FALSE,
    conveff.final.allparams = TRUE),
  control.MCMC = list(confModel.expression.toadd = NULL, sampler = expression(hmc()),
    warmup = 1000, n.adapt = -1, RNG.names = c("base::Wichmann-Hill",
    "base::Marsaglia-Multicarry", "base::Super-Duper", "base::Mersenne-Twister"), n_cores
    = NULL, showCompilerOutput = FALSE, buildDerivs = FALSE, resetMV = FALSE, parallelize
    = FALSE, parallelizeInitExpr = expression(if (MCMC_language == "Nimble") {
    
    library("nimble")
     if (control.MCMC$APT) {
         library("nimbleAPT")
     }

    } else {
     NULL
 }), useConjugacy = FALSE, WAIC = FALSE, WAIC.Nsamples = 2000,
    WAIC.control = list(online = TRUE, dataGroups = NULL, marginalizeNodes = NULL,
    niterMarginal = 1000, convergenceSet = c(0.25, 0.5, 0.75), thin = TRUE, nburnin_extra
    = 0), APT = FALSE, APT.NTemps = 7, APT.initTemps = NULL, APT.tuneTemps = c(10, 0.7),
    APT.thinPrintTemps = expression(niter/5), includeAllStochNodes = FALSE,
    monitorAllStochNodes = FALSE, saveAllStochNodes = FALSE, includeParentNodes = FALSE,
    monitorParentNodes = FALSE, saveParentNodes = FALSE, 
     extraCalculations = NULL)
)
}
\arguments{
\item{code}{R object: code for the model that will be used to build the MCMC when \code{MCMC_language} is "Nimble" or "Jags". If "Nimble", must be the name (in R) of the object which is the result of the function \code{nimbleCode}. If "Jags", should be either: (i) a character string which is the name of a txt file that contains the code of the model (as used in the function jags.model): should then end up by ".txt"; or (ii) a character string that contains the text of the Jags code.}

\item{data}{R list: a list that will contain the data when \code{MCMC_language} is "Nimble" or "Jags". If "Nimble", will be sent to the \code{data} argument of the \code{nimbleModel} function in \code{nimble} package, i.e. the data that have a random distribution in the model. If \code{MCMC_language} is "Greta", can be used just to document the summary of data in the output.}

\item{constants}{R list: a list that will contain the rest of the data (in addition to data) when \code{MCMC_language} is "Nimble". Will be sent to the \code{constants} argument of the \code{nimbleModel} function in \code{nimble} package, i.e. the data that do not have a random distribution in the model. If \code{MCMC_language} is "Greta", can be used just to document the summary of other data in the output.}

\item{model}{R object: should be the result of the \code{model} command of Greta when \code{MCMC_language} is "Greta".}

\item{MCMC_language}{character value: designates the \code{MCMC_language} used to write & fit the Bayesian model in R. Current choices are "Nimble" - the default-, "Greta" or "Jags". Note that in case it is "Nimble", package \code{nimble} should be loaded in your search list if not using parallelization.}

\item{Nchains}{integer value : the number of Markov chains to run in the MCMC.}

\item{inits}{either NULL, a function or an R list, with \code{Nchains} components. Each component is then a named list that contains the initial values of the parameters for the MCMC.
In case \code{MCMC_language=="Greta"}, each component should be the result of the \code{initials} function in \code{greta} package.
If a function, it will generate values for one chain.}

\item{params}{character vector: contains the names of the parameters to save at the end of the MCMC and to monitor for convergence and effective sample size;
inactive for convergence/effective sample size if \code{params.conv} is specified;
inactive for saving if \code{params.save} is specified.}

\item{params.conv}{character vector: contains the names of the parameters to monitor for convergence and effective sample size.}

\item{params.save}{character vector: contains the names of the parameters to be saved at the end of the MCMC.}

\item{niter.min}{integer value: the minimum number of iterations in each chain of the MCMC.}

\item{niter.max}{integer value: the maximum number of iterations in each chain of the MCMC. Will stop the MCMC once the number of iterations will reach this limit.}

\item{nburnin.min}{integer value: the minimum number of burn-in (=transitory) iterations in each chain of the MCMC.}

\item{nburnin.max}{integer value: the maximum number of burn-in (=transitory) iterations in each chain of the MCMC. Will stay at this burn-in value once this limit is reached.}

\item{thin.min}{integer value: the minimum value of the thin parameter of the MCMC.}

\item{thin.max}{integer value: the maximum value of the thin parameter of the MCMC. Will stay at this thin value once this limit is reached.}

\item{neff.min}{positive real number: minimum effective sample size - over parameters used to diagnose convergence & effective sample size- , as calculated with \code{neff.method} (specified in \code{Control}).
The algorithm will not stop if the minimum number of effective values is not above this value (unless another limit - e.g. \code{niter.max} - is reached).}

\item{neff.med}{positive real number: median effective sample size - over parameters used to diagnose convergence & effective sample size- , as calculated with \code{neff.method} (specified in \code{Control}).
The algorithm will not stop if the median number of effective values is not above this value (unless another limit - e.g. \code{niter.max} - is reached).}

\item{neff.mean}{positive real number: mean effective sample size - over parameters used to diagnose convergence & effective sample size-, as calculated with \code{neff.method} (specified in \code{Control}).
The algorithm will not stop if the mean number of effective values is not above this value (unless another limit - e.g. \code{niter.max} - is reached).}

\item{conv.max}{positive real number: maximum - over parameters used to diagnose convergence & effective sample size - convergence diagnostic, as calculated with \code{convtype} method (specified in \code{Control}).
The algorithm will not stop if the maximum convergence diagnostic is not below this value (unless another limit - e.g. \code{niter.max} - is reached).}

\item{conv.med}{positive real number: median - over parameters used to diagnose convergence & effective sample size - convergence diagnostic, as calculated with \code{convtype} method(specified in \code{Control}).
The algorithm will not stop if the median convergence diagnostic is not below this value (unless another limit - e.g. \code{niter.max} - is reached).}

\item{conv.mean}{positive real number: mean - over parameters used to diagnose convergence & effective sample size - convergence diagnostic, as calculated with \code{convtype} method (specified in \code{Control}).
The algorithm will not stop if the mean convergence diagnostic is not below this value (unless another limit - e.g. \code{niter.max} - is reached).}

\item{control}{list of \code{runMCMC_btadjust} control parameters: with the following components:\cr
\itemize{ \item \code{time.max}: positive number (units: seconds): maximum time of the process in seconds, not including WAIC calculations and extra.calculations; the program will organize itself to stop (not including WAIC calculation and extra calculations as specified in \code{control.MCMC$extraCalculations}) before roughly \code{time.max}. Default to NULL, corresponding to no time constraint except in case of parallelization (see below).
 \item \code{check.convergence}: logical value: should the program check convergence at all? Default to TRUE. See Details.
 \item \code{check.convergence.firstrun}: logical value: should we check convergence after the first run? Default to NULL in which case will depend on \code{MCMC_language}: if MCMC_language!="Greta" & control.MCMC$n.adapt<=0, will be FALSE otherwise will be TRUE because warmup phase separated from the rest in case of Greta.
 \item \code{recheck.convergence}: logical value: should the algorithm recheck convergence once convergence has been found in a previous run? Default to TRUE.
 \item \code{convtype}: character or NULL value: specifies the type of convergence diagnostic used. Currently implemented: "Gelman" for original Gelman-Rubin diagnostic (only possible if \code{Nchains>=2}), "Gelman_new" for the version of the Gelman-Rubin diagnostic in the second version of "Bayesian Data Analysis" (Gelman, Carlin, Stern and Rubin)(only possible if \code{Nchains>=2}), "Geweke" for Geweke diagnostic (at present applied only in case \code{Nchains==1}) and "Heidelberger" for the reciprocal of Heidelberger-Welch first part of convergence diagnostic based on the Cramer-von Mises test statistic. If NULL (the default), chooses "Geweke" in case Nchains==1 and "Gelman" in case Nchains>1.
 \item \code{convtype.Gelman}: integer value: when \code{convtype=="Gelman"}, do we target the Point estimate diagnostic (value 1) or the Upper C.I. diagnostic (value 2). Default to 2.
 \item \code{convtype.Geweke}: real vector with two components between 0 and 1: (i) the fraction of samples to consider as the beginning of the chain (frac1 in geweke.diag); (ii) the fraction of samples to consider as the end of the chain (frac2 in \code{gewke.diag}). Default to c(0.1,0.5) as in \code{geweke.diag}.
 \item \code{convtype.alpha}: real value between 0 and 1: significance level used in case \code{convtype=="Gelman"} and \code{convtype.Gelman==2}, or \code{convtype=="Heidelberger"}
 \item \code{props.conv}: numeric vector: quantiles of number of iterations removed to recheck convergence and number of effective values (if not converged before or \code{conv.thorough.check} is TRUE). Values should be between 0 and 1.
 \item \code{ip.nc}: real value: inflexion point for log(scaleconvs-1); if (very) negative will tend to double the number of iterations in case of non convergence (i.e. add niter iterations) ; if (very) positive will tend to add niter/Ncycles iterations. Default to 0.
 \item \code{conv.thorough.check}: logical value: whether one goes through all the \code{props.conv} proportions to find the best one in terms first of convergence (just in terms of respecting the criterion) and then of neffs (if TRUE) or stops at the first \code{props.conv} corresponding to convergence or does not change if there is no convergence (if FALSE, the default).
 \item \code{neff.method}: character value: method used to calculate the effective sample sizes. Current choice between "Stan" (the default) and "Coda". If "Stan", uses the function \code{monitor} in package \code{rstan}. If "Coda", uses the function \code{effectiveSize} in package \code{coda}.
 \item \code{Ncycles.target}: integer value: targeted number of MCMC runs. Default to 2.
 \item \code{min.Nvalues}: integer value or expression giving an integer value or NULL: minimum number of values to diagnose convergence or level of autocorrelation. Default to NULL in which case is will be given the value neff.max calculated within the program.
 \item \code{round.thinmult}: logical value: should the thin multiplier be rounded to the nearest integer so that past values are precisely positioned on the modified iteration sequence? Default to TRUE. Value of FALSE may not be rigorous or may not converge well.
 \item \code{thinmult.in.resetMV.temporary}: logical value: should the thin multiplier be taken into account in resetting parameter collection in case MCMC_language is "Nimble". Important mainly if control.MCMC$WAIC is TRUE. If TRUE, resetting will be more frequent, and WAIC calculation will be longer and more rigorous. Default to TRUE.
 \item \code{check.thinmult}: integer value between 1, 2 and 3: how should we check thinmult value after thinmult calculation? If 3, it is tested whether thinmult meets specific criteria -relative to convergence reaching (i.e. if no convergence, no change), number of effective value reaching conservation, minimum number of output values - min.Nvalues- and proportional reduction of number of effective values - and if not decreased values are tested with the same criteria. If 2, the same checkings are done except the one on proportional reduction of effective values. If 1, only the min.Nvalues criterion is taken into account. Default to 2. A value of 3 should produce shorter MCMCs, more values in the output, with more autocorrelation, than a value of 1.
 \item \code{decrease.thinmult.multiplier}: positive number below 1: when adapting the proposed multiplier of thin (thinmult), the multiplier of the current thinmult used to propose a new - smaller - thinmult value, provided thinmult is above decrease.thinmult.threshold. Default to 0.8.
 \item \code{decrease.thinmult.threshold}: positive number above 3: when adapting the proposed multiplier of thin (thinmult), the threshold value for thinmult below which decreases of proposed thinmult are substractions of one unit. Default to 20.
 \item \code{only.final.adapt.thin}: logical value: should the thin parameter be adapted only at the end - so that during running of the MCMC we conserve a sufficient number of values - esp. with respect to min.Nvalues. Default to FALSE.
 \item \code{min.thinmult}: numeric value: minimum value of thin multiplier: if diagnostics suggest to multiply by less than this, this multiplication is not done. Default to 1.1.
 \item \code{force.niter.max}: logical value: if TRUE, the number of iterations is forced to go to niter.max - except for time.max constraints. Default to FALSE.
 \item \code{force.time.max}: logical value: if TRUE, the number of iterations is forced to go to approximately time.max seconds - except for niter.max constraints. Default to FALSE.
 \item \code{time.max.turns.off.niter.max}: logical value: if TRUE, and if time.max is finite, the number of iterations can be greater than niter.max (except if we are only in the phase of force.niter.max or force.time.max). Default to FALSE.
 \item \code{seed}: integer number or NULL value: seed for the pseudo-random number generator inside runMCMC_btadjust. Default to NULL in which case here is no control of this seed.
 \item \code{identifier.to.print}: character string: printed each time an MCMC update is ran to identify the model (esp. if multiple successive calls to \code{runMCMC_btadjust} are made).
 \item \code{safemultiplier.Nvals}: positive number: number bigger than 1 used to multiply the targeted number of effective values in calculations of additional number of iterations. Default to 1.2.
 \item \code{max.prop.decr.neff}: number between 0 and 1 used, if check.thinmult==3, to decide if we accept dimension reduction - through augmentation of thin parameter with thinmult -: maximum acceptable Proportional Decrease of the number of effective values: guaranties that at least (1-max.prop.decr.neff) times the number of effective values estimated prior to dimension reduction are kept. Default to 0.1.
 \item \code{print.diagnostics}: logical value: should diagnostics be printed each time they are calculated? Default to FALSE.
 \item \code{print.thinmult}: logical value: should the raw multiplier of thin be printed each time it is calculated? Default to TRUE.
 \item \code{innerprint}: logical value: should printings be done inside the function \code{monitor} of \code{rstan} in case \code{neff.method=="Stan"}? Default to FALSE.
 \item \code{remove.fixedchains}: logical value: should we remove Markov chains that do not vary during the first cycle (i.e. whose all parameters have zero variances)? Default to TRUE. If so, those chains are removed from the diagnostics, from the output values and associated chains are no longer updated if MCMC_language=="Nimble".
 \item \code{check.installation}: logical value: should the function check installation of packages and programs? Default to TRUE.
 \item \code{save.data}: logical value: should the program save the entire data in the call.params section of the attributes? Default to FALSE, in which case only a summary of data is saved.
 \item \code{conveff.final.allparams}: logical value: should the final convergence/number of effective values calculations in final.diags be done on all parameters? Default to TRUE.
  }}

\item{control.MCMC}{list of MCMC control parameters: with the following components - that depend on \code{MCMC_language}:
 \itemize{ \item \code{confModel.expression.toadd} (only for \code{MCMC_language=="Nimble"}): expression to add to \code{confModel} to specify samplers, remove nodes... \code{confModel} should be referred to by \code{confModel[[i]]}. See Details for an example.
\item \code{sampler} (only for \code{MCMC_language=="Greta"}): expression used to specify the sampler used.
\item \code{warmup} (only for \code{MCMC_language=="Greta"}): integer value used as warmup parameter in the mcmc.
\item \code{n.adapt} (only for \code{MCMC_language=="Jags"} or \code{MCMC_language=="Nimble"}): integer value: number of iterations used for adaptation (in function \code{jags.model} in \code{rjags} package in case \code{MCMC_language=="Jags"} and otherwise in Nimble - added to burnin: first iterations that are not saved).
\item \code{RNG.names} (only for \code{MCMC_language=="Jags"}): character vector: name of pseudo-random number generators for each chain. Each component of the vector should be among "base::Wichmann-Hill", "base::Marsaglia-Multicarry", "base::Super-Duper", "base::Mersenne-Twister". If less values than \code{Nchains} are provided, they are specified periodically.
\item \code{n_cores} (only for \code{MCMC_language=="Greta"}): integer or NULL: maximum number of cores to use by each sampler.
\item \code{showCompilerOutput} (only for \code{MCMC_language=="Nimble"}): logical value indicating whether details of C++ compilation should be printed. Default to FALSE to get default printings of limited size.
\item \code{buildDerivs} (only for \code{MCMC_language=="Nimble"}): logical value indicating whether derivatives should be prepared when preparing Nimble model (will esp. allow to use HMC sampler). Default to FALSE.
\item \code{resetMV} (only for \code{MCMC_language=="Nimble"}): logical value to be passed to $run specifying whether previous parameter samples should be reset or not. Default to FALSE to speed up WAIC calculations. You can turn it to TRUE if you wish to speed up runs of MCMC (cf. https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI/m/Su40lgNRBgAJ).
\item \code{parallelize} (only for \code{MCMC_language=="Jags"} and \code{MCMC_language=="Nimble"}): logical value specifying whether the MCMC should be parallelized within the \code{runMCMC_btadjust} function with the \code{parallel} package (and for the moment default settings of this package). Default to FALSE. If TRUE, library \code{parallel} should be loaded. If TRUE and \code{control$time.max} is unspecified or infinite, each parallelized process will have a maximum duration of 30 days. In case \code{MCMC_language=="Greta"}, parallelization is managed directly by Greta.
\item \code{parallelizeInitExpr} (only for \code{MCMC_language=="Jags"} and \code{MCMC_language=="Nimble"}): expression to add in each cluster created by parallelization. Default to \code{expression(if(MCMC_language=="Nimble"){library(nimble);if(control.MCMC$APT) {library(nimbleAPT)}} else {NULL})}.
\item \code{useConjugacy} (only for \code{MCMC_language=="Nimble"}): logical value specifying whether Nimble should search for conjugate priors in the model. Default to FALSE. If TRUE, can render model configuration shorter (https://groups.google.com/g/nimble-users/c/a6DFCefYfjU/m/kqUWx9UXCgAJ) at the expense of not allowing any conjugate sampler
\item \code{WAIC} (only for \code{MCMC_language=="Nimble"}): logical value specifying whether WAIC should be calculated online within Nimble. Default to FALSE.
\item \code{WAIC.Nsamples} (only for \code{MCMC_language=="Nimble"}): integer number: number of (nearly independent) samples of parameters in the posterior distribution over which to calculate WAIC. If WAIC is calculated on more than one chain, this number will be required for each of these chains. Default to 2000.
\item \code{WAIC.control} (only for \code{MCMC_language=="Nimble"}): named list or list of such named lists: list (or lists) specifying the control parameters to calculate WAIC online within Nimble. Default to list(online=TRUE,dataGroups=NULL,marginalizeNodes=NULL,niterMarginal=1000,convergenceSet=c(0.25,0.5,0.75),thin=TRUE,nburnin_extra=0). Given the way WAIC is here calculated (after convergence and over the last sample outputs for one chain by named list), components thin will be turned to TRUE and nburnin_extra to 0. If several lists are used, only at most the first Nchains - reduced by removal of fixed chains - lists will be taken into account to calculate WAICs in different ways on different Markov Chains.
\item \code{APT} (only for \code{MCMC_language=="Nimble"}): logical value specifying whether \code{NimbleAPT} should be used to allow Adaptive Parallel Tempering for at least a subset of the model parameters. If so, the sampler_RW_tempered sampler is declared for all parameters. If wishing to change sampler, use \code{confModel.expression.toadd}. Default to FALSE.
\item \code{APT.NTemps} (only for \code{MCMC_language=="Nimble"}): integer number: number of temperatures for \code{NimbleAPT}. Default to 7.
\item \code{APT.initTemps} (only for \code{MCMC_language=="Nimble"}): NULL or double vector of length \code{APT.NTemps}: initial temperatures for Nimble APT. Default to NULL in which case initial temperatures will be 1:APT.NTemps. The values should be increasing with a first value of 1.
\item \code{APT.tuneTemps} (only for \code{MCMC_language=="Nimble"}): numerical vector of length 2: values to feed the parameters \code{tuneTemper1} and \code{tuneTemper2} in \code{NimbleAPT}. See documentation of \code{NimbleAPT}. Default to c(10,0.7).
\item \code{APT.thinPrintTemps} (only for \code{MCMC_language=="Nimble"}): expression or numerical value : thinning parameter for printing temperatures in case APT. Default to expression(niter/5).
\item \code{includeAllStochNodes} (only for \code{MCMC_language=="Nimble"}): logical value specifying whether all stochastic nodes should be made available in \code{runMCMC_btadjust}. Default to FALSE. Maybe useful for the \code{extraCalculations} component of \code{control.MCMC}. Can include more parameters than \code{includeParentNodes}. Can be useful for example to calculate goodness-of-fit p-values.
\item \code{monitorAllStochNodes} (only for \code{MCMC_language=="Nimble"}): logical value specifying whether all stochastic nodes should be monitored for convergence and number of effective values in \code{runMCMC_btadjust}. Default to FALSE.
\item \code{saveAllStochNodes} (only for \code{MCMC_language=="Nimble"}): logical value specifying whether all stochastic nodes should be made available in \code{runMCMC_btadjust} and kept in the result of \code{runMCMC_btadjust}. Default to FALSE. Note that if TRUE, will change the content of \code{params_saved}.
\item \code{includeParentNodes} (only for \code{MCMC_language=="Nimble"}): logical value specifying whether parent stochastic nodes of data should be made available in \code{runMCMC_btadjust}. Default to FALSE. Maybe useful for the \code{extraCalculations} component of \code{control.MCMC} - for example for "offline" calculations of WAIC or DIC (Deviance information criterion).
\item \code{monitorParentNodes} (only for \code{MCMC_language=="Nimble"}): logical value specifying whether parent stochastic nodes of data should be monitored for convergence and number of effective values in \code{runMCMC_btadjust}. Default to FALSE.
\item \code{saveParentNodes} (only for \code{MCMC_language=="Nimble"}): logical value specifying whether parent stochastic nodes of data should be made available in \code{runMCMC_btadjust} and kept in the result of \code{runMCMC_btadjust}. Default to FALSE. Note that if TRUE, will change the content of \code{params_saved}.
\item \code{extraCalculations} (mainly useful for \code{MCMC_language=="Nimble"} but can be used with other languages as well): NULL value or expression that will be evaluated at the end of \code{runMCMC_btadjust}. The value of this expression will be saved in the \code{extraResults} component of the \code{final.params} component of the attributes of the result of \code{runMCMC_btadjust}. See the vignette devoted to the use of it. Default to NULL.
       }}
}
\value{
a \code{mcmc.list} object with attributes with the following components:
      \itemize{ \item \code{call.params}: a list containing most of the important arguments of the \code{runMCMC_btadjust} call as well as either a summary of dimensions/lengths and mean of components of \code{data} and \code{constants} arguments or their entire values.
      \item \code{final.params}: a list with the parameters of the MCMC at the end of fitting:
            \itemize{ \item \code{converged}: logical: TRUE if model has converged when stopping the MCMC, FALSE otherwise
            \item \code{neffs.reached}: logical: TRUE if model has converged and reached the objectives in terms of effective sample size, FALSE otherwise
            \item \code{final.Nchains}: number of Markov chains finally retained - since some chains could be excluded if invariable.
            \item \code{removed.chains}: identity of the Markov chains finally removed - those that would have been invariable during the first cycle.
            \item \code{burnin}: number of iterations of the transient (burn-in) period
            \item \code{thin}: number of iterations used for thinning the final output
            \item \code{niter.tot}: total number of iterations (of each MCMC chain)
            \item \code{Nvalues}: number of saved values (over all the MCMC chains)
            \item \code{neff.min}: minimum number of effective values in the final MCMC (over params.conv)
            \item \code{neff.median}: median number of effective values (over params.conv)
            \item \code{WAIC}: results of the calculus of online WAIC (only if control.MCMC$WAIC and \code{MCMC_language=="Nimble"}). One component per component of control.MCMC$control.WAIC. Each component has a WAIC component and then a WAICDetails component.
            \item \code{extraResults}: results of the implementation of control.MCMC$extraCalculations. Can be any kind of R object.
            \item \code{Temps}: results of the series of Temperatures met in APT (only if control.MCMC$APT and \code{MCMC_language=="Nimble"}). One line for each end of $run (=MCMC run).
            \item \code{duration}: total duration (elapsed time) of the fit (in seconds)
            \item \code{duration.MCMC.preparation}: duration (elapsed time) of MCMC preparation (in seconds)
            \item \code{duration.MCMC.transient}: duration (elapsed time) of the MCMC transient (burn-in) phase (in seconds)
            \item \code{duration.MCMC.asymptotic}: duration (elapsed time) of the MCMC asymptotic phase (in seconds)
            \item \code{duration.MCMC.after}: duration (elapsed time) of the MCMC phase after the asymptotic phase of sampling (e.g. to calculate WAIC) (in seconds)
            \item \code{duration.btadjust}: duration (elapsed time) outside MCMC preparation & fitting (in seconds)
            \item \code{CPUduration}: total CPU duration (user+system, self+child if not NA - otherwise only self) of the fit (in seconds)
            \item \code{CPUduration.MCMC.preparation}: CPU duration (user+system, self+child if not NA - otherwise only self) of MCMC preparation (in seconds)
            \item \code{CPUduration.MCMC.transient}: CPU duration (user+system, self+child if not NA - otherwise only self) of the MCMC transient (burn-in) phase (in seconds)
            \item \code{CPUduration.MCMC.asymptotic}: CPU duration (user+system, self+child if not NA - otherwise only self) of the MCMC asymptotic phase (in seconds)
            \item \code{CPUduration.MCMC.after}: CPU duration (user+system, self+child if not NA - otherwise only self) of the MCMC phase after the asymptotic phase of sampling (e.g. to calculate WAIC) (in seconds)
            \item \code{CPUduration.btadjust}: CPU duration (user+system, self+child if not NA - otherwise only self) outside MCMC preparation & fitting (in seconds)
            \item \code{childCPUduration}: total child CPU duration (user+system) of the fit (in seconds)
            \item \code{childCPUduration.MCMC.preparation}: child CPU duration (user+system) of MCMC preparation (in seconds)
            \item \code{childCPUduration.MCMC.transient}: child CPU duration (user+system) of the MCMC transient (burn-in) phase (in seconds)
            \item \code{childCPUduration.MCMC.asymptotic}: child CPU duration (user+system) of the MCMC asymptotic phase (in seconds)
            \item \code{childCPUduration.MCMC.after}: child CPU duration (user+system) of the MCMC phase after the asymptotic phase of sampling (e.g. to calculate WAIC) (in seconds)
            \item \code{childCPUduration.btadjust}: child CPU duration (user+system) outside MCMC preparation & fitting (in seconds)
            \item \code{time}: time (from Sys.time) at the end of model fitting
            }
       \item \code{final.diags}: a list with final diagnostics of the fit:
            \itemize{ \item \code{params}: parameters of the MCMC (burn-in, thin, niter...)
            \item \code{conv_synth}: synthetic output of convergence diagnostics
            \item \code{neff_synth}: synthetic output for calculations of effective sample sizes
            \item \code{conv}: raw convergence values for all the parameters being diagnosed if control$conveff.final.allparams is FALSE and all the parameters otherwise
            \item \code{neff}: raw effective sample size values for all the parameters being diagnosed if control$conveff.final.allparams is FALSE and all the parameters otherwise
            }
       \item \code{sessionInfo}: a list containing the result of the call to sessionInfo() function at the end of runMCMC_btadjust function; contains info on the platform, versions of packages, R version...;
       \item \code{warnings}: a list of the warning messages issued during fitting; unsure it still works with this version
       \item \code{error}: a list with the error messages issued during fitting; unsure it still works with this version
       }
}
\description{
returns a mcmc.list object which is the output of a Markov Chain Monte Carlo obtained after adjusting burn-in & thinning parameters to meet pre-specified criteria in terms of convergence & effective sample size - i.e. sample size adjusted for autocorrelation - of the MCMC output
}
\details{
Recap: \cr
If \code{MCMC_language=="Nimble"}, the code, data and constants arguments should be specified according to the requirements of \code{nimble} package. \cr
If \code{MCMC_language=="Jags"}, the code and  data arguments need to be specified as required by \code{rjags} package.  \cr
If \code{MCMC_language=="Greta"}, the model argument must be specified and should be the result of the \code{model} command in \code{greta} package.  \cr

Details on \code{check.convergence}: \cr
If FALSE, no check of convergence at all, after \code{nburnin.min} (& \code{recheck.convergence} is put to FALSE & \code{check.convergence.firstrun} is dominated by \code{check.convergence}). \cr
If TRUE, the convergence behavior is governed by \code{check.convergence.firstrun} & \code{recheck.convergence}.\cr

Example for \code{confModel.expression.toadd} component of \code{control.MCMC}:\cr
		\code{confModel.expression.toadd<-expression({ConfModel[[i]]$removeSamplers(c("alpha","dzetad","beta","exper_bias[2]","exper_bias[3]","exper_precision[2]","exper_precision[3]"))
			ConfModel[[i]]$addSampler(target = c("alpha","dzetad","beta"),type = "RW_block")
			ConfModel[[i]]$addSampler(target = c("exper_bias[2]","exper_bias[3]"),type = "RW_block")
			ConfModel[[i]]$addSampler(target = c("exper_precision[2]","exper_precision[3]"),type = "RW_block")
			})}

Remark for \code{params, params.conv, params.save}:\cr
in cases of parameters that are vectors, matrices... the \code{params} vector can contain only the name of the vector or matrix... in which case all its components will be used. It can also contain the names of individual components.
}
\examples{
# for examples with Nimble or Greta, see the Presentation Vignette.
\dontrun{
#generating data
set.seed(1)
y1000<-rnorm(n=1000,mean=600,sd=30)
ModelData <-list(mass = y1000,nobs = length(y1000))

#writing the Jags code as a character chain in R
modeltotransfer<-"model {

# Priors
population.mean ~ dunif(0,5000)
population.sd ~ dunif(0,100)

# Precision = 1/variance: Normal distribution parameterized by precision in Jags
population.variance <- population.sd * population.sd
precision <- 1 / population.variance

# Likelihood
for(i in 1:nobs){
  mass[i] ~ dnorm(population.mean, precision)
 }
 }"

#specifying the initial values
ModelInits <- function()
{list (population.mean = rnorm(1,600,90), population.sd = runif(1, 1, 30))}
params <- c("population.mean", "population.sd", "population.variance")
K<-3
set.seed(1)
Inits<-lapply(1:K,function(x){ModelInits()})

# running runMCMC_btadjust with MCMC_language="Jags":
set.seed(1)
out.mcmc.Coda<-runMCMC_btadjust(MCMC_language="Jags", code=modeltotransfer,
data=ModelData,
Nchains=K, params=params, inits=Inits,
niter.min=1000, niter.max=300000,
nburnin.min=100, nburnin.max=200000,
thin.min=1, thin.max=1000,
neff.min=1000, conv.max=1.05,
control=list(print.diagnostics=TRUE, neff.method="Coda"))

summary(out.mcmc.Coda)
}
}
