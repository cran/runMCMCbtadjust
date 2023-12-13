# runMCMCbtadjust 1.0.5

- Correction of a bug in the section 2.1 section of the code for the calculation of the new number of iterations: the number of iterations was too small when thin.max was reached (acknowledged with thin.max=1): replaced thin by thin.theoretical (which has been introduced at the beginning or section 2.1); replaced (index.conv-1) by ifelse(converged,index.conv-1,size.samplesList) in final calculations of durations. Otherwise, erroneous in case of non convergence..
- Made corrections regarding how control$seed was treated: (i) added the possibility that control$seed could be NULL - in which case there is no control of seeds for MCMCs; (ii) put the default to NULL instead of 1; (iii) in case length(control$seed) was 1, added set.seed(NULL) at the end to not constrain the seed of the R environment after the end of the function by the seed used in the function.
- Added neffs.reached in the final.params component of the attributes of the result; added help for final.params$converged and final.params$neffs.reached
- Allowed control.MCMC$n.adapt to also apply in case Nimble was the MCMC.language; changed the default value accordingly (so that it is in the program 1000 in case Jags is used and 0 in case Nimble is used); corrected an error in the inclusion of n.adapt in duration.MCMC.transient in case Jags is used.

---

# runMCMCbtadjust 1.0.4

- Replaced (index.conv-1) by ifelse(converged,index.conv-1,size.samplesList) in final calculations of durations. Otherwise, erroneous in case of non convergence..

---

# runMCMCbtadjust 1.0.3

- Added buildDerivs component to control.MCMC parameter.

---

# runMCMCbtadjust 1.0.2

- Changed default value for round.thinmult in parameter control to TRUE, which is more rigorous.

---

# runMCMCbtadjust 1.0.1

- Added a `NEWS.md` file to track changes to the package.
- Changed title of package - referring to 'JAGS', 'nimble' or 'greta' explicitly
- Changed function runMCMC_btadjust(): placing neff.max at the beginning, replacing neffmax by neff.max in functions scale.available.neffs & calculate.thinmult.target, adding warnings if niter.max is finite and niter.max/thin.max is greater than 3*neff.max
- Replaced niter>0 by niter>=thin in sections 2 and 2.2 in function runMCMC_btadjust(), to avoid a bug
- Slight changes to the vignette named 'runMCMCbtadjust_Presentation'

---
# runMCMCbtadjust 1.0.0
