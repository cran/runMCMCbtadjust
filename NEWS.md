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
