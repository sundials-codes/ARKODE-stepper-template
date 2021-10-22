/* -----------------------------------------------------------------------------
 * Programmer(s): John C. Doe @ State University
 * -----------------------------------------------------------------------------
 * ARKODE-stepper-template Copyright Start
 * Copyright (c) 2021, State University.
 * All rights reserved.
 *
 * See the top-level LICENSE file for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * ARKODE-stepper-template Copyright End
 * -----------------------------------------------------------------------------
 * This is the implementation file for the optional input and output functions
 * for the ARKODE CustomStep time stepper module.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_customstep_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM "Lg"
#else
#define RSYM "g"
#endif


/*===============================================================
  CustomStep Optional input functions (wrappers for generic ARKODE
  utility routines).  All are documented in arkode_io.c.
  ===============================================================*/
int CustomStepSetInterpolantType(void *arkode_mem,
                                 int itype) {
  return(arkSetInterpolantType(arkode_mem, itype)); }
int CustomStepSetInterpolantDegree(void *arkode_mem,
                                   int degree) {
  if (degree < 0) degree = ARK_INTERP_MAX_DEGREE;
  return(arkSetInterpolantDegree(arkode_mem, degree)); }
int CustomStepSetErrHandlerFn(void *arkode_mem,
                              ARKErrHandlerFn ehfun,
                              void *eh_data) {
  return(arkSetErrHandlerFn(arkode_mem, ehfun, eh_data)); }
int CustomStepSetErrFile(void *arkode_mem,
                         FILE *errfp) {
  return(arkSetErrFile(arkode_mem, errfp)); }
int CustomStepSetDiagnostics(void *arkode_mem,
                             FILE *diagfp) {
  return(arkSetDiagnostics(arkode_mem, diagfp)); }
int CustomStepSetMaxNumSteps(void *arkode_mem,
                             long int mxsteps) {
  return(arkSetMaxNumSteps(arkode_mem, mxsteps)); }
int CustomStepSetMaxHnilWarns(void *arkode_mem,
                              int mxhnil) {
  return(arkSetMaxHnilWarns(arkode_mem, mxhnil)); }
int CustomStepSetStopTime(void *arkode_mem,
                          realtype tstop) {
  return(arkSetStopTime(arkode_mem, tstop)); }
int CustomStepSetRootDirection(void *arkode_mem,
                               int *rootdir) {
  return(arkSetRootDirection(arkode_mem, rootdir)); }
int CustomStepSetNoInactiveRootWarn(void *arkode_mem) {
  return(arkSetNoInactiveRootWarn(arkode_mem)); }


/*---------------------------------------------------------------
  These wrappers for ARKLs module 'set' routines all are
  documented in arkode_customstep.h.
  ---------------------------------------------------------------*/
int CustomStepSetLinearSolver(void *arkode_mem,
                              SUNLinearSolver LS,
                              SUNMatrix A) {
  return(arkLSSetLinearSolver(arkode_mem, LS, A)); }
int CustomStepSetJacFn(void *arkode_mem,
                       ARKLsJacFn jac) {
  return(arkLSSetJacFn(arkode_mem, jac)); }
int CustomStepSetJacEvalFrequency(void *arkode_mem,
                                  long int msbj) {
  return(arkLSSetJacEvalFrequency(arkode_mem, msbj)); }
int CustomStepSetEpsLin(void *arkode_mem,
                        realtype eplifac) {
  return(arkLSSetEpsLin(arkode_mem, eplifac)); }
int CustomStepSetLSNormFactor(void *arkode_mem,
                              realtype nrmfac) {
  return(arkLSSetNormFactor(arkode_mem, nrmfac)); }
int CustomStepSetPreconditioner(void *arkode_mem,
                                ARKLsPrecSetupFn psetup,
                                ARKLsPrecSolveFn psolve) {
  return(arkLSSetPreconditioner(arkode_mem, psetup, psolve)); }
int CustomStepSetJacTimes(void *arkode_mem,
                          ARKLsJacTimesSetupFn jtsetup,
                          ARKLsJacTimesVecFn jtimes) {
  return(arkLSSetJacTimes(arkode_mem, jtsetup, jtimes)); }
int CustomStepSetJacTimesRhsFn(void *arkode_mem,
                               ARKRhsFn jtimesRhsFn) {
  return(arkLSSetJacTimesRhsFn(arkode_mem, jtimesRhsFn)); }
int CustomStepSetLinSysFn(void *arkode_mem,
                          ARKLsLinSysFn linsys) {
  return(arkLSSetLinSysFn(arkode_mem, linsys)); }


/*===============================================================
  CustomStep Optional output functions (wrappers for generic ARKODE
  utility routines).  All are documented in arkode_io.c.
  ===============================================================*/
int CustomStepGetNumSteps(void *arkode_mem,
                          long int *nssteps) {
  return(arkGetNumSteps(arkode_mem, nssteps)); }
int CustomStepGetLastStep(void *arkode_mem,
                          realtype *hlast) {
  return(arkGetLastStep(arkode_mem, hlast)); }
int CustomStepGetCurrentTime(void *arkode_mem,
                             realtype *tcur) {
  return(arkGetCurrentTime(arkode_mem, tcur)); }
int CustomStepGetCurrentState(void *arkode_mem,
                              N_Vector *state) {
  return(arkGetCurrentState(arkode_mem, state)); }
int CustomStepGetTolScaleFactor(void *arkode_mem,
                                realtype *tolsfact) {
  return(arkGetTolScaleFactor(arkode_mem, tolsfact)); }
int CustomStepGetErrWeights(void *arkode_mem,
                            N_Vector eweight) {
  return(arkGetErrWeights(arkode_mem, eweight)); }
int CustomStepGetWorkSpace(void *arkode_mem,
                           long int *lenrw,
                           long int *leniw) {
  return(arkGetWorkSpace(arkode_mem, lenrw, leniw)); }
int CustomStepGetNumGEvals(void *arkode_mem,
                           long int *ngevals) {
  return(arkGetNumGEvals(arkode_mem, ngevals)); }
int CustomStepGetRootInfo(void *arkode_mem,
                          int *rootsfound) {
  return(arkGetRootInfo(arkode_mem, rootsfound)); }
char *CustomStepGetReturnFlagName(long int flag) {
  return(arkGetReturnFlagName(flag)); }

/*---------------------------------------------------------------
  These wrappers for ARKLs module 'get' routines all are
  documented in arkode_customstep.h.
  ---------------------------------------------------------------*/
int CustomStepGetLinWorkSpace(void *arkode_mem,
                              long int *lenrwLS,
                              long int *leniwLS) {
  return(arkLSGetWorkSpace(arkode_mem, lenrwLS, leniwLS)); }
int CustomStepGetNumJacEvals(void *arkode_mem,
                             long int *njevals) {
  return(arkLSGetNumJacEvals(arkode_mem, njevals)); }
int CustomStepGetNumPrecEvals(void *arkode_mem,
                              long int *npevals) {
  return(arkLSGetNumPrecEvals(arkode_mem, npevals)); }
int CustomStepGetNumPrecSolves(void *arkode_mem,
                               long int *npsolves) {
  return(arkLSGetNumPrecSolves(arkode_mem, npsolves)); }
int CustomStepGetNumLinIters(void *arkode_mem,
                             long int *nliters) {
  return(arkLSGetNumLinIters(arkode_mem, nliters)); }
int CustomStepGetNumLinConvFails(void *arkode_mem,
                                 long int *nlcfails) {
  return(arkLSGetNumConvFails(arkode_mem, nlcfails)); }
int CustomStepGetNumJTSetupEvals(void *arkode_mem,
                                 long int *njtsetups) {
  return(arkLSGetNumJTSetupEvals(arkode_mem, njtsetups)); }
int CustomStepGetNumJtimesEvals(void *arkode_mem,
                                long int *njvevals) {
  return(arkLSGetNumJtimesEvals(arkode_mem, njvevals)); }
int CustomStepGetNumLinRhsEvals(void *arkode_mem,
                                long int *nfevalsLS) {
  return(arkLSGetNumRhsEvals(arkode_mem, nfevalsLS)); }
int CustomStepGetLastLinFlag(void *arkode_mem,
                             long int *flag) {
  return(arkLSGetLastFlag(arkode_mem, flag)); }
char *CustomStepGetLinReturnFlagName(long int flag) {
  return(arkLSGetReturnFlagName(flag)); }



/*===============================================================
  CustomStep optional input functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  CustomStepSetUserData:

  Wrapper for generic arkSetUserData and arkLSSetUserData
  routines.
  ---------------------------------------------------------------*/
int CustomStepSetUserData(void *arkode_mem, void *user_data)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "CustomStepSetUserData",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* set user_data in ARKode mem */
  retval = arkSetUserData(arkode_mem, user_data);
  if (retval != ARK_SUCCESS) return(retval);

  /* set user data in ARKodeLS mem */
  if (step_mem->lmem != NULL) {
    retval = arkLSSetUserData(arkode_mem, user_data);
    if (retval != ARKLS_SUCCESS) return(retval);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  CustomStepSetFixedStep:

  Wrapper for generic arkSetFixedStep routine.  Additionally
  enforces current CustomStep constraint for fixed time-stepping.
  ---------------------------------------------------------------*/
int CustomStepSetFixedStep(void *arkode_mem, realtype hsfixed)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::CustomStep",
                    "CustomStepSetFixedStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (hsfixed == ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::CustomStep",
                    "CustomStepSetFixedStep",
                    "CustomStep does not support adaptive steps at this time.");
    return(ARK_ILL_INPUT);
  }

  /* call generic routine for remaining work */
  return(arkSetFixedStep(ark_mem, hsfixed));
}


/*---------------------------------------------------------------
  CustomStepSetLSetupFrequency:

  Specifies the user-provided linear setup decision constant
  msbp.  Positive values give the frequency for calling lsetup;
  negative values imply recomputation of lsetup at each nonlinear
  solve; a zero value implies a reset to the default.
  ---------------------------------------------------------------*/
int CustomStepSetLSetupFrequency(void *arkode_mem, int msbp)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "CustomStepSetLSetupFrequency",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if argument legal set it, otherwise set default */
  if (msbp == 0) {
    step_mem->msbp = MSBP;
  } else {
    step_mem->msbp = msbp;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  CustomStepSetPredictorMethod:

  Specifies the method to use for predicting implicit solutions.
  Non-default choices are {1,2,3,4}, all others will use default
  (trivial) predictor.
  ---------------------------------------------------------------*/
int CustomStepSetPredictorMethod(void *arkode_mem, int pred_method)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "CustomStepSetPredictorMethod",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set parameter */
  step_mem->predictor = pred_method;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  CustomStepSetMaxNonlinIters:

  Specifies the maximum number of nonlinear iterations during
  one solve.  A non-positive input implies a reset to the
  default value.
  ---------------------------------------------------------------*/
int CustomStepSetMaxNonlinIters(void *arkode_mem, int maxcor)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "CustomStepSetMaxNonlinIters",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* argument <= 0 sets default, otherwise set input */
  if (maxcor <= 0) {
    step_mem->maxcor = MAXCOR;
  } else {
    step_mem->maxcor = maxcor;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  CustomStepSetNonlinConvCoef:

  Specifies the coefficient in the nonlinear solver convergence
  test.  A non-positive input implies a reset to the default value.
  ---------------------------------------------------------------*/
int CustomStepSetNonlinConvCoef(void *arkode_mem, realtype nlscoef)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "CustomStepSetNonlinConvCoef",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* argument <= 0 sets default, otherwise set input */
  if (nlscoef <= ZERO) {
    step_mem->nlscoef = NLSCOEF;
  } else {
    step_mem->nlscoef = nlscoef;
  }

  return(ARK_SUCCESS);
}


/*===============================================================
  CustomStep optional output functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  CustomStepGetNumRhsEvals:

  Returns the current number of calls to f
  ---------------------------------------------------------------*/
int CustomStepGetNumRhsEvals(void *arkode_mem, long int *nf_evals)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "CustomStepGetNumRhsEvals",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* get number of f evals from step_mem */
  *nfs_evals = step_mem->nf;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  CustomStepGetNumLinSolvSetups:

  Returns the current number of calls to the lsetup routine
  ---------------------------------------------------------------*/
int CustomStepGetNumLinSolvSetups(void *arkode_mem, long int *nlinsetups)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "CustomStepGetNumLinSolvSetups",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* get value from step_mem */
  *nlinsetups = step_mem->nsetups;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  CustomStepGetNumNonlinSolvIters:

  Returns the current number of nonlinear solver iterations
  ---------------------------------------------------------------*/
int CustomStepGetNumNonlinSolvIters(void *arkode_mem, long int *nniters)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "CustomStepGetNumNonlinSolvIters",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  *nniters = step_mem->nls_iters;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  CustomStepGetNumNonlinSolvConvFails:

  Returns the current number of nonlinear solver convergence fails
  ---------------------------------------------------------------*/
int CustomStepGetNumNonlinSolvConvFails(void *arkode_mem, long int *nncfails)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "CustomStepGetNumNonlinSolvConvFails",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set output from step_mem */
  *nncfails = ark_mem->ncfn;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  CustomStepGetNonlinSolvStats:

  Returns nonlinear solver statistics
  ---------------------------------------------------------------*/
int CustomStepGetNonlinSolvStats(void *arkode_mem, long int *nniters,
                              long int *nncfails)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "CustomStepGetNonlinSolvStats",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  *nniters  = step_mem->nls_iters;
  *nncfails = ark_mem->ncfn;

  return(ARK_SUCCESS);
}


/*===============================================================
  CustomStep parameter output
  ===============================================================*/

/*---------------------------------------------------------------
  CustomStepWriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int CustomStepWriteParameters(void *arkode_mem, FILE *fp)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "CustomStepWriteParameters",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* output ARKODE infrastructure parameters first */
  retval = arkWriteParameters(arkode_mem, fp);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE::CustomStep",
                    "CustomStepWriteParameters",
                    "Error writing ARKODE infrastructure parameters");
    return(retval);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
