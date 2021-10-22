/* -----------------------------------------------------------------
 * Programmer(s): John C. Doe @ State University
 * -----------------------------------------------------------------
 * ARKODE-stepper-template Copyright Start
 * Copyright (c) 2021, State University.
 * All rights reserved.
 *
 * See the top-level LICENSE file for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * ARKODE-stepper-template Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the ARKODE CustomStep module.
 * -----------------------------------------------------------------*/

#ifndef _CUSTOMSTEP_H
#define _CUSTOMSTEP_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_nonlinearsolver.h>
#include <arkode/arkode.h>
#include <arkode/arkode_ls.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* --------------------------------------------------------------
   Public CustomStep Constants

   Insert any "#define NAME VALUE" constants here as shortcuts
   for users
   -------------------------------------------------------------- */





/* --------------------------------------------------------------
   CustomStep User-Supplied Function Types

   Insert any "typedef type (*NameFn)(type a, type b);" function
   types that users should supply for this stepper (beyond those
   already in arkode.h)
   -------------------------------------------------------------- */




/* ---------------------------------------------------------------
   Custom data structures and associated utility routines

   If the custom stepper needs to define any user-facing data
   structures (e.g., that are used to store method coefficients)
   then those data structures and corresponding utility routines
   can be inserted here (or in a separate header file that is
   #include-d above)
   --------------------------------------------------------------- */




/* ---------------------------------------------------------------
   Exported Functions

   Functions that define the user interface for this stepper
   belong below.  Some of these with our "standard" naming
   conventions are provided.
   --------------------------------------------------------------- */

/* Create, Resize, Reinitialization and Reset functions */
SUNDIALS_EXPORT void* CustomStepCreate(ARKRhsFn f,
                                       realtype t0,
                                       N_Vector y0);

SUNDIALS_EXPORT int CustomStepResize(void *arkode_mem,
                                     N_Vector ynew,
                                     realtype t0,
                                     ARKVecResizeFn resize,
                                     void *resize_data);

SUNDIALS_EXPORT int CustomStepReInit(void* arkode_mem,
                                     ARKRhsFn f,
                                     realtype t0,
                                     N_Vector y0);

SUNDIALS_EXPORT int CustomStepReset(void* arkode_mem,
                                    realtype tR,
                                    N_Vector yR);

/* Tolerance input functions */
SUNDIALS_EXPORT int CustomStepSStolerances(void *arkode_mem,
                                        realtype reltol,
                                        realtype abstol);
SUNDIALS_EXPORT int CustomStepSVtolerances(void *arkode_mem,
                                        realtype reltol,
                                        N_Vector abstol);
SUNDIALS_EXPORT int CustomStepWFtolerances(void *arkode_mem,
                                        ARKEwtFn efun);

/* Linear solver set function */
SUNDIALS_EXPORT int CustomStepSetLinearSolver(void *arkode_mem,
                                           SUNLinearSolver LS,
                                           SUNMatrix A);

/* Rootfinding initialization */
SUNDIALS_EXPORT int CustomStepRootInit(void *arkode_mem,
                                       int nrtfn,
                                       ARKRootFn g);


/* Optional input functions -- must be called AFTER CustomStepCreate

   These are typically implemented in arkode_customstep_io.c, many
   of which just wrap corresponding functions from the underlying
   ARKODE infrastructure */

SUNDIALS_EXPORT int CustomStepSetInterpolantType(void *arkode_mem,
                                                 int itype);
SUNDIALS_EXPORT int CustomStepSetInterpolantDegree(void *arkode_mem,
                                                   int degree);
SUNDIALS_EXPORT int CustomStepSetNonlinearSolver(void *arkode_mem,
                                                 SUNNonlinearSolver NLS);
SUNDIALS_EXPORT int CustomStepSetNlsRhsFn(void *arkode_mem,
                                          ARKRhsFn nls_fs);
SUNDIALS_EXPORT int CustomStepSetMaxNumSteps(void *arkode_mem,
                                             long int mxsteps);
SUNDIALS_EXPORT int CustomStepSetLSetupFrequency(void *arkode_mem,
                                                 int msbp);
SUNDIALS_EXPORT int CustomStepSetPredictorMethod(void *arkode_mem,
                                                 int method);
SUNDIALS_EXPORT int CustomStepSetMaxNonlinIters(void *arkode_mem,
                                                int maxcor);
SUNDIALS_EXPORT int CustomStepSetNonlinConvCoef(void *arkode_mem,
                                                realtype nlscoef);
SUNDIALS_EXPORT int CustomStepSetMaxHnilWarns(void *arkode_mem,
                                              int mxhnil);
SUNDIALS_EXPORT int CustomStepSetStopTime(void *arkode_mem,
                                          realtype tstop);
SUNDIALS_EXPORT int CustomStepSetFixedStep(void *arkode_mem,
                                           realtype hsfixed);
SUNDIALS_EXPORT int CustomStepSetRootDirection(void *arkode_mem,
                                               int *rootdir);
SUNDIALS_EXPORT int CustomStepSetNoInactiveRootWarn(void *arkode_mem);
SUNDIALS_EXPORT int CustomStepSetErrHandlerFn(void *arkode_mem,
                                              ARKErrHandlerFn ehfun,
                                              void *eh_data);
SUNDIALS_EXPORT int CustomStepSetErrFile(void *arkode_mem,
                                         FILE *errfp);
SUNDIALS_EXPORT int CustomStepSetUserData(void *arkode_mem,
                                          void *user_data);
SUNDIALS_EXPORT int CustomStepSetDiagnostics(void *arkode_mem,
                                             FILE *diagfp);

/* Optional linear solver interface input functions -- must be called
   AFTER CustomStepSetLinearSolver

   These are typically implemented in arkode_customstep_io.c */

SUNDIALS_EXPORT int CustomStepSetJacFn(void *arkode_mem, ARKLsJacFn jac);
SUNDIALS_EXPORT int CustomStepSetJacEvalFrequency(void *arkode_mem,
                                                  long int msbj);
SUNDIALS_EXPORT int CustomStepSetEpsLin(void *arkode_mem,
                                        realtype eplifac);
SUNDIALS_EXPORT int CustomStepSetLSNormFactor(void *arkode_mem,
                                              realtype nrmfac);
SUNDIALS_EXPORT int CustomStepSetPreconditioner(void *arkode_mem,
                                                ARKLsPrecSetupFn psetup,
                                                ARKLsPrecSolveFn psolve);
SUNDIALS_EXPORT int CustomStepSetJacTimes(void *arkode_mem,
                                          ARKLsJacTimesSetupFn jtsetup,
                                          ARKLsJacTimesVecFn jtimes);
SUNDIALS_EXPORT int CustomStepSetJacTimesRhsFn(void *arkode_mem,
                                               ARKRhsFn jtimesRhsFn);
SUNDIALS_EXPORT int CustomStepSetLinSysFn(void *arkode_mem,
                                          ARKLsLinSysFn linsys);


/* User-facing routines to be called repeatedly throughout
   integration of the problem */

/* Main time-integration routine -- integrates the IVP over an interval in t */
SUNDIALS_EXPORT int CustomStepEvolve(void *arkode_mem,
                                     realtype tout,
                                     N_Vector yout,
                                     realtype *tret,
                                     int itask);

/* Compute the kth derivative of the y function at time t */
SUNDIALS_EXPORT int CustomStepGetDky(void *arkode_mem,
                                     realtype t,
                                     int k,
                                     N_Vector dky);

/* Utility function to retrieve the current internal solution */
SUNDIALS_EXPORT int CustomStepGetCurrentState(void *arkode_mem,
                                              N_Vector *state);


/* Optional output functions

   These are typically implemented in arkode_customstep_io.c */

SUNDIALS_EXPORT int CustomStepGetNumRhsEvals(void *arkode_mem,
                                             long int *nf_evals);
SUNDIALS_EXPORT int CustomStepGetNumLinSolvSetups(void *arkode_mem,
                                                  long int *nlinsetups);
SUNDIALS_EXPORT int CustomStepGetWorkSpace(void *arkode_mem,
                                           long int *lenrw,
                                           long int *leniw);
SUNDIALS_EXPORT int CustomStepGetNumSteps(void *arkode_mem,
                                          long int *nssteps);
SUNDIALS_EXPORT int CustomStepGetLastStep(void *arkode_mem,
                                          realtype *hlast);
SUNDIALS_EXPORT int CustomStepGetCurrentTime(void *arkode_mem,
                                             realtype *tcur);
SUNDIALS_EXPORT int CustomStepGetTolScaleFactor(void *arkode_mem,
                                                realtype *tolsfac);
SUNDIALS_EXPORT int CustomStepGetErrWeights(void *arkode_mem,
                                            N_Vector eweight);
SUNDIALS_EXPORT int CustomStepGetNumGEvals(void *arkode_mem,
                                           long int *ngevals);
SUNDIALS_EXPORT int CustomStepGetRootInfo(void *arkode_mem,
                                          int *rootsfound);
SUNDIALS_EXPORT char *CustomStepGetReturnFlagName(long int flag);
SUNDIALS_EXPORT int CustomStepWriteParameters(void *arkode_mem,
                                              FILE *fp);

/* Nonlinear solver optional output functions */
SUNDIALS_EXPORT int CustomStepGetNumNonlinSolvIters(void *arkode_mem,
                                                    long int *nniters);
SUNDIALS_EXPORT int CustomStepGetNumNonlinSolvConvFails(void *arkode_mem,
                                                        long int *nncfails);
SUNDIALS_EXPORT int CustomStepGetNonlinSolvStats(void *arkode_mem,
                                                 long int *nniters,
                                                 long int *nncfails);

/* Linear solver optional output functions */
SUNDIALS_EXPORT int CustomStepGetLinWorkSpace(void *arkode_mem,
                                              long int *lenrwLS,
                                              long int *leniwLS);
SUNDIALS_EXPORT int CustomStepGetNumJacEvals(void *arkode_mem,
                                             long int *njevals);
SUNDIALS_EXPORT int CustomStepGetNumPrecEvals(void *arkode_mem,
                                              long int *npevals);
SUNDIALS_EXPORT int CustomStepGetNumPrecSolves(void *arkode_mem,
                                               long int *npsolves);
SUNDIALS_EXPORT int CustomStepGetNumLinIters(void *arkode_mem,
                                             long int *nliters);
SUNDIALS_EXPORT int CustomStepGetNumLinConvFails(void *arkode_mem,
                                                 long int *nlcfails);
SUNDIALS_EXPORT int CustomStepGetNumJTSetupEvals(void *arkode_mem,
                                                 long int *njtsetups);
SUNDIALS_EXPORT int CustomStepGetNumJtimesEvals(void *arkode_mem,
                                                long int *njvevals);
SUNDIALS_EXPORT int CustomStepGetNumLinRhsEvals(void *arkode_mem,
                                                long int *nfevalsLS);
SUNDIALS_EXPORT int CustomStepGetLastLinFlag(void *arkode_mem,
                                             long int *flag);

SUNDIALS_EXPORT char *CustomStepGetLinReturnFlagName(long int flag);


/* Free function */
SUNDIALS_EXPORT void CustomStepFree(void **arkode_mem);

/* Output the CustomStep memory structure (useful when debugging) */
SUNDIALS_EXPORT void CustomStepPrintMem(void* arkode_mem,
                                        FILE* outfile);


#ifdef __cplusplus
}
#endif

#endif
