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
 * This is the implementation file for the ARKODE CustomStep module.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"
#include "arkode_customstep_impl.h"
#include "arkode_interp_impl.h"
#include <sundials/sundials_math.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM ".32Lg"
#else
#define RSYM ".16g"
#endif

/* constants */
#define ZERO   RCONST(0.0)
#define ONE    RCONST(1.0)


/*===============================================================
  CustomStep Exported functions -- Required
  ===============================================================*/

/*---------------------------------------------------------------
  Create CustomStep integrator memory struct
  ---------------------------------------------------------------*/
void* CustomStepCreate(ARKRhsFn f,
                       realtype t0,
                       N_Vector y0)
{
  ARKodeMem           ark_mem;         /* outer ARKode memory   */
  ARKODECustomStepMem step_mem;        /* outer stepper memory  */
  booleantype         nvectorOK;
  int                 retval;

  /* Check that f is supplied */
  if (f == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::CustomStep",
                    "CustomStepCreate", MSG_ARK_NULL_F);
    return(NULL);
  }

  /* Check that y0 is supplied */
  if (y0 == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::CustomStep",
                    "CustomStepCreate", MSG_ARK_NULL_Y0);
    return(NULL);
  }

  /* Test if all required vector operations are implemented */
  nvectorOK = customStep_CheckNVector(y0);
  if (!nvectorOK) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::CustomStep",
                    "CustomStepCreate", MSG_ARK_BAD_NVECTOR);
    return(NULL);
  }

  /* Create ark_mem structure and set default values */
  ark_mem = arkCreate();
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::CustomStep",
                    "CustomStepCreate", MSG_ARK_NO_MEM);
    return(NULL);
  }

  /* Attach step_mem structure and function pointers to ark_mem */
  ark_mem->step_attachlinsol   = customStep_AttachLinsol;
  ark_mem->step_disablelsetup  = customStep_DisableLSetup;
  ark_mem->step_getlinmem      = customStep_GetLmem;
  ark_mem->step_getimplicitrhs = customStep_GetImplicitRHS;
  ark_mem->step_getgammas      = customStep_GetGammas;
  ark_mem->step_init           = customStep_Init;
  ark_mem->step_fullrhs        = customStep_FullRHS;
  ark_mem->step                = customStep_TakeStep;
  ark_mem->step_mem            = (void*) step_mem;

  /* Set default values for CustomStep optional inputs */


  /* Allocate the general stepper vectors using y0 as a template */
  /* NOTE: method-specific vectors will be allocated later on */

  /* Copy the RHS function into stepper memory */
  step_mem->f = f;

  /* Update the ARKode workspace requirements */
  ark_mem->liw += 2;  /* fcn/data ptr, int, long int, sunindextype, booleantype */
  ark_mem->lrw += 0;

  /* Set the linear solver addresses to NULL (we check != NULL later) */
  step_mem->linit  = NULL;
  step_mem->lsetup = NULL;
  step_mem->lsolve = NULL;
  step_mem->lfree  = NULL;
  step_mem->lmem   = NULL;

  /* Initialize all the counters */
  step_mem->nf        = 0;
  step_mem->nsetups   = 0;

  /* Initialize main ARKODE infrastructure (allocates vectors) */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKODE::CustomStep", "CustomStepCreate",
                    "Unable to initialize main ARKODE infrastructure");
    CustomStepFree((void**) &ark_mem);  return(NULL);
  }

  /* return ARKODE memory */
  return((void*) ark_mem);
}


/*---------------------------------------------------------------
  CustomStepResize:

  This routine resizes the memory within the CustomStep module.
  It first resizes the main ARKODE infrastructure memory, and
  then resizes its own data.
  ---------------------------------------------------------------*/
int CustomStepResize(void *arkode_mem,
                     N_Vector y0,
                     realtype t0,
                     ARKVecResizeFn resize,
                     void *resize_data)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  sunindextype lrw1, liw1, lrw_diff, liw_diff;
  int retval, i;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "CustomStepResize",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Determing change in vector sizes */
  lrw1 = liw1 = 0;
  if (y0->ops->nvspace != NULL)
    N_VSpace(y0, &lrw1, &liw1);
  lrw_diff = lrw1 - ark_mem->lrw1;
  liw_diff = liw1 - ark_mem->liw1;
  ark_mem->lrw1 = lrw1;
  ark_mem->liw1 = liw1;

  /* resize ARKODE infrastructure memory (use hscale = 1.0) */
  retval = arkResize(ark_mem, y0, RCONST(1.0), t0, resize, resize_data);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKODE::CustomStep", "CustomStepResize",
                    "Unable to resize main ARKODE infrastructure");
    return(retval);
  }

  /* Resize the internal vectors */
  /* if (step_mem->sdata != NULL) */
  /*   if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, */
  /*                     liw_diff, y0, &step_mem->sdata)) { */
  /*     arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE::CustomStep", */
  /*                     "CustomStepResize", "Unable to resize vector"); */
  /*     return(ARK_MEM_FAIL); */
  /*   } */

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  CustomStepReInit:

  This routine re-initializes the CustomStep module to solve a new
  problem of the same size as was previously solved (all counter
  values are set to 0).
  ---------------------------------------------------------------*/
int CustomStepReInit(void* arkode_mem,
                     ARKRhsFn f,
                     realtype t0,
                     N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "CustomStepReInit",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE::CustomStep",
                    "CustomStepReInit", MSG_ARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check that f is supplied */
  if (f == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::CustomStep",
                    "CustomStepReInit", MSG_ARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  /* Check that y0 is supplied */
  if (y0 == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::CustomStep",
                    "CustomStepReInit", MSG_ARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* ReInitialize main ARKODE infrastructure */
  retval = arkInit(arkode_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKODE::CustomStep", "CustomStepReInit",
                    "Unable to reinitialize main ARKODE infrastructure");
    return(retval);
  }

  /* Copy the input parameters into ARKODE state */
  step_mem->f = f;

  /* Initialize all the counters */
  step_mem->nf = 0;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  CustomStepReset:

  This routine resets the CustomStep module state to solve the same
  problem from the given time with the input state (all counter
  values are retained).
  ---------------------------------------------------------------*/
int CustomStepReset(void* arkode_mem,
                    realtype tR,
                    N_Vector yR)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "CustomStepReset",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, tR, yR, RESET_INIT);

  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKODE::CustomStep", "CustomStepReset",
                    "Unable to initialize main ARKODE infrastructure");
    return(retval);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  CustomStepSStolerances, CustomStepSVtolerances, CustomStepWFtolerances:

  These routines set integration tolerances (wrappers for general
  ARKODE utility routines)
  ---------------------------------------------------------------*/
int CustomStepSStolerances(void *arkode_mem,
                           realtype reltol,
                           realtype abstol)
{
  /* unpack ark_mem, call arkSStolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::CustomStep",
                    "CustomStepSStolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSStolerances(ark_mem, reltol, abstol));
}

int CustomStepSVtolerances(void *arkode_mem,
                           realtype reltol,
                           N_Vector abstol)
{
  /* unpack ark_mem, call arkSVtolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::CustomStep",
                    "CustomStepSVtolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSVtolerances(ark_mem, reltol, abstol));
}

int CustomStepWFtolerances(void *arkode_mem,
                           ARKEwtFn efun)
{
  /* unpack ark_mem, call arkWFtolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::CustomStep",
                    "CustomStepWFtolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkWFtolerances(ark_mem, efun));
}


/*---------------------------------------------------------------
  CustomStepRootInit:

  Initialize (attach) a rootfinding problem to the stepper
  (wrappers for general ARKODE utility routine)
  ---------------------------------------------------------------*/
int CustomStepRootInit(void *arkode_mem,
                       int nrtfn,
                       ARKRootFn g)
{
  /* unpack ark_mem, call arkRootInit, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::CustomStep",
                    "CustomStepRootInit", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkRootInit(ark_mem, nrtfn, g));
}


/*---------------------------------------------------------------
  CustomStepEvolve:

  This is the main time-integration driver (wrappers for general
  ARKODE utility routine)
  ---------------------------------------------------------------*/
int CustomStepEvolve(void *arkode_mem,
                     realtype tout,
                     N_Vector yout,
                     realtype *tret,
                     int itask)
{
  /* unpack ark_mem, call arkEvolve, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::CustomStep",
                    "CustomStepEvolve", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkEvolve(ark_mem, tout, yout, tret, itask));
}


/*---------------------------------------------------------------
  CustomStepGetDky:

  This returns interpolated output of the solution or its
  derivatives over the most-recently-computed step (wrapper for
  generic ARKODE utility routine)
  ---------------------------------------------------------------*/
int CustomStepGetDky(void *arkode_mem,
                     realtype t,
                     int k,
                     N_Vector dky)
{
  /* unpack ark_mem, call arkGetDky, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::CustomStep",
                    "CustomStepGetDky", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetDky(ark_mem, t, k, dky));
}


/*---------------------------------------------------------------
  CustomStepFree frees all CustomStep memory, and then calls an
  ARKODE utility routine to free the ARKODE infrastructure memory.
  ---------------------------------------------------------------*/
void CustomStepFree(void **arkode_mem)
{
  int j;
  sunindextype Cliw, Clrw;
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;

  /* nothing to do if arkode_mem is already NULL */
  if (*arkode_mem == NULL)  return;

  /* conditional frees on non-NULL CustomStep module */
  ark_mem = (ARKodeMem) (*arkode_mem);
  if (ark_mem->step_mem != NULL) {

    step_mem = (ARKODECustomStepMem) ark_mem->step_mem;

    /* free the linear solver memory */
    if (step_mem->lfree != NULL) {
      step_mem->lfree((void *) ark_mem);
      step_mem->lmem = NULL;
    }

    /* free the local vectors */
    /* if (step_mem->sdata != NULL) { */
    /*   arkFreeVec(ark_mem, &step_mem->sdata); */
    /*   step_mem->sdata = NULL; */
    /* } */

    /* free the time stepper module itself */
    free(ark_mem->step_mem);
    ark_mem->step_mem = NULL;
  }

  /* free memory for overall ARKODE infrastructure */
  arkFree(arkode_mem);
}


/*---------------------------------------------------------------
  CustomStepPrintMem:

  This routine outputs the memory from the CustomStep structure and
  the main ARKODE infrastructure to a specified file pointer
  (useful when debugging).
  ---------------------------------------------------------------*/
void CustomStepPrintMem(void* arkode_mem,
                        FILE* outfile)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int i, retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "CustomStepPrintMem",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return;

  /* if outfile==NULL, set it to stdout */
  if (outfile == NULL)  outfile = stdout;

  /* output data from main ARKODE infrastructure */
  fprintf(outfile,"CustomStep ARKODE infrastructure mem:\n");
  arkPrintMem(ark_mem, outfile);

  /* output stepper-specific quantities */

  /* output long integer quantities */
  fprintf(outfile,"CustomStep: nf = %li\n", step_mem->nf);
  fprintf(outfile,"CustomStep: nsetups = %li\n", step_mem->nsetups);

#ifdef SUNDIALS_DEBUG_PRINTVEC
  /* output vector quantities */
  /* fprintf(outfile, "CustomStep: sdata:\n"); */
  /* N_VPrintFile(step_mem->sdata, outfile); */
#endif

  return;
}



/*===============================================================
  CustomStep Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  Interface routines supplied to ARKODE
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  customStep_AttachLinsol:

  This routine attaches the various set of system linear solver
  interface routines, data structure, and solver type to the
  CustomStep module.
  ---------------------------------------------------------------*/
int customStep_AttachLinsol(void* arkode_mem,
                            ARKLinsolInitFn linit,
                            ARKLinsolSetupFn lsetup,
                            ARKLinsolSolveFn lsolve,
                            ARKLinsolFreeFn lfree,
                            SUNLinearSolver_Type lsolve_type,
                            void *lmem)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "customStep_AttachLinsol",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* free any existing system solver */
  if (step_mem->lfree != NULL)  step_mem->lfree(arkode_mem);

  /* Attach the provided routines, data structure and solve type */
  step_mem->linit  = linit;
  step_mem->lsetup = lsetup;
  step_mem->lsolve = lsolve;
  step_mem->lfree  = lfree;
  step_mem->lmem   = lmem;

  /* Reset all linear solver counters */
  step_mem->nsetups = 0;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  customStep_DisableLSetup:

  This routine NULLifies the lsetup function pointer in the
  CustomStep module.
  ---------------------------------------------------------------*/
void customStep_DisableLSetup(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "customStep_DisableLSetup",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return;

  /* nullify the lsetup function pointer */
  step_mem->lsetup = NULL;
}


/*---------------------------------------------------------------
  customStep_GetLmem:

  This routine returns the system linear solver interface memory
  structure, lmem.
  ---------------------------------------------------------------*/
void* customStep_GetLmem(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure, and return lmem */
  retval = customStep_AccessStepMem(arkode_mem, "customStep_GetLmem",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(NULL);
  return(step_mem->lmem);
}


/*---------------------------------------------------------------
  customStep_GetImplicitRHS:

  This routine returns the implicit RHS function pointer, f.
  ---------------------------------------------------------------*/
ARKRhsFn customStep_GetImplicitRHS(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure, and return fs */
  retval = customStep_AccessStepMem(arkode_mem, "customStep_GetImplicitRHS",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(NULL);
  return(step_mem->f);
}


/*---------------------------------------------------------------
  customStep_GetGammas:

  This routine fills the current value of gamma, and states
  whether the gamma ratio fails the dgmax criteria.
  ---------------------------------------------------------------*/
int customStep_GetGammas(void* arkode_mem,
                         realtype *gamma,
                         realtype *gamrat,
                         booleantype **jcur,
                         booleantype *dgamma_fail)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "customStep_GetGammas",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set outputs */
  *gamma  = step_mem->gamma;
  *gamrat = step_mem->gamrat;
  *jcur = &step_mem->jcur;
  *dgamma_fail = (SUNRabs(*gamrat - ONE) >= step_mem->dgmax);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  customStep_Init:

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup.

  With initialization types FIRST_INIT this routine:
  - sets/checks the method coefficients to be used
  - allocates any memory that depends on the number of
    stages, method order, or solver options
  - sets the call_fullrhs flag

  With other initialization types, this routine does nothing.
  ---------------------------------------------------------------*/
int customStep_Init(void* arkode_mem,
                    int init_type)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval, j;
  booleantype reset_efun;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "customStep_Init",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* immediately return if reset */
  if (init_type == RESET_INIT) return(ARK_SUCCESS);

  /* initializations/checks for (re-)initialization call */
  if (init_type == FIRST_INIT) {

    /* enforce fixed outer step size */
    if (!ark_mem->fixedstep) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::CustomStep", "customStep_Init",
                      "Adaptive time stepping is not currently supported");
      return(ARK_ILL_INPUT);
    }

    /* allocate sdata, ... */
    /* if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->sdata))) */
    /*   return(ARK_MEM_FAIL); */

    /* Limit interpolant degree based on method order (use negative
       argument to specify update instead of overwrite) */
    if (ark_mem->interp != NULL) {
      retval = arkInterpSetDegree(ark_mem, ark_mem->interp, -(step_mem->q-1));
      if (retval != ARK_SUCCESS) {
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::CustomStep", "customStep_Init",
                        "Unable to update interpolation polynomial degree");
        return(ARK_ILL_INPUT);
      }
    }

  }

  /* Call linit (if it exists) */
  if (step_mem->linit) {
    retval = step_mem->linit(ark_mem);
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_LINIT_FAIL, "ARKODE::CustomStep", "customStep_Init",
                      MSG_ARK_LINIT_FAIL);
      return(ARK_LINIT_FAIL);
    }
  }

  /* Signal to shared arkode module that fullrhs is required after each step */
  ark_mem->call_fullrhs = SUNTRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  customStep_FullRHS:

  This is just a wrapper to call the user-supplied RHS function,
  f(t,y).

  This will be called in one of three 'modes':
    ARK_FULLRHS_START -> called at the beginning of a simulation
                         or after post processing at step
    ARK_FULLRHS_END   -> called at the end of a successful step
    ARK_FULLRHS_OTHER -> called elsewhere (e.g. for dense output)

  In this template, no shortcuts are applied.
  ---------------------------------------------------------------*/
int customStep_FullRHS(void* arkode_mem,
                       realtype t,
                       N_Vector y,
                       N_Vector f,
                       int mode)
{
  ARKodeMem ark_mem;
  ARKODECustomStepMem step_mem;
  int retval;

  /* access ARKODECustomStepMem structure */
  retval = customStep_AccessStepMem(arkode_mem, "customStep_FullRHS",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  retval = step_mem->f(t, y, f, ark_mem->user_data);
  step_mem->nf++;
  if (retval != 0) {
    arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::CustomStep",
                    "customStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
    return(ARK_RHSFUNC_FAIL);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  customStep_TakeStep:

  This routine serves the primary purpose of the CustomStep module:
  it performs a single time step (with embedding, if possible).

  The output variable dsmPtr should contain estimate of the
  weighted local error if an embedding is present; otherwise it
  should be 0.

  The input/output variable nflagPtr is used to gauge convergence
  of any algebraic solvers within the step.  At the start of a new
  time step, this will initially have the value FIRST_CALL.  On
  return from this function, nflagPtr should have a value:
            0 => algebraic solve completed successfully
           >0 => solve did not converge at this step size
                 (but may with a smaller stepsize)
           <0 => solve encountered an unrecoverable failure

  The return value from this routine is:
            0 => step completed successfully
           >0 => step encountered recoverable failure;
                 reduce step and retry (if possible)
           <0 => step encountered unrecoverable failure
  ---------------------------------------------------------------*/
int customStep_TakeStep(void* arkode_mem,
                        realtype *dsmPtr,
                        int *nflagPtr)
{
  ARKodeMem ark_mem; /* ARKODE infrastructure memory */
  ARKODECustomStepMem step_mem; /* custom stepper memory */
  int is;     /* current stage index        */
  int retval; /* reusable return flag       */

  /* initialize algebraic solver convergence flag to success;
     error estimate to zero */
  *nflagPtr = ARK_SUCCESS;
  *dsmPtr = ZERO;

  /* access the CustomStep mem structure */
  retval = customStep_AccessStepMem(arkode_mem, "customStep_TakeStep",
                                    &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* various entities from ark_mem that may be of use include:

     ark_mem->nst: current time step index
     ark_mem->h: current time step size
     ark_mem->tn: time value at start of step
     ark_mem->tcur: current time value (start of step)
         this should be updated as the step proceeds
     ark_mem->ycur: current state (start of step)
         this should be updated as the step proceeds
     ark_mem->user_data: data structure for calls to user-supplied fcns
     ark_mem->tempv1: temporary vector
     ark_mem->tempv2: temporary vector
     ark_mem->tempv3: temporary vector
     ark_mem->tempv4: temporary vector

   */


  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  Internal utility routines
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  customStep_AccessStepMem:

  Shortcut routine to unpack ark_mem and step_mem structures from
  void* pointer.  If either is missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int customStep_AccessStepMem(void* arkode_mem,
                             const char *fname,
                             ARKodeMem *ark_mem,
                             ARKODECustomStepMem *step_mem)
{

  /* access ARKodeMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::CustomStep",
                    fname, MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *ark_mem = (ARKodeMem) arkode_mem;
  if ((*ark_mem)->step_mem==NULL) {
    arkProcessError(*ark_mem, ARK_MEM_NULL, "ARKODE::CustomStep",
                    fname, MSG_MRISTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *step_mem = (ARKODECustomStepMem) (*ark_mem)->step_mem;
  return(ARK_SUCCESS);
}



/*---------------------------------------------------------------
  customStep_CheckNVector:

  This routine checks if all required vector operations are
  present.  If any of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
booleantype customStep_CheckNVector(N_Vector tmpl)
{
  if ( (tmpl->ops->nvclone     == NULL) ||
       (tmpl->ops->nvdestroy   == NULL) ||
       (tmpl->ops->nvlinearsum == NULL) ||
       (tmpl->ops->nvconst     == NULL) ||
       (tmpl->ops->nvscale     == NULL) ||
       (tmpl->ops->nvwrmsnorm  == NULL) )
    return(SUNFALSE);
  return(SUNTRUE);
}


/*===============================================================
  EOF
  ===============================================================*/
