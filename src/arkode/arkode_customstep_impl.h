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
 * Implementation header file for the ARKODE CustomStep module.
 * ---------------------------------------------------------------------------*/

#ifndef _ARKODE_CUSTOMSTEP_IMPL_H
#define _ARKODE_CUSTOMSTEP_IMPL_H

#include "arkode/arkode_customstep.h"
#include "arkode_impl.h"
#include "arkode_ls_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/* --------------------------------------------------------------
   Private CustomStep Constants

   Insert any "#define NAME VALUE" constants here as
   implementation shortcuts
   -------------------------------------------------------------- */


/*===============================================================
  Custom time step module data structure
  ===============================================================*/

/*---------------------------------------------------------------
  The type ARKODECustomStepMem is type pointer to struct
  ARKODECustomStepMemRec. This structure contains fields to
  perform a Custom time step.
  ---------------------------------------------------------------*/
typedef struct ARKODECustomStepMemRec {

  /* problem specification */
  ARKRhsFn f;    /* y' = f(t,y) */

  /* method data and parameters */

  /* algebraic solver data and parameters */

  /* linear solver interface function/data pointers */
  ARKLinsolInitFn    linit;
  ARKLinsolSetupFn   lsetup;
  ARKLinsolSolveFn   lsolve;
  ARKLinsolFreeFn    lfree;
  void              *lmem;

  /* work counters */
  long int nf;        /* num f calls                   */
  long int nsetups;   /* num linear solver setup calls */

} *ARKODECustomStepMem;


/*===============================================================
  Custom time step module private function prototypes
  ===============================================================*/

/* Interface routines supplied to ARKODE infrastructure */
int customStep_AttachLinsol(void* arkode_mem,
                            ARKLinsolInitFn linit,
                            ARKLinsolSetupFn lsetup,
                            ARKLinsolSolveFn lsolve,
                            ARKLinsolFreeFn lfree,
                            SUNLinearSolver_Type lsolve_type,
                            void *lmem);
void customStep_DisableLSetup(void* arkode_mem);
int customStep_Init(void* arkode_mem,
                    int init_type);
void* customStep_GetLmem(void* arkode_mem);
ARKRhsFn customStep_GetImplicitRHS(void* arkode_mem);
int customStep_GetGammas(void* arkode_mem,
                         realtype *gamma,
                         realtype *gamrat,
                         booleantype **jcur,
                         booleantype *dgamma_fail);
int customStep_FullRHS(void* arkode_mem,
                       realtype t,
                       N_Vector y,
                       N_Vector f,
                       int mode);
int customStep_TakeStep(void* arkode_mem,
                        realtype *dsmPtr,
                        int *nflagPtr);

/* Internal utility routines */
int customStep_AccessStepMem(void* arkode_mem,
                             const char *fname,
                             ARKodeMem *ark_mem,
                             ARKODECustomStepMem *step_mem);
booleantype customStep_CheckNVector(N_Vector tmpl);


/*===============================================================
  Reusable CustomStep Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_CUSTOMSTEP_NO_MEM    "Time step module memory is NULL."

#ifdef __cplusplus
}
#endif

#endif
