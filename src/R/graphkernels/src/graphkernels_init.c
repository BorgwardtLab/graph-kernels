#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP graphkernels_CalculateGraphletKernelCpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP graphkernels_CalculateKernelCpp(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"graphkernels_CalculateGraphletKernelCpp", (DL_FUNC) &graphkernels_CalculateGraphletKernelCpp, 4},
    {"graphkernels_CalculateKernelCpp",         (DL_FUNC) &graphkernels_CalculateKernelCpp,         3},
    {NULL, NULL, 0}
};

void R_init_graphkernels(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
