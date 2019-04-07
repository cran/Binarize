#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP binarizeBASCA(SEXP, SEXP, SEXP);
extern SEXP binarizeBASCB(SEXP, SEXP, SEXP, SEXP);
extern SEXP TASCA(SEXP, SEXP, SEXP);
extern SEXP TASCA_min(SEXP, SEXP, SEXP);
extern SEXP TASCB(SEXP, SEXP, SEXP, SEXP);
extern SEXP TASCB_min(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"binarizeBASCA", (DL_FUNC) &binarizeBASCA, 3},
    {"binarizeBASCB", (DL_FUNC) &binarizeBASCB, 4},
    {"TASCA",         (DL_FUNC) &TASCA,         3},
    {"TASCA_min",     (DL_FUNC) &TASCA_min,     3},
    {"TASCB",         (DL_FUNC) &TASCB,         4},
    {"TASCB_min",     (DL_FUNC) &TASCB_min,     4},
    {NULL, NULL, 0}
};

void R_init_Binarize(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

