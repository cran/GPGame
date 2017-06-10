#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP GPGame_getPoffs(SEXP, SEXP, SEXP, SEXP);
extern SEXP GPGame_getPoffsCross(SEXP, SEXP, SEXP, SEXP);
extern SEXP GPGame_PSNE_sparseMat(SEXP, SEXP, SEXP);
extern SEXP GPGame_PSNE_sparseMat_cross(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GPGame_PSNE_sparseMat_sorted(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"GPGame_getPoffs",              (DL_FUNC) &GPGame_getPoffs,              4},
    {"GPGame_getPoffsCross",         (DL_FUNC) &GPGame_getPoffsCross,         4},
    {"GPGame_PSNE_sparseMat",        (DL_FUNC) &GPGame_PSNE_sparseMat,        3},
    {"GPGame_PSNE_sparseMat_cross",  (DL_FUNC) &GPGame_PSNE_sparseMat_cross,  5},
    {"GPGame_PSNE_sparseMat_sorted", (DL_FUNC) &GPGame_PSNE_sparseMat_sorted, 3},
    {NULL, NULL, 0}
};

void R_init_GPGame(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
