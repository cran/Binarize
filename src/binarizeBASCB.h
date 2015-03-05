/* 
 * File:   binarizationScalespace.h
 * Author: Stefan Mundus
 *
 * Created on 5. Juli 2010, 16:42
 */

#ifndef _BINARIZEBASCB_H
#define	_BINARIZEBASCB_H

#ifdef	__cplusplus
extern "C" {
#endif
    #include "common.h"
    #if !DEBUG_MODE
    #include <Rinternals.h>
    #include <R.h>
    #endif
    

    extern dbl_matrix* b;
    extern int_matrix* b_returned;

    typedef struct mgs_result{
        dbl_matrix* smoothed;
        int_matrix* zerocrossing;
        dbl_array* deriv;
    }
    mgs_result;

    typedef struct quant_result{
        int_matrix* steps;
        int greatest_steps_row;
        int greatest_steps_col;
        int_array* index;
        int greatest_index_ind;
    }
    quant_result;

    typedef struct calc_V_result{
        int_array* v;
        dbl_matrix* smoothedX;
        dbl_matrix* meanlist;
    }
    calc_V_result;

    double cost_Scalespace(dbl_array* vect, int a, int b, double y);
    #if SAVE_BESSEL_VALUES && !DEBUG_MODE
    void save_bessel_values();
    #endif
    #if DEBUG_MODE
    void restore_bessel_values();
    #endif
    void calc_V_Scalespace(calc_V_result* result, mgs_result* mgs_res, quant_result* q_res, dbl_matrix* H, dbl_array* vect);
    void revert_int_matrix(int_matrix* mat);
    void getQuantizations(quant_result* qr, mgs_result* mr);
    void mgs(mgs_result* result, dbl_array* vect, dbl_array* sigma);
    #if !DEBUG_MODE
    SEXP binarizeBASCB(SEXP vect, SEXP tau, SEXP numberofSamples, SEXP sigma);
    #endif


#ifdef	__cplusplus
}
#endif

#endif	/* _BINARIZEBASCB_H */

