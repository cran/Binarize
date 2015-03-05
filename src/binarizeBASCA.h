/* 
 * File:   binarizationPredict.h
 * Author: Stefan Mundus
 *
 * Created on 2. Juli 2010, 09:37
 */

#ifndef _BINARIZEBASCA_H
#define	_BINARIZEBASCA_H

#ifdef	__cplusplus
extern "C" {
#endif
    #include "common.h"
    #if !DEBUG_MODE
    #include <Rinternals.h>
    #include <R.h>
    #endif
    

    void calc_First_Cost_Matrix_Line(dbl_matrix* result, dbl_array* vect);
    double mean_BASCA(double* values, int a, int b);
    double costs(double* values, int a, int b);
    void calc_RestCc_and_Ind_Matrices(dbl_matrix* cc, int_matrix* ind, dbl_array* vect);
    void calc_P_Matrix(int_matrix* P, int_matrix* Ind);
    double calc_jump_height(int_matrix* P, dbl_array* vect, int i, int j);
    double calc_error(int_matrix* P, dbl_array* vect, int i, int j);
    double calc_score(int_matrix* P, dbl_matrix* H, dbl_array* vect, int i, int j);
    void calc_V(int_array* v, dbl_array* Q_max, dbl_matrix* Q, dbl_matrix* H, int_matrix* P, dbl_array* vect);

    #if !DEBUG_MODE
    SEXP binarizeBASCA(SEXP vect, SEXP tau, SEXP numberofsamples);
    #endif
    //SEXP calc_p_value_R(SEXP vect, SEXP tau, SEXP numberofsamples);

    int alloc_Accelerator_Memory(int value_count);
    void free_Accelerator_Memory();



#ifdef	__cplusplus
}
#endif

#endif	/* _BINARIZEBASCA_H */

