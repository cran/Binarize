/* 
 * File:   binarizationScanStatistic.h
 * Author: stefan
 *
 * Created on 12. Juli 2010, 17:36
 */

#ifndef _BINARIZESCANSTATISTIC_H
#define	_BINARIZESCANSTATISTIC_H

#ifdef	__cplusplus
extern "C" {
#endif
    #include "common.h"
    #if !DEBUG_MODE
    #include <Rinternals.h>
    #include <R.h>
    #endif

    typedef struct scan_result{
        int_array* binarized_vector;
        double* threshold;
        int* reject;
        double* p_value;
    }
    scan_result;

    typedef struct nk_stack_element{
        int n;
        int k;
        struct nk_stack_element* next;
    }
    nk_stack_element;

    double n_over_k(int n, int k);
    double Gb(int k, int N, double w);
    double b_val(int k, int N, double w);
    double probability(int k, int N, double w);
    void calc_probability(dbl_matrix* solution, dbl_array* vect_sorted, double sign_level, double windowsize, double w);
    void solutionListAdapt(dbl_matrix* list, double start, double stop, double p, int start_ind, int stop_ind);
    int getSolutionListLength(dbl_matrix* list);
    void scanStatistic(scan_result* result, dbl_array* vect, double w, double sign_level);
    #if !DEBUG_MODE
    SEXP binarizeScanStatistic(SEXP vect, SEXP w, SEXP sign_level);
    #endif



#ifdef	__cplusplus
}
#endif

#endif	/* _BINARIZESCANSTATISTIC_H */

