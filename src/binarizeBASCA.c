/*
 * File:   main.c
 * Author: Stefan Mundus
 *
 * Created on 30. Juni 2010, 11:28
 */

#include "binarizeBASCA.h"
#if DEBUG_MODE
#include "mersennetwister.h"
#else
#include <Rmath.h>
#endif
#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <time.h>

//the computed costs will be stored in c (so they don't have to be computed again)
dbl_matrix* c = 0;
//in returned is stored how often a certain cost-value was returned
//if the c_returned[a][b] == 0 the cost-value haven't been computed so far
int_matrix* c_returned = 0;

//sames as storing the compted costs, here for the mean values
dbl_matrix* m = 0;
int_matrix* m_returned = 0;

//sames as storing the compted costs, here for the approximation error values
dbl_array* e = 0;
int_array* e_returned = 0;

//helper function(s) and interfaces ---------------------------------------------------------------------------------------------------------


//main functions --------------------------------------------------------------------------------------------------------------------------

/*
 * calculates the values for the first line of the Cc Matrix
 */
void calc_First_Cost_Matrix_Line(dbl_matrix* result, dbl_array* vect)
{
    int i,j;
    double dist;
    //holds the current mean-value
    double mean_val = 0;

    for(i = 0; i < vect->length; i++)
    {
        //update the current mean value with a new value
        mean_val += ((vect->values[vect->length - i - 1] - mean_val) / ((double)i + 1.0));
        if(m && m_returned)
        {
            m->values[vect->length - i - 1][vect->length - 1] = mean_val;
            m_returned->values[vect->length - i - 1][vect->length - 1]++;
        }

        if(i > 0)
        {
            //calculate the quadratic distance of the last i values to the current mean value
            for(j = vect->length - i - 1; j < vect->length; j++)
            {
                dist = vect->values[j] - mean_val;
                result->values[0][vect->length - i - 1] += dist * dist;
            }

            //store the calculated costs in c
            if(c && c_returned)
            {
                c->values[vect->length - i - 1][vect->length - 1] = result->values[0][vect->length - i - 1];
                c_returned->values[vect->length - i - 1][vect->length - 1]++;
            }
        }
        else
        {
            //store the calculated costs in c
            if(c && c_returned)
            {
                c->values[vect->length-1][vect->length-1] = 0.0;
                c_returned->values[vect->length-1][vect->length-1]++;
            }
        }
    }
}

/*
 * Calculates the mean values in the interval [a,b]
 */
double mean_BASCA(double* values, int a, int b)
{
    int i;
    double mean_val;

    //check if already computed
    if(m && m_returned && m_returned->values[a][b])
    {
        m_returned->values[a][b]++;
        return m->values[a][b];
    }

    //compute the mean
    mean_val = 0.0;
    for(i = a; i <= b; i++){
        mean_val += values[i];
    }
    mean_val /= ((double)b-(double)a+1.0);

    //store the computed value
    if(m && m_returned)
    {
        m_returned->values[a][b]++;
        m->values[a][b] = mean_val;
    }

    return mean_val;
}

/*
 * calcutes the costs(quadratic error to the mean) in the interval [a,b]
 */
double costs(double* values, int a, int b)
{
    //variable definitions
    double mean_val,costs,costs_root;
    int i;

    //check if the costs for parameter a,b have been previously computed.
    //if that's true give back the old result and increment the respective "returned"-value by one
    if(c && c_returned && c_returned->values[a][b])
    {
        c_returned->values[a][b]++;
        return c->values[a][b];
    }

    //calculate the mean value for the interval [a,b]
    mean_val = mean_BASCA(values, a, b);

    //calculate the costs for interval [a,b] (which means the quadratic error of
    //the function values in this interval to their mean value
    costs = 0.0;
    for(i = a; i <= b; i++){
        costs_root = values[i] - mean_val;
        costs += costs_root * costs_root;
    }

    //give back the result and increment the respective "returned"-value by one
    if(c && c_returned)
    {
        c->values[a][b] = costs;
        c_returned->values[a][b]++;
    }

    return costs;
}

/*
 * calculates the values for all lines of the Cc Matrix, except of the first line,
 * and all the values of the Ind-Matrix
 */
void calc_RestCc_and_Ind_Matrices(dbl_matrix* cc, int_matrix* ind, dbl_array* vect)
{
    //variable definitions
    int i,j,d,min_costs_ind;
    double min_costs,cur_costs;
    int cols = cc->cols;

    //Implemention of Algorithm 1 of the paper (except line 1 which is done in function "calc_First_Cost_Matrix_Line")
    for(j = 0; j < vect->length - 2; j++)
    {
        for(i = 0; i < vect->length - j - 1; i++)
        {
            //find the minimal costs and the corresponding index for inserting j steps starting at index i
            min_costs = DBL_MAX;
            min_costs_ind = -1;
            for(d = i; d < vect->length - j - 1; d++)
            {
                //the query is required, because of the fact I omitted the last column of the Cc matrix, which is always 0
                //And one time during this algorithm this 0 is needed
                cur_costs = d + 1 < cols ? costs(vect->values, i, d) + cc->values[j][d+1]/*[MAT2ARR_CC(j - 1, d + 1, cols)]*/ : costs(vect->values, i, d);
                //save the new minimal value and index
                if(cur_costs < min_costs)
                {
                    min_costs = cur_costs;
                    min_costs_ind = d;
                }
            }

            //store the computed minimum and the corresponding index
            cc->values[j+1][i]/*[MAT2ARR_CC(j, i, cols)]*/ = min_costs;
            ind->values[j][i]/*[MAT2ARR(j, i, cols)]*/ = min_costs_ind + 1;
        }
    }
}

/*
 * Calculates the P-Matrix, where the indices for the optimal step-functions for each number
 * of steps are stored.
 */
void calc_P_Matrix(int_matrix* P, int_matrix* Ind)
{
    int i,j,z;

    //Implemention of Algorithm 2 from the paper
    for(j = 0; j < P->rows; j++)
    {
        z = j;
        P->values[j][0] = Ind->values[z--][0];
        if(j > 0)
        {
            for(i = 1; i <= j; i++)
            {
                P->values[j][i] = Ind->values[z--][P->values[j][i-1]];
            }
        }
    }
}

/*
 * Calculates the jump height for the step at index P_i(j)
 */
double calc_jump_height(int_matrix* P, dbl_array* vect, int i, int j)
{
    double height;

    //compute the jump height according to the definition from the paper
    if(i == 0 && j == 0)
    {
        height = mean_BASCA(vect->values, P->values[j][i], vect->length - 1);
        height -= mean_BASCA(vect->values, 0, P->values[j][i]-1);
    }
    else if(i == 0 && j > 0)
    {
        height = mean_BASCA(vect->values, P->values[j][i], P->values[j][i+1]-1);
        height -= mean_BASCA(vect->values, 0, P->values[j][i]-1);
    }
    else if(i == j && i > 0)
    {
        height = mean_BASCA(vect->values, P->values[j][i], vect->length - 1);
        height -= mean_BASCA(vect->values, P->values[j][i-1], P->values[j][i]-1);
    }
    else
    {
        height = mean_BASCA(vect->values, P->values[j][i], P->values[j][i+1]-1);
        height -= mean_BASCA(vect->values, P->values[j][i-1], P->values[j][i]-1);
    }

    return height;
}

/*
 * Calculates the approximation error for the step at index P_i(j)
 */
double calc_error(int_matrix* P, dbl_array* vect, int i, int j)
{
    double z,err,e_root;
    int d, index;

    //check if value was computed before
    index = P->values[j][i]-1;
    if(e && e_returned && e_returned->values[index])
    {
        e_returned->values[index]++;
        return e->values[index];
    }

    //compute the approximation error according to the definition from the paper
    z = (vect->values[P->values[j][i]-1] + vect->values[P->values[j][i]]) * 0.5;

    err = 0.0;
    for(d = 0; d < vect->length; d++)
    {
        e_root = vect->values[d] - z;
        err += e_root * e_root;
    }

    //store computed value
    if(e && e_returned)
    {
        e_returned->values[index]++;
        e->values[index] = err;
    }

    return err;
}

/*
 * Calculates the score for a certain step
 */
double calc_score(int_matrix* P, dbl_matrix* H, dbl_array* vect, int i, int j)
{
    double score;

    score = calc_jump_height(P, vect, i, j);
    H->values[j][i] = score;
    score /= calc_error(P, vect, i, j);

    return score;
}

/*
 * Calculates the Vector V which contains the indices of the strongest discontinuities for each amount of steps
 */
void calc_V(int_array* v, dbl_array* Q_max, dbl_matrix* Q, dbl_matrix* H, int_matrix* P, dbl_array* vect)
{
    int i,j,max_ind;
    double max_val,cur_val;

    for(j = 0; j < P->rows; j++)
    {
        max_val = -1.0;
        max_ind = -1;
        for(i = 0; i <= j; i++)
        {
            cur_val = calc_score(P, H, vect, i, j);
            Q->values[j][i] = cur_val;
            if(cur_val > max_val)
            {
                max_val = cur_val;
                max_ind = i;
            }
        }
        v->values[j] = P->values[j][max_ind];
        Q_max->values[j] = max_val;
    }
}

/*
 * allocates memory for storing all the computed intermediate values, that could be reused
 */
int alloc_Accelerator_Memory(int value_count)
{
    int bytes = 0;
    
    c = init_dbl_matrix(0, value_count, value_count, 0);
    bytes += value_count * value_count * sizeof(double);

    c_returned = init_int_matrix(0, value_count, value_count, 0);
    bytes += value_count * value_count * sizeof(int);
    
    m = init_dbl_matrix(0, value_count, value_count, 0);
    bytes += value_count * value_count * sizeof(double);
    m_returned = init_int_matrix(0, value_count, value_count, 0);
    bytes += value_count * value_count * sizeof(int);
    
    e = init_dbl_array(0, value_count - 1, 0);
    bytes += value_count * sizeof(double);
    e_returned = init_int_array(0, value_count - 1, 0);
    bytes += value_count * sizeof(int);

    return bytes;
}

/*
 * Free the memory used for accelerating computations. This method also sets all pointers to 0.
 */
void free_Accelerator_Memory(void)
{
    destroy_dbl_matrix(c);
    destroy_int_matrix(c_returned);
    destroy_dbl_matrix(m);
    destroy_int_matrix(m_returned);
    destroy_dbl_array(e);
    destroy_int_array(e_returned);
    c_returned = 0;
    c = 0;
    m_returned = 0;
    m = 0;
    e_returned = 0;
    e = 0;
}

#if !DEBUG_MODE
/*
 * Interface Function for calling direct from R. This Method manages the conversions from R to C objects and from C to R objects,
 * and calls the necessary functions for computing the binarized vector.
 */
SEXP binarizeBASCA(SEXP vect, SEXP tau, SEXP numberofsamples)
{
    //get the lengths
    //int i,j,sum,sum_tot;
    //int bytes = 0;
    int value_count = length(vect);
    int vc_m1 = value_count - 1;
    int vc_m2 = vc_m1 - 1;
    dbl_array *vector, *vect_sorted, *Q_Max_vec;
    dbl_matrix *Cc_Mat, *Q_Mat, *H_Mat;
    int_matrix *Ind_Mat, *P_Mat;
    int_array *v_vec;
    final_result f_res;

    //sort the vect into vect_sorted
    vector = init_dbl_array(REAL(vect), value_count, 1);

    vect_sorted = init_dbl_array(0, value_count, 0);
    memcpy(vect_sorted->values, vector->values, vect_sorted->length * sizeof(double));
    qsort(vect_sorted->values, vect_sorted->length, sizeof(double), comp);

    //name the required SEXP Objects
    SEXP result, binarized_vector, threshold, p_value, other_results, Cc, Ind, P, Q, H, Q_max, v, Names;

    //allocate memory for saving calculated values
    alloc_Accelerator_Memory(value_count);
    
    //allocate the final result and set the names of the entries
    PROTECT(result = allocVector(VECSXP, 4));
    PROTECT(Names = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(Names,0, mkChar("binarized_vector"));
    SET_VECTOR_ELT(Names,1, mkChar("threshold"));
    SET_VECTOR_ELT(Names,2, mkChar("p_value"));
    SET_VECTOR_ELT(Names,3, mkChar("other_results"));
    setAttrib(result, R_NamesSymbol, Names);
    UNPROTECT(1);

    PROTECT(other_results = allocVector(VECSXP, 7));
    PROTECT(Names = allocVector(VECSXP, 7));
    SET_VECTOR_ELT(Names,0, mkChar("Cc"));
    SET_VECTOR_ELT(Names,1, mkChar("Ind"));
    SET_VECTOR_ELT(Names,2, mkChar("P_Mat"));
    SET_VECTOR_ELT(Names,3, mkChar("Q_Mat"));
    SET_VECTOR_ELT(Names,4, mkChar("H_Mat"));
    SET_VECTOR_ELT(Names,5, mkChar("maximal_Qs"));
    SET_VECTOR_ELT(Names,6, mkChar("v_vec"));
    setAttrib(other_results, R_NamesSymbol, Names);
    UNPROTECT(1);


    //allocate memory for result matrices and vectors and set the matrix values to zero (because they aren't
    //all overwritten by the functions)
    PROTECT(binarized_vector = allocVector(INTSXP, value_count));
    PROTECT(threshold = allocVector(REALSXP, 1));
    PROTECT(p_value = allocVector(REALSXP, 1));
    PROTECT(Cc = allocMatrix(REALSXP, vc_m1, vc_m1));
    PROTECT(Ind = allocMatrix(INTSXP, vc_m1, vc_m2));
    PROTECT(P = allocMatrix(INTSXP, vc_m2, vc_m2));
    PROTECT(Q = allocMatrix(REALSXP, vc_m2, vc_m2));
    PROTECT(H = allocMatrix(REALSXP, vc_m2, vc_m2));
    PROTECT(Q_max = allocVector(REALSXP, vc_m2));
    PROTECT(v = allocVector(INTSXP, vc_m2));

    Cc_Mat = init_dbl_matrix(REAL(Cc), vc_m1, vc_m1, 0);
    Ind_Mat = init_int_matrix(INTEGER(Ind), vc_m2, vc_m1, 0);
    P_Mat = init_int_matrix(INTEGER(P), vc_m2, vc_m2, 0);
    v_vec = init_int_array(INTEGER(v), vc_m2, 0);
    Q_Max_vec = init_dbl_array(REAL(Q_max), vc_m2, 0);
    Q_Mat = init_dbl_matrix(REAL(Q), vc_m2, vc_m2, 0);
    H_Mat = init_dbl_matrix(REAL(H), vc_m2, vc_m2, 0);
    f_res.binarized_vector = init_int_array(INTEGER(binarized_vector), value_count, 0);
    f_res.p = REAL(p_value);
    f_res.threshold = REAL(threshold);

    //start the computation of the entries of all matrices
    calc_First_Cost_Matrix_Line(Cc_Mat, vect_sorted);
    calc_RestCc_and_Ind_Matrices(Cc_Mat, Ind_Mat, vect_sorted);
    calc_P_Matrix(P_Mat, Ind_Mat);
    calc_V(v_vec, Q_Max_vec, Q_Mat, H_Mat, P_Mat, vect_sorted);
    
    //free the memory for calculated values
    free_Accelerator_Memory();

    //calculate the final three results
    calc_final_results(&f_res, v_vec, vector, vect_sorted, *REAL(tau), *INTEGER(numberofsamples));

    //free(vect_sorted);
    destroy_dbl_array(vector);
    destroy_dbl_array(vect_sorted);
    destroy_dbl_matrix(Cc_Mat);
    destroy_int_matrix(Ind_Mat);
    destroy_int_matrix(P_Mat);
    destroy_int_array(v_vec);
    destroy_dbl_array(Q_Max_vec);
    destroy_dbl_matrix(Q_Mat);
    destroy_dbl_matrix(H_Mat);
    destroy_int_array(f_res.binarized_vector);

    //assign the computed elements to the final result
    SET_VECTOR_ELT(other_results,0, Cc);
    SET_VECTOR_ELT(other_results,1, Ind);
    SET_VECTOR_ELT(other_results,2, P);
    SET_VECTOR_ELT(other_results,3, Q);
    SET_VECTOR_ELT(other_results,4, H);
    SET_VECTOR_ELT(other_results,5, Q_max);
    SET_VECTOR_ELT(other_results,6, v);

    SET_VECTOR_ELT(result,0, binarized_vector);
    SET_VECTOR_ELT(result,1, threshold);
    SET_VECTOR_ELT(result,2, p_value);
    SET_VECTOR_ELT(result,3, other_results);

    UNPROTECT(12);
    return result;
}
#endif
