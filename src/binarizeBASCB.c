#include "common.h"
#if DEBUG_MODE
#include "mersennetwister.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#else
#include <Rmath.h>
#endif
#include "binarizeBASCB.h"

dbl_matrix* b = 0;
int_matrix* b_returned = 0;

int enableWarnings(int val)
{
  SEXP s, t;
  PROTECT(t = s = allocList(2));
  SET_TYPEOF(s, LANGSXP);
  SETCAR(t, install("options")); 
  t = CDR(t);
  SETCAR(t,allocVector(INTSXP, 1));
  INTEGER(CAR(t))[0] = val;
  SET_TAG(t, install("warn"));
  SEXP oldStatus;
  PROTECT(oldStatus = coerceVector(eval(s, R_GlobalEnv),INTSXP));
  UNPROTECT(2);
  return INTEGER(oldStatus)[0];
} 

/*
 * Inverts the order of the written rows of a int_matrix. It works only for matrices which has written lines
 * at the beginning and the unwritten lines at the end.
 */
void revert_int_matrix(int_matrix* mat)
{
    int i, j, numbytes, *temp_row;

    numbytes = mat->cols * sizeof(int);
    temp_row = malloc(numbytes);

    for(j=0, i = mat->rows - 1; i >= 0 && i > j; i--)
    {
        if(mat->values[i][0])
        {
            memcpy(temp_row, mat->values[i], numbytes);
            memcpy(mat->values[i], mat->values[j], numbytes);
            memcpy(mat->values[j], temp_row, numbytes);
            j++;
        }
    }

    free(temp_row);
}

void mgs(mgs_result* result, dbl_array* vect, dbl_array* sigma)
{
    int i, j, k, is_zerocrossing, **zerocrossing;
    double deriv_tot, s_times_two, e_pow_mstt, sum, **smoothed, bessel, *deriv;

    int oldStatus = enableWarnings(-1);

    smoothed = result->smoothed->values;
    deriv = result->deriv->values;
    zerocrossing = result->zerocrossing->values;

    for(i = 0, deriv_tot = 0.0; i < result->deriv->length; i++)
    {
        deriv[i] = vect->values[i+1] - vect->values[i];
        deriv_tot += deriv[i];
    }

    for(i = 0; i < sigma->length; i++)
    {
        s_times_two = 2.0 * sigma->values[i];
        e_pow_mstt = exp(-s_times_two);
        for(k = 0; k < result->deriv->length; k++)
        {
            for(j = 0, sum = 0.0; j < result->deriv->length; j++)
            {
                if(b && b_returned)
                {
                    if(b_returned->values[i][abs(k-j)])
                    {
                        bessel = b->values[i][abs(k-j)];
                        b_returned->values[i][abs(k-j)]++;
                    }
                    else
                    {
                        bessel = bessel_i(s_times_two, (double)(k-j), 1.0);
                        b->values[i][abs(k-j)] = bessel;
                        b_returned->values[i][abs(k-j)]++;
                    }
                }
                else
                {
                    bessel = bessel_i(s_times_two, (double)(k-j), 1.0);
                }

                sum += deriv[j] * e_pow_mstt * bessel;
            }
            smoothed[i][k] = sum / deriv_tot;
        }

        for(j = 0, k = 0; j < result->smoothed->cols; j++)
        {
            is_zerocrossing = 0;
            is_zerocrossing += (int)(j == 0 && smoothed[i][j] > smoothed[i][j+1]);
            is_zerocrossing += (int)(j == result->smoothed->cols - 1 && smoothed[i][j-1] < smoothed[i][j]);
            is_zerocrossing += (int)(j > 0 && j < result->smoothed->cols - 1 && smoothed[i][j-1] < smoothed[i][j] && smoothed[i][j] > smoothed[i][j+1]);
            if(is_zerocrossing)
                zerocrossing[i][k++] = j + 1;
        }
        
        if (k == 0)
        // fix for linear functions (no maxima in slope):
        // all positions are considered as maxima
        {
          for(j = 0; j < result->zerocrossing->cols; j++)
          {
            zerocrossing[i][j] = j + 1;
          }
        }

    }
    enableWarnings(oldStatus);
}

void getQuantizations(quant_result* qr, mgs_result* mr)
{
    int cross, current, numbytes, qr_pos, i;

    //print_int_matrix(mr->zerocrossing, "zerocrossing");

    numbytes = mr->zerocrossing->cols * sizeof(int);

    for(cross = 0, current = -1, qr_pos = 0; cross < mr->zerocrossing->rows; cross++)
    {
        if(cross == 0 || memcmp(mr->zerocrossing->values[current], mr->zerocrossing->values[cross], numbytes))
        {
            memcpy(qr->steps->values[qr_pos], mr->zerocrossing->values[cross], numbytes);
            qr->index->values[qr_pos++] = cross + 1;
            current = cross;
            
            //break if there's only 1 value in current row
            if(!mr->zerocrossing->values[current][1])
                cross = mr->zerocrossing->rows;
        }
    }

    qr->greatest_steps_row = qr_pos;
    qr->greatest_index_ind = qr_pos;

    qsort(qr->index->values, qr->index->length, sizeof(int), comp_desc_int);
    revert_int_matrix(qr->steps);

    //print_int_matrix(qr->steps, "steps");

    for(i = 0; i <= mr->zerocrossing->cols; i++)
    {
        if(i == mr->zerocrossing->cols || !mr->zerocrossing->values[0][i])
        {
            qr->greatest_steps_col = i;
            i = mr->zerocrossing->cols + 1;
        }
    }
}

void calc_V_Scalespace(calc_V_result* result, mgs_result* mgs_res, quant_result* q_res, dbl_matrix* H, dbl_array* vect)
{
    dbl_array *smoothed, *step_heights;
    double *smoothedSlopes, mn, max_quot, cur_quot;
    int i, j, idx, max_quot_ind;
    
    smoothed = init_dbl_array(0, mgs_res->smoothed->cols + 1, 0);
    step_heights = init_dbl_array(0, result->meanlist->cols - 1, 0);

    memcpy(result->meanlist->values[result->meanlist->rows-1], vect->values, vect->length * sizeof(double));

    for(i = 0; i < result->v->length; i++)
    {
        smoothedSlopes = mgs_res->smoothed->values[q_res->index->values[i]-1];

        smoothed->values[0] = vect->values[0];
        for(j = 1; j < smoothed->length; j++)
        {
            smoothed->values[j] = smoothed->values[j-1] + smoothedSlopes[j-1];
        }

        memcpy(result->smoothedX->values[i], smoothed->values, smoothed->length * sizeof(double));

        for(j=0; j < q_res->steps->cols && q_res->steps->values[i][j]; j++)
        {
            if(j == 0)
            {
                result->meanlist->values[i][j] = mean(smoothed->values, 1, q_res->steps->values[i][j]);
            }
            else
            {
                result->meanlist->values[i][j] = mean(smoothed->values, q_res->steps->values[i][j-1]+1, q_res->steps->values[i][j]);
                step_heights->values[j-1] = result->meanlist->values[i][j] - result->meanlist->values[i][j-1];
                H->values[i][j-1] = step_heights->values[j-1];
            }
        }
        //now j is the index of the first 0 value or j == length(quant$steps[[z]]) + 1
        if(j <= q_res->steps->cols)
        {
            result->meanlist->values[i][j] = mean(smoothed->values, q_res->steps->values[i][j-1]+1, smoothed->length);
            step_heights->values[j-1] = result->meanlist->values[i][j] - result->meanlist->values[i][j-1];
            H->values[i][j-1] = step_heights->values[j-1];
        }

        max_quot = -1.0;
        max_quot_ind = -1;
        for(j = 0; j < q_res->steps->cols && q_res->steps->values[i][j]; j++)
        {
            idx = q_res->steps->values[i][j];
            mn = (smoothed->values[idx] + smoothed->values[idx - 1]) * 0.5;
            cur_quot = step_heights->values[j] / cost_Scalespace(smoothed, 0, smoothed->length - 1, mn);
            if(cur_quot > max_quot)
            {
                max_quot = cur_quot;
                max_quot_ind = j;
            }
        }

        result->v->values[i] = q_res->steps->values[i][max_quot_ind];
    }

    destroy_dbl_array(smoothed);
    destroy_dbl_array(step_heights);
}

double cost_Scalespace(dbl_array* vect, int a, int b, double y)
{
    int i;
    double cost, cost_root;

    cost = 0.0;
    for(i = a; i <= b; i++)
    {
        cost_root = vect->values[i] - y;
        cost += cost_root * cost_root;
    }

    return cost;
}

SEXP binarizeBASCB(SEXP vect, SEXP tau, SEXP numberofSamples, SEXP sigma)
{
    int value_count, sigma_count;
    mgs_result mgs_res;
    quant_result q_res;
    dbl_array *vector, *vect_sorted, *s;
    dbl_matrix *H_Mat;
    calc_V_result v_res;
    final_result f_res;
    SEXP result, binarized_vector, threshold, p_value, other_results, names;
    SEXP smoothed, zerocrossing, deriv, steps, H, index, v_vec, meanlist, smoothedX;

    value_count = length(vect);
    sigma_count = length(sigma);

    PROTECT(result = allocVector(VECSXP, 4));
    PROTECT(names = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(names,0, mkChar("binarized_vector"));
    SET_VECTOR_ELT(names,1, mkChar("threshold"));
    SET_VECTOR_ELT(names,2, mkChar("p_value"));
    SET_VECTOR_ELT(names,3, mkChar("other_results"));
    setAttrib(result, R_NamesSymbol, names);
    UNPROTECT(1);

    PROTECT(other_results = allocVector(VECSXP, 9));
    PROTECT(names = allocVector(VECSXP, 9));
    SET_VECTOR_ELT(names,0, mkChar("smoothed"));
    SET_VECTOR_ELT(names,1, mkChar("zerocrossing"));
    SET_VECTOR_ELT(names,2, mkChar("deriv"));
    SET_VECTOR_ELT(names,3, mkChar("steps"));
    SET_VECTOR_ELT(names,4, mkChar("H_Mat"));
    SET_VECTOR_ELT(names,5, mkChar("index"));
    SET_VECTOR_ELT(names,6, mkChar("v_vec"));
    SET_VECTOR_ELT(names,7, mkChar("smoothedX"));
    SET_VECTOR_ELT(names,8, mkChar("meanlist"));
    setAttrib(other_results, R_NamesSymbol, names);
    UNPROTECT(1);

    PROTECT(smoothed = allocMatrix(REALSXP, value_count - 1, sigma_count));
    PROTECT(zerocrossing = allocMatrix(INTSXP, (int)ceil((double)value_count / 2.0), sigma_count));
    PROTECT(deriv = allocVector(REALSXP, value_count - 1));
    
    mgs_res.smoothed = init_dbl_matrix(REAL(smoothed), sigma_count, value_count - 1, 0);
    mgs_res.zerocrossing = init_int_matrix(INTEGER(zerocrossing), sigma_count, (int)ceil((double)value_count / 2.0), 0);
    mgs_res.deriv = init_dbl_array(REAL(deriv), value_count - 1, 0);
    vector = init_dbl_array(REAL(vect), value_count, 1);
    vect_sorted = init_dbl_array(0, value_count, 0);
    s = init_dbl_array(REAL(sigma), sigma_count, 1);
    b = init_dbl_matrix(0, sigma_count, mgs_res.deriv->length, 0);
    b_returned = init_int_matrix(0, sigma_count, mgs_res.deriv->length, 0);

    memcpy(vect_sorted->values, REAL(vect), vect_sorted->length * sizeof(double));
    qsort(vect_sorted->values, vect_sorted->length, sizeof(double), comp);

    mgs(&mgs_res, vect_sorted, s);

    q_res.steps = init_int_matrix(0, sigma_count, (int)ceil((double)value_count / 2.0), 0);
    q_res.index = init_int_array(0, sigma_count, 0);
    q_res.greatest_index_ind = 0;
    q_res.greatest_steps_col = 0;
    q_res.greatest_steps_row = 0;

    getQuantizations(&q_res, &mgs_res);
    
    PROTECT(steps = allocMatrix(INTSXP, q_res.greatest_steps_col, q_res.greatest_steps_row));
    PROTECT(H = allocMatrix(REALSXP, q_res.greatest_steps_col, q_res.greatest_steps_row));
    PROTECT(index = allocVector(INTSXP, q_res.greatest_index_ind));

    cut_int_matrix(q_res.steps, INTEGER(steps), 0, q_res.greatest_steps_row - 1, 0, q_res.greatest_steps_col - 1);
    cut_int_array(q_res.index, INTEGER(index), 0, q_res.greatest_index_ind - 1);
    H_Mat = init_dbl_matrix(REAL(H), q_res.greatest_steps_row, q_res.greatest_steps_col, 0);

    PROTECT(v_vec = allocVector(INTSXP, q_res.index->length));
    PROTECT(smoothedX = allocMatrix(REALSXP, mgs_res.smoothed->cols + 1, q_res.index->length));
    PROTECT(meanlist = allocMatrix(REALSXP, vect_sorted->length, q_res.index->length + 1));
    v_res.v = init_int_array(INTEGER(v_vec), q_res.index->length, 0);
    
    v_res.meanlist = init_dbl_matrix(REAL(meanlist), q_res.index->length + 1, vect_sorted->length, 0);
    v_res.smoothedX = init_dbl_matrix(REAL(smoothedX), q_res.index->length, mgs_res.smoothed->cols + 1, 0);

    calc_V_Scalespace(&v_res, &mgs_res, &q_res, H_Mat, vect_sorted);

    PROTECT(binarized_vector = allocVector(INTSXP, value_count));
    PROTECT(threshold = allocVector(REALSXP, 1));
    PROTECT(p_value = allocVector(REALSXP, 1));

    f_res.binarized_vector = init_int_array(INTEGER(binarized_vector), value_count, 0);
    f_res.p = REAL(p_value);
    f_res.threshold = REAL(threshold);

    calc_final_results(&f_res, v_res.v, vector, vect_sorted, *REAL(tau), *INTEGER(numberofSamples));

    SET_VECTOR_ELT(other_results, 0, smoothed);
    SET_VECTOR_ELT(other_results, 1, zerocrossing);
    SET_VECTOR_ELT(other_results, 2, deriv);
    SET_VECTOR_ELT(other_results, 3, steps);
    SET_VECTOR_ELT(other_results, 4, H);
    SET_VECTOR_ELT(other_results, 5, index);
    SET_VECTOR_ELT(other_results, 6, v_vec);
    SET_VECTOR_ELT(other_results, 7, smoothedX);
    SET_VECTOR_ELT(other_results, 8, meanlist);
    
    SET_VECTOR_ELT(result, 0, binarized_vector);
    SET_VECTOR_ELT(result, 1, threshold);
    SET_VECTOR_ELT(result, 2, p_value);
    SET_VECTOR_ELT(result, 3, other_results);

    destroy_dbl_matrix(mgs_res.smoothed);
    destroy_int_matrix(mgs_res.zerocrossing);
    destroy_dbl_array(mgs_res.deriv);
    destroy_dbl_array(vector);
    destroy_dbl_array(vect_sorted);
    destroy_dbl_array(s);
    destroy_dbl_matrix(b);
    destroy_int_matrix(b_returned);
    b = 0;
    b_returned = 0;
    destroy_dbl_matrix(H_Mat);
    destroy_int_matrix(q_res.steps);
    destroy_int_array(q_res.index);
    destroy_int_array(v_res.v);
    destroy_dbl_matrix(v_res.meanlist);
    destroy_dbl_matrix(v_res.smoothedX);
    destroy_int_array(f_res.binarized_vector);

    UNPROTECT(14);

    return result;
}
