
#include "common.h"
#include <Rinternals.h>
#include <R.h>
#include "binarizeBASCA.h"
#include "binarizeBASCB.h"
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <time.h>

//sames as storing the compted costs, here for the approximation error values
dbl_array* e_tri_min = 0;
int_array* e_returned_tri_min = 0;

double calc_jump_height_tri_min(int_matrix* P, dbl_array* vect, int i, int k, int j){

	double height = 0.0;

	if(k == -1){
		//compute the jump height according to the definition from the paper
		if(i == 0){
			height = mean_BASCA(vect->values, P->values[j][i],  P->values[j][i+1]-1);
			height -= mean_BASCA(vect->values, 0, P->values[j][i]-1);
		}else{
			height = mean_BASCA(vect->values, P->values[j][i], P->values[j][i+1]-1);
			height -= mean_BASCA(vect->values, P->values[j][i-1], P->values[j][i]-1);
		}
	}else if(i == -1){
		if(k == j){
			height = mean_BASCA(vect->values, P->values[j][k], vect->length - 1);
			height -= mean_BASCA(vect->values, P->values[j][k-1], P->values[j][k]-1);
		}else{
			height = mean_BASCA(vect->values, P->values[j][k], P->values[j][k+1]-1);
			height -= mean_BASCA(vect->values, P->values[j][k-1], P->values[j][k]-1);
		}
	}
	return height;
}

/*
 * Calculates the approximation error for the step at index P_i(j)
 */
double calc_error_tri_min(int_matrix* P, dbl_array* vect, int i, int k, int j){

	double z1, z2, err, e_root, tmp, tmp2;
	int d, index;

	//check if value was computed before
	index = P->values[j][i]-1;
	if(e_tri_min && e_returned_tri_min && e_returned_tri_min->values[index])
	{
		e_returned_tri_min->values[index]++;
		return e_tri_min->values[index];
	}

	//compute the approximation error
	z1 = (vect->values[P->values[j][i]-1] + vect->values[P->values[j][i]]) * 0.5;
	z2 = (vect->values[P->values[j][k]-1] + vect->values[P->values[j][k]]) * 0.5;

	err = 0.0;
	for(d = 0; d < vect->length; d++){

		if(d <= i){
			e_root = vect->values[d] - z1;
			err += e_root * e_root;
		}else if(d > k){
			e_root = vect->values[d] - z2;
			err += e_root * e_root;
		}else{
			e_root = vect->values[d] - z1;
			tmp = e_root*e_root;
			e_root = vect->values[d] - z2;
			tmp2 = e_root*e_root;
			
			if(tmp < tmp2){
				err += tmp;
			}else{
				err += tmp2;
			}
		}
	}

	//store computed value
	if(e_tri_min && e_returned_tri_min)
	{
		e_returned_tri_min->values[index]++;
		e_tri_min->values[index] = err;
	}

	return err;
}

/*
 * Calculates the score for a certain step
 */
double calc_score_tri_min(int_matrix* P, dbl_matrix* H1, dbl_matrix* H2, dbl_array* vect, int i, int k, int j){

	double score,score1,score2;

	score1 = calc_jump_height_tri_min(P, vect, i, -1, j);
	score2 = calc_jump_height_tri_min(P, vect, -1, k, j);
	H1->values[j][i] = score1;
	H2->values[j][k] = score2;
	score = score1+score2;
	score /= 2*calc_error_tri_min(P, vect, i, k, j);

	return score;
}

/*
 * Calculates the Vector V which contains the indices of the strongest discontinuities for each amount of steps
 */
void calc_V_tri_min(int_matrix* v, dbl_array* Q_max, dbl_matrix* Q, dbl_matrix* H1, dbl_matrix* H2, int_matrix* P, dbl_array* vect){

	int i,k,j,max_ind1,max_ind2;
	double max_val,cur_val;

	for(j = 1; j < P->rows; j++){

		max_val = -1.0;
		max_ind1 = -1;
		max_ind2 = -1;
		
		for(i = 0; i <= (j-1); i++){

			for(k = (i+1); k <= j; k++){
			
				cur_val = calc_score_tri_min(P, H1, H2, vect, i, k, j);
				Q->values[j][i] = cur_val;
				if(cur_val > max_val)
				{
					max_val = cur_val;
					max_ind1 = i;
					max_ind2 = k;
				}
			}
		}

		v->values[j][0] = P->values[j][max_ind1];
		v->values[j][1] = P->values[j][max_ind2];

		Q_max->values[j] = max_val;
	}
}

/*
 * Calculates the final three results (binarized_vector, threshold, p-value).
 */
void calc_final_results_tri_min(final_result_tri* result, int_matrix* v, dbl_array* vect, dbl_array* vect_sorted, double tau, int numberofsamples)
{
	int i, k, c;
	int_array* samples, * v_col1, * v_col2;
	double nom, nom1, nom2, t_zero, t_star, mdm, mdm1, mdm2;

	v_col1 = init_int_array(0, v->rows-1, 0);
	v_col2 = init_int_array(0, v->rows-1, 0);

	for(c = 1; c < v->rows; c++){
		v_col1->values[c-1] = v->values[c][0];
		v_col2->values[c-1] = v->values[c][1];
	}

	//calculate the threshold and the binarized vector:
	i = (int)floor(median(v_col1));
	*(result->threshold1) = (vect_sorted->values[i] + vect_sorted->values[i-1]) * 0.5;

	k = (int)floor(median(v_col2));
	*(result->threshold2) = (vect_sorted->values[k] + vect_sorted->values[k-1]) * 0.5;

	for(i = 0; i < result->binarized_vector->length; i++)
	{
		result->binarized_vector->values[i] = (int)(vect->values[i] > *(result->threshold1));
		if(vect->values[i] > *(result->threshold2)){
			result->binarized_vector->values[i] = 2;
		}
	}

	//calculate the statistics:
	samples = init_int_array(0, v->rows-1, 0);

	//nom = MDM(v')
	nom1 = normabsmedian(v_col1, vect_sorted);
	nom2 = normabsmedian(v_col2, vect_sorted);

	nom = (nom1 + nom2) / 2;

//    print_int_array(v_col1, "col1");
//    print_int_array(v_col2, "col2");

//    printf("norm med1: %f ; norm med2: %f ; nom: %f\n", nom1, nom2, nom);

	t_zero = tau - nom;

	*(result->p) = 1.0;

	if(v->rows < 3)
	{
		warning("Too few members in the vector of strongest discontinuities of the optimal step functions. The computed p-value may not be reliable.");
	}

	#if !DEBUG_MODE
	GetRNGstate();
	#endif
	for(i = 0; i < numberofsamples; i++)
	{
		//resample the values and calc t(v*) = MDM(v') - MDM(v*)
		blockwiseboot(samples, v_col1);
		mdm1 = normabsmedian(samples, vect_sorted);

		blockwiseboot(samples, v_col2);
		mdm2 = normabsmedian(samples, vect_sorted);
		
		mdm = (mdm1 + mdm2)/2;

//        printf("norm block med1: %f ; norm block med2: %f ; mdm: %f\n", mdm1, mdm2, mdm);
		
		t_star = nom - mdm;

		//sum up the number of t(v*) >= t0
		*(result->p) += (double)(t_star >= t_zero);
	}
	#if !DEBUG_MODE
	PutRNGstate();
	#endif

	//divide p by the number of samples, which is the maximal possible value for p
	//so p is in interval [0,1]
	*(result->p) /= (((double)numberofsamples) +1.0);

	destroy_int_array(samples);
	destroy_int_array(v_col1);
	destroy_int_array(v_col2);
}

/*
 * allocates memory for storing all the computed intermediate values, that could be reused
 */
int alloc_Accelerator_Memory_tri_min(int value_count){

	int bytes = 0;

	e_tri_min = init_dbl_array(0, value_count - 1, 0);
	bytes += value_count * sizeof(double);
	e_returned_tri_min = init_int_array(0, value_count - 1, 0);
	bytes += value_count * sizeof(int);

	return bytes;
}

/*
 * Free the memory used for accelerating computations. This method also sets all pointers to 0.
 */
void free_Accelerator_Memory_tri_min(void){

	destroy_dbl_array(e_tri_min);
	destroy_int_array(e_returned_tri_min);
	e_returned_tri_min = 0;
	e_tri_min = 0;
}

SEXP TASCA_min(SEXP vect, SEXP tau, SEXP numberofsamples){
	
	//get the lengths
	//int i,j,sum,sum_tot;
	//int bytes = 0;
	int value_count = length(vect);
	int vc_m1 = value_count - 1;
	int vc_m2 = vc_m1 - 1;
	dbl_array *vector, *vect_sorted, *Q_Max_vec;
	dbl_matrix *Cc_Mat, *Q_Mat, *H_Mat1, *H_Mat2;
	int_matrix *Ind_Mat, *P_Mat, *v_vec;
	final_result_tri f_res;

	//sort the vect into vect_sorted
	vector = init_dbl_array(REAL(vect), value_count, 1);

	vect_sorted = init_dbl_array(0, value_count, 0);
	memcpy(vect_sorted->values, vector->values, vect_sorted->length * sizeof(double));
	qsort(vect_sorted->values, vect_sorted->length, sizeof(double), comp);

	//name the required SEXP Objects
	SEXP result, binarized_vector, threshold1, threshold2, p_value, other_results, Cc, Ind, P, Q, H1, H2, Q_max, v, Names;

	//allocate memory for saving calculated values
	alloc_Accelerator_Memory(value_count);
	alloc_Accelerator_Memory_tri_min(value_count);
	
	//allocate the final result and set the names of the entries
	PROTECT(result = allocVector(VECSXP, 5));
	PROTECT(Names = allocVector(VECSXP, 5));
	SET_VECTOR_ELT(Names,0, mkChar("binarized_vector"));
	SET_VECTOR_ELT(Names,1, mkChar("threshold1"));
	SET_VECTOR_ELT(Names,2, mkChar("threshold2"));
	SET_VECTOR_ELT(Names,3, mkChar("p_value"));
	SET_VECTOR_ELT(Names,4, mkChar("other_results"));
	setAttrib(result, R_NamesSymbol, Names);
	UNPROTECT(1);

	PROTECT(other_results = allocVector(VECSXP, 8));
	PROTECT(Names = allocVector(VECSXP, 8));
	SET_VECTOR_ELT(Names,0, mkChar("Cc"));
	SET_VECTOR_ELT(Names,1, mkChar("Ind"));
	SET_VECTOR_ELT(Names,2, mkChar("P_Mat"));
	SET_VECTOR_ELT(Names,3, mkChar("Q_Mat"));
	SET_VECTOR_ELT(Names,4, mkChar("H_Mat1"));
	SET_VECTOR_ELT(Names,5, mkChar("H_Mat2"));
	SET_VECTOR_ELT(Names,6, mkChar("maximal_Qs"));
	SET_VECTOR_ELT(Names,7, mkChar("v_vec"));
	setAttrib(other_results, R_NamesSymbol, Names);
	UNPROTECT(1);


	//allocate memory for result matrices and vectors and set the matrix values to zero (because they aren't
	//all overwritten by the functions)
	PROTECT(binarized_vector = allocVector(INTSXP, value_count));
	PROTECT(threshold1 = allocVector(REALSXP, 1));
	PROTECT(threshold2 = allocVector(REALSXP, 1));
	PROTECT(p_value = allocVector(REALSXP, 1));
	PROTECT(Cc = allocMatrix(REALSXP, vc_m1, vc_m1));
	PROTECT(Ind = allocMatrix(INTSXP, vc_m1, vc_m2));
	PROTECT(P = allocMatrix(INTSXP, vc_m2, vc_m2));
	PROTECT(Q = allocMatrix(REALSXP, vc_m2, vc_m2));
	PROTECT(H1 = allocMatrix(REALSXP, vc_m2, vc_m2));
	PROTECT(H2 = allocMatrix(REALSXP, vc_m2, vc_m2));
	PROTECT(Q_max = allocVector(REALSXP, vc_m2));
	PROTECT(v = allocMatrix(INTSXP, vc_m2, 2));

	Cc_Mat = init_dbl_matrix(REAL(Cc), vc_m1, vc_m1, 0);
	Ind_Mat = init_int_matrix(INTEGER(Ind), vc_m2, vc_m1, 0);
	P_Mat = init_int_matrix(INTEGER(P), vc_m2, vc_m2, 0);
	v_vec = init_int_matrix(INTEGER(v), vc_m2, 2, 0);
	Q_Max_vec = init_dbl_array(REAL(Q_max), vc_m2, 0);
	Q_Mat = init_dbl_matrix(REAL(Q), vc_m2, vc_m2, 0);
	H_Mat1 = init_dbl_matrix(REAL(H1), vc_m2, vc_m2, 0);
	H_Mat2 = init_dbl_matrix(REAL(H2), vc_m2, vc_m2, 0);
	f_res.binarized_vector = init_int_array(INTEGER(binarized_vector), value_count, 0);
	f_res.p = REAL(p_value);
	f_res.threshold1 = REAL(threshold1);
	f_res.threshold2 = REAL(threshold2);

	//start the computation of the entries of all matrices
	calc_First_Cost_Matrix_Line(Cc_Mat, vect_sorted);
	calc_RestCc_and_Ind_Matrices(Cc_Mat, Ind_Mat, vect_sorted);
	calc_P_Matrix(P_Mat, Ind_Mat);
	calc_V_tri_min(v_vec, Q_Max_vec, Q_Mat, H_Mat1, H_Mat2, P_Mat, vect_sorted);
	
	//free the memory for calculated values
	free_Accelerator_Memory();
	free_Accelerator_Memory_tri_min();

	//calculate the final three results
	calc_final_results_tri_min(&f_res, v_vec, vector, vect_sorted, *REAL(tau), *INTEGER(numberofsamples));

	//free(vect_sorted);
	destroy_dbl_array(vector);
	destroy_dbl_array(vect_sorted);
	destroy_dbl_matrix(Cc_Mat);
	destroy_int_matrix(Ind_Mat);
	destroy_int_matrix(P_Mat);
	destroy_int_matrix(v_vec);
	destroy_dbl_array(Q_Max_vec);
	destroy_dbl_matrix(Q_Mat);
	destroy_dbl_matrix(H_Mat1);
	destroy_dbl_matrix(H_Mat2);
	destroy_int_array(f_res.binarized_vector);

	//assign the computed elements to the final result
	SET_VECTOR_ELT(other_results,0, Cc);
	SET_VECTOR_ELT(other_results,1, Ind);
	SET_VECTOR_ELT(other_results,2, P);
	SET_VECTOR_ELT(other_results,3, Q);
	SET_VECTOR_ELT(other_results,4, H1);
	SET_VECTOR_ELT(other_results,5, H2);
	SET_VECTOR_ELT(other_results,6, Q_max);
	SET_VECTOR_ELT(other_results,7, v);

	SET_VECTOR_ELT(result,0, binarized_vector);
	SET_VECTOR_ELT(result,1, threshold1);
	SET_VECTOR_ELT(result,2, threshold2);
	SET_VECTOR_ELT(result,3, p_value);
	SET_VECTOR_ELT(result,4, other_results);

	UNPROTECT(14);
	return result;
}



/*
*
*
*
*
*	TASC-B
*
*
*
*
*/

void calc_V_Scalespace_tri_min(calc_V_result_tri* result, mgs_result* mgs_res, quant_result* q_res, dbl_matrix* H1, dbl_matrix* H2, dbl_array* vect){
	dbl_array *smoothed,*step_heights,*step_heights2;
	double *smoothedSlopes, mn1, mn2, max_quot, cur_quot, tmp, tmp2, tmp3;
	int i, j, idx1, idx2, max_quot_ind1, max_quot_ind2, k;

	smoothed = init_dbl_array(0, mgs_res->smoothed->cols + 1, 0);
	step_heights = init_dbl_array(0, result->meanlist->cols - 1, 0);
	step_heights2 = init_dbl_array(0, result->meanlist->cols - 1, 0);

	memcpy(result->meanlist->values[result->meanlist->rows-1], vect->values, vect->length * sizeof(double));

	for(i = 1; i < result->v->rows; i++)
	{

		smoothedSlopes = mgs_res->smoothed->values[q_res->index->values[i]-1];

		smoothed->values[0] = vect->values[0];
		for(j = 1; j < smoothed->length; j++)
		{
			smoothed->values[j] = smoothed->values[j-1] + smoothedSlopes[j-1];
		}

		memcpy(result->smoothedX->values[i], smoothed->values, smoothed->length * sizeof(double));

		for(j=0; j < q_res->steps->cols && q_res->steps->values[i][j]; j++){
			
			if(j == 0){
				result->meanlist->values[i][j] = mean(smoothed->values, 1, q_res->steps->values[i][j]);
			}

			for(k = (j+1); k < q_res->steps->cols && q_res->steps->values[i][k]; k++){
				//printf("i: %d, j: %d ,k: %d\n", i, j, k);
				if(k ==1){
					result->meanlist->values[i][k] = mean(smoothed->values, q_res->steps->values[i][k-1]+1, q_res->steps->values[i][k]);
				}else{
					result->meanlist->values[i][k] = mean(smoothed->values, q_res->steps->values[i][k-1]+1, q_res->steps->values[i][k]);
					step_heights2->values[k-1] = result->meanlist->values[i][k] - result->meanlist->values[i][k-1];
					
					H2->values[i][k-1] = step_heights2->values[k-1];
				}

			}

			if(j > 0){
				step_heights->values[j-1] = result->meanlist->values[i][j] - result->meanlist->values[i][j-1];
				H1->values[i][j-1] = step_heights->values[j-1];
			}
			//now k is the index of the first 0 value or j == length(quant$steps[[z]]) + 1
			
			//printf("last k: %d <= %d\n", k, q_res->steps->cols);
			if(k <= q_res->steps->cols){
				result->meanlist->values[i][k] = mean(smoothed->values, q_res->steps->values[i][k-1]+1, smoothed->length);
				step_heights2->values[k-1] = result->meanlist->values[i][k] - result->meanlist->values[i][k-1];
				H2->values[i][k-1] = step_heights2->values[k-1];
			}
		}

		max_quot = -1.0;
		max_quot_ind1 = -1;
		max_quot_ind2 = -1;
/*
		print_dbl_matrix(H1, "H1");
		print_dbl_matrix(H2, "H2");

		print_dbl_array(step_heights, "step_heights");
		print_dbl_array(step_heights2, "step_heights2");

		print_int_matrix(q_res->steps, "qres steps");
*/
		for(j = 0; j < q_res->steps->cols && q_res->steps->values[i][j]; j++){
			
			for(k = (j+1); k < q_res->steps->cols && q_res->steps->values[i][k]; k++){
				
				idx1 = q_res->steps->values[i][j];
				mn1 = (smoothed->values[idx1] + smoothed->values[idx1 - 1]) * 0.5;
				idx2 = q_res->steps->values[i][k];
				mn2 = (smoothed->values[idx2] + smoothed->values[idx2 - 1]) * 0.5;
				
				cur_quot = (H1->values[i][j] + H2->values[i][k]);
				
				tmp = (cost_Scalespace(smoothed, 0, (j-1), mn1)
					+ cost_Scalespace(smoothed, k, smoothed->length - 1, mn2));
				tmp2 = cost_Scalespace(smoothed, j, k - 1, mn1);
				tmp3 = cost_Scalespace(smoothed, j, k - 1, mn2);

				if(tmp2 < tmp3){
					cur_quot /= 2*(tmp+tmp2);
				}else{
					cur_quot /= 2*(tmp+tmp3);
				}
				
//				printf("i: %d, j: %d ,k: %d, idx1: %d, idx2: %d, cur_quot: %f\n", i, j, k,idx1, idx2,cur_quot);
				if(cur_quot > max_quot){
					max_quot = cur_quot;
					max_quot_ind1 = j;
					max_quot_ind2 = k;
				}
			}
		}
//		printf("max_quot_ind1: %d; max_quot_ind2: %d\n", max_quot_ind1, max_quot_ind2);
		result->v->values[i][0] = q_res->steps->values[i][max_quot_ind1];
		result->v->values[i][1] = q_res->steps->values[i][max_quot_ind2];
	}

/*	print_dbl_array(vect, "vect");
	print_dbl_array(smoothed, "smoothed");
	print_dbl_matrix(result->meanlist, "result meanlist");
*/
	destroy_dbl_array(smoothed);
	destroy_dbl_array(step_heights);
	destroy_dbl_array(step_heights2);
}

SEXP TASCB_min(SEXP vect, SEXP tau, SEXP numberofSamples, SEXP sigma){
	int value_count, sigma_count;
	mgs_result mgs_res;
	quant_result q_res;
	dbl_array *vector, *vect_sorted, *s;
	dbl_matrix *H_Mat1, *H_Mat2;
	calc_V_result_tri v_res;
	final_result_tri f_res;
	SEXP result, binarized_vector, threshold1, threshold2, p_value, other_results, names;
	SEXP smoothed, zerocrossing, deriv, steps, H1, H2, index, v_vec, meanlist, smoothedX;

	value_count = length(vect);
	sigma_count = length(sigma);

	PROTECT(result = allocVector(VECSXP, 5));
	PROTECT(names = allocVector(VECSXP, 5));
	SET_VECTOR_ELT(names,0, mkChar("binarized_vector"));
	SET_VECTOR_ELT(names,1, mkChar("threshold1"));
	SET_VECTOR_ELT(names,2, mkChar("threshold2"));
	SET_VECTOR_ELT(names,3, mkChar("p_value"));
	SET_VECTOR_ELT(names,4, mkChar("other_results"));
	setAttrib(result, R_NamesSymbol, names);
	UNPROTECT(1);

	PROTECT(other_results = allocVector(VECSXP, 10));
	PROTECT(names = allocVector(VECSXP, 10));
	SET_VECTOR_ELT(names,0, mkChar("smoothed"));
	SET_VECTOR_ELT(names,1, mkChar("zerocrossing"));
	SET_VECTOR_ELT(names,2, mkChar("deriv"));
	SET_VECTOR_ELT(names,3, mkChar("steps"));
	SET_VECTOR_ELT(names,4, mkChar("H_Mat1"));
	SET_VECTOR_ELT(names,5, mkChar("H_Mat2"));
	SET_VECTOR_ELT(names,6, mkChar("index"));
	SET_VECTOR_ELT(names,7, mkChar("v_vec"));
	SET_VECTOR_ELT(names,8, mkChar("smoothedX"));
	SET_VECTOR_ELT(names,9, mkChar("meanlist"));
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
	PROTECT(H1 = allocMatrix(REALSXP, q_res.greatest_steps_col, q_res.greatest_steps_row));
	PROTECT(H2 = allocMatrix(REALSXP, q_res.greatest_steps_col, q_res.greatest_steps_row));
	PROTECT(index = allocVector(INTSXP, q_res.greatest_index_ind));

	cut_int_matrix(q_res.steps, INTEGER(steps), 0, q_res.greatest_steps_row - 1, 0, q_res.greatest_steps_col - 1);
	cut_int_array(q_res.index, INTEGER(index), 0, q_res.greatest_index_ind - 1);
	H_Mat1 = init_dbl_matrix(REAL(H1), q_res.greatest_steps_row, q_res.greatest_steps_col, 0);
	H_Mat2 = init_dbl_matrix(REAL(H2), q_res.greatest_steps_row, q_res.greatest_steps_col, 0);

	PROTECT(v_vec = allocMatrix(INTSXP, q_res.index->length, 2));
	PROTECT(smoothedX = allocMatrix(REALSXP, mgs_res.smoothed->cols + 1, q_res.index->length));
	PROTECT(meanlist = allocMatrix(REALSXP, vect_sorted->length, q_res.index->length + 1));
	v_res.v = init_int_matrix(INTEGER(v_vec), q_res.index->length, 2, 0);
	
	v_res.meanlist = init_dbl_matrix(REAL(meanlist), q_res.index->length + 1, vect_sorted->length, 0);
	v_res.smoothedX = init_dbl_matrix(REAL(smoothedX), q_res.index->length, mgs_res.smoothed->cols + 1, 0);

	calc_V_Scalespace_tri_min(&v_res, &mgs_res, &q_res, H_Mat1, H_Mat2, vect_sorted);

	PROTECT(binarized_vector = allocVector(INTSXP, value_count));
	PROTECT(threshold1 = allocVector(REALSXP, 1));
	PROTECT(threshold2 = allocVector(REALSXP, 1));
	PROTECT(p_value = allocVector(REALSXP, 1));

	f_res.binarized_vector = init_int_array(INTEGER(binarized_vector), value_count, 0);
	f_res.p = REAL(p_value);
	f_res.threshold1 = REAL(threshold1);
	f_res.threshold2 = REAL(threshold2);

	calc_final_results_tri_min(&f_res, v_res.v, vector, vect_sorted, *REAL(tau), *INTEGER(numberofSamples));

	SET_VECTOR_ELT(other_results, 0, smoothed);
	SET_VECTOR_ELT(other_results, 1, zerocrossing);
	SET_VECTOR_ELT(other_results, 2, deriv);
	SET_VECTOR_ELT(other_results, 3, steps);
	SET_VECTOR_ELT(other_results, 4, H1);
	SET_VECTOR_ELT(other_results, 5, H2);
	SET_VECTOR_ELT(other_results, 6, index);
	SET_VECTOR_ELT(other_results, 7, v_vec);
	SET_VECTOR_ELT(other_results, 8, smoothedX);
	SET_VECTOR_ELT(other_results, 9, meanlist);
	
	SET_VECTOR_ELT(result, 0, binarized_vector);
	SET_VECTOR_ELT(result, 1, threshold1);
	SET_VECTOR_ELT(result, 2, threshold2);
	SET_VECTOR_ELT(result, 3, p_value);
	SET_VECTOR_ELT(result, 4, other_results);

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
	destroy_dbl_matrix(H_Mat1);
	destroy_dbl_matrix(H_Mat2);
	destroy_int_matrix(q_res.steps);
	destroy_int_array(q_res.index);
	destroy_int_matrix(v_res.v);
	destroy_dbl_matrix(v_res.meanlist);
	destroy_dbl_matrix(v_res.smoothedX);
	destroy_int_array(f_res.binarized_vector);

	UNPROTECT(16);

	return result;
}


