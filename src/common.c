#include "common.h"
#if DEBUG_MODE
    #include "mersennetwister.h"
#else
    #include <Rinternals.h>
    #include <R.h>
    #include <Rmath.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

//int is_called_from_R = 1;

/*
 * Comparator Function which is necessary for applying qsort on a double array
 */
int comp(const void* a, const void* b)
{
    double val = (double)(*(double*)a) - (double)(*(double*)b);
    if(val < 0.0)return -1;
    else if(val == 0.0)return 0;
    else return 1;
}

/*
 * Comparator Function which is necessary for applying qsort on a integer array
 */
int comp_int(const void* a, const void* b)
{
    return (int)(*(int*)a) - (int)(*(int*)b);
}

int comp_desc(const void* a, const void* b)
{
    return comp(b, a);
}

int comp_desc_int(const void* a, const void* b)
{
    return comp_int(b, a);
}

int indexOf(double value, dbl_array* array)
{
    int i;
    for(i=0;i<array->length;i++)
    {
        if(array->values[i] == value)
            return i;
    }
    return -1;
}

int indexOf_int(int value, int_array* array)
{
    int i;
    for(i=0;i<array->length;i++)
    {
        if(array->values[i] == value)
            return i;
    }
    return -1;
}
dbl_matrix* init_dbl_matrix(double* values, int rows, int cols, int keep)
{
    int i,val_ind;
    double *vals;
    dbl_matrix* mat;

    mat = malloc(sizeof(dbl_matrix));
    mat->values = malloc(rows * sizeof(double*));

    vals = values ? values : malloc(rows * cols * sizeof(double));
    
    if (vals == 0)
      Rf_error("Could not allocate memory!");

    if(!values || !keep)
        memset(vals, 0, rows * cols * sizeof(double));

    mat->allocated_values = !values;

    for(i = 0, val_ind = 0; i < rows; i++, val_ind += cols)
    {
        mat->values[i] = &vals[val_ind];
    }
    
    mat->cols = cols;
    mat->rows = rows;

    return mat;
}

/*dbl_matrix* init_dbl_matrix_keep(double* values, int rows, int cols)
{
    int i,val_ind;
    double *vals;
    dbl_matrix* mat;

    mat = malloc(sizeof(dbl_matrix));
    mat->values = malloc(rows * sizeof(double*));

    vals = values ? values : malloc(rows * cols * sizeof(double));

    mat->allocated_values = !values;

    for(i = 0, val_ind = 0; i < rows; i++, val_ind += cols)
    {
        mat->values[i] = &vals[val_ind];
    }

    mat->cols = cols;
    mat->rows = rows;

    return mat;
}*/

int_matrix* init_int_matrix(int* values, int rows, int cols, int keep)
{
    int i,val_ind, *vals;
    int_matrix* mat;

    mat = malloc(sizeof(int_matrix));
    mat->values = malloc(rows * sizeof(int*));

    vals = values ? values : malloc(rows * cols * sizeof(int));

    if (vals == 0)
      Rf_error("Could not allocate memory!");

    if(!values || !keep)
        memset(vals, 0, rows * cols * sizeof(int));

    mat->allocated_values = !values;

    for(i = 0, val_ind = 0; i < rows; i++, val_ind += cols)
    {
        mat->values[i] = &vals[val_ind];
    }

    mat->cols = cols;
    mat->rows = rows;

    return mat;
}

/*int_matrix* init_int_matrix_keep(int* values, int rows, int cols)
{
    int i,val_ind, *vals;
    int_matrix* mat;

    mat = malloc(sizeof(int_matrix));
    mat->values = malloc(rows * sizeof(int*));

    vals = values ? values : malloc(rows * cols * sizeof(int));

    mat->allocated_values = !values;

    for(i = 0, val_ind = 0; i < rows; i++, val_ind += cols)
    {
        mat->values[i] = &vals[val_ind];
    }

    mat->cols = cols;
    mat->rows = rows;

    return mat;
}*/

dbl_array* init_dbl_array(double* values, int length, int keep)
{
    dbl_array* arr;

    arr = malloc(sizeof(dbl_array));
    arr->values = values ? values : malloc(length * sizeof(double));
    
    if (arr->values == 0)
      Rf_error("Could not allocate memory!");
    
    if(!values || !keep)
        memset(arr->values, 0, length * sizeof(double));
    arr->allocated_values = !values;
    arr->length = length;

    return arr;
}

/*dbl_array* init_dbl_array_keep(double* values, int length)
{
    dbl_array* arr;

    arr = malloc(sizeof(dbl_array));
    arr->values = values ? values : malloc(length * sizeof(double));
    arr->allocated_values = !values;
    arr->length = length;

    return arr;
}*/

int_array* init_int_array(int* values, int length, int keep)
{
    int_array* arr;

    arr = malloc(sizeof(int_array));
    arr->values = values ? values : malloc(length * sizeof(int));
    
    if (arr->values == 0)
      Rf_error("Could not allocate memory!");
    
    if(!values || !keep)
        memset(arr->values, 0, length * sizeof(int));
    arr->allocated_values = !values;
    arr->length = length;

    return arr;
}

/*int_array* init_int_array_keep(int* values, int length)
{
    int_array* arr;

    arr = malloc(sizeof(int_array));
    arr->values = values ? values : malloc(length * sizeof(int));
    arr->allocated_values = !values;
    arr->length = length;

    return arr;
}*/

void cut_dbl_matrix(dbl_matrix* mat, double* values, int row_begin, int row_end, int col_begin, int col_end)
{
    int new_rows, new_cols, i;
    double** new_values;

    new_rows = row_end - row_begin + 1;
    new_cols = col_end - col_begin + 1;

    new_values = malloc(new_rows * sizeof(double*));
    new_values[0] = values ? values : malloc(new_cols * new_rows * sizeof(double));

    for(i=0;i<new_rows;i++)
    {
        if(i>0)
            new_values[i] = &new_values[0][i*new_cols];
        memcpy(new_values[i], &mat->values[row_begin+i][col_begin], new_cols * sizeof(double));
    }

    if(mat->allocated_values)
        free(mat->values[0]);
    free(mat->values);

    mat->values = new_values;
    mat->cols = new_cols;
    mat->rows = new_rows;
    mat->allocated_values = !values;
}

void cut_int_matrix(int_matrix* mat, int* values, int row_begin, int row_end, int col_begin, int col_end)
{
    int new_rows, new_cols, i, **new_values;

    new_rows = row_end - row_begin + 1;
    new_cols = col_end - col_begin + 1;

    new_values = malloc(new_rows * sizeof(int*));
    new_values[0] = values ? values : malloc(new_cols * new_rows * sizeof(int));

    for(i=0;i<new_rows;i++)
    {
        if(i>0)
            new_values[i] = &new_values[0][i*new_cols];
        memcpy(new_values[i], &mat->values[row_begin+i][col_begin], new_cols * sizeof(int));
    }

    if(mat->allocated_values)
        free(mat->values[0]);
    free(mat->values);

    mat->values = new_values;
    mat->cols = new_cols;
    mat->rows = new_rows;
    mat->allocated_values = !values;
}

void cut_dbl_array(dbl_array* arr, double* values, int begin, int end)
{
    int new_length;
    double* new_values;

    new_length = end - begin + 1;
    new_values = values ? values : malloc(new_length * sizeof(double));

    memcpy(new_values, &arr->values[begin], new_length * sizeof(double));

    if(arr->allocated_values)
        free(arr->values);

    arr->values = new_values;
    arr->length = new_length;
    arr->allocated_values = !values;
}

void cut_int_array(int_array* arr, int* values, int begin, int end)
{
    int new_length, *new_values;

    new_length = end - begin + 1;
    new_values = values ? values : malloc(new_length * sizeof(int));

    memcpy(new_values, &arr->values[begin], new_length * sizeof(int));

    if(arr->allocated_values)
        free(arr->values);

    arr->values = new_values;
    arr->length = new_length;
    arr->allocated_values = !values;
}

void destroy_dbl_matrix(dbl_matrix* mat)
{
    if(mat->allocated_values)
        free(mat->values[0]);
    free(mat->values);
    free(mat);
}

void destroy_int_matrix(int_matrix* mat)
{
    if(mat->allocated_values)
        free(mat->values[0]);
    free(mat->values);
    free(mat);
}

void destroy_dbl_array(dbl_array* arr)
{
    if(arr->allocated_values)
        free(arr->values);
    free(arr);
}

void destroy_int_array(int_array* arr)
{
    if(arr->allocated_values)
        free(arr->values);
    free(arr);
}


void print_dbl_array(dbl_array* arr, const char* name)
{
    int i;
    PRINT("%s = [%f", name, arr->values[0]);
    for(i=1;i<arr->length;i++)
        PRINT(", %f", arr->values[i]);
    PRINT("]\n");
}

void print_int_array(int_array* arr, const char* name)
{
    int i;
    PRINT("%s = [%d", name, arr->values[0]);
    for(i=1;i<arr->length;i++)
        PRINT(", %d", arr->values[i]);
    PRINT("]\n");
}

void print_dbl_matrix(dbl_matrix* mat, const char* name)
{
    int i,j;
    PRINT("%s = \n", name);
    for(i = 0;i < mat->rows; i++)
    {
        for(j = 0;j < mat->cols; j++)
        {
            PRINT("%f ", mat->values[i][j]);
        }
        PRINT("\n");
    }
}

void print_int_matrix(int_matrix* mat, const char* name)
{
    int i,j;
    PRINT("%s = \n", name);
    for(i = 0;i < mat->rows; i++)
    {
        for(j = 0;j < mat->cols; j++)
        {
            PRINT("%d ", mat->values[i][j]);
        }
        PRINT("\n");
    }
}


double mean(double* values, int a, int b)
{
    int i, begin, end;
    double mean_val;

    begin = a - 1;
    end = b - 1;

    //compute the mean
    mean_val = 0.0;
    for(i = begin; i <= end; i++){
        mean_val += values[i];
    }
    mean_val /= ((double)end-(double)begin+1.0);

    return mean_val;
}

/*
 * Calculates the median of an array
 */
double median(int_array* vect)
{
    double result;
    int index;
    int_array* v_sorted;

    //allocate memory for sorting the values, then copy the values and sort them
    v_sorted = init_int_array(0, vect->length, 0);
    memcpy(v_sorted->values, vect->values, v_sorted->length * sizeof(int));
    qsort(v_sorted->values, v_sorted->length, sizeof(int), comp_int);

    //The median of the array is in the middle of the sorted array. If there are an even number of elements
    //the median is the mean of the 2 elements in the middle of the sorted array
    index = v_sorted->length / 2;
    if(v_sorted->length == 1)
        result = (double)v_sorted->values[0];
    else if((v_sorted->length % 2) == 1)
        result = (double)v_sorted->values[index];
    else
        result = ((double)v_sorted->values[index] + (double)v_sorted->values[index - 1]) * 0.5;

    //free the allocated memory
    destroy_int_array(v_sorted);
    return result;
}

/*
 * Computes the average deviation for an array from its median
 */
double normabsmedian(int_array* vect,dbl_array* vector_sorted)
{
    int i;
    double median_val, result, mean_val;
    dbl_array* mem;

    //allocate memory
    mem = init_dbl_array(0, vect->length, 0);

    //compute the median of the array
    median_val = median(vect);

    //sum the deviations of the elements from the median and divide it by the number of elements plus one(for normalizing the values
    //to maximal 1)
    for(i = 0; i < mem->length; i++)
    {
        mem->values[i] = fabs(median_val - (double)vect->values[i]);
    }

    mean_val = mean(mem->values, 1, mem->length);
    
    result = mean_val / (vector_sorted->length - 1);

    //free the memory
    destroy_dbl_array(mem);

    return result;
}

/*
 * Resamples the original values. It takes ceil(#elements / bl) blocks of length bl = round((#elements)^(0.25))+1 out of the original
 * values and concatenates them to a new vector
 */
void blockwiseboot(int_array* result, int_array* originals)
{
    int bl, i, j, sample_count, index;
    double temp, random, max;

    //bl = round(#elements^0.25) + 1
    temp = sqrt(sqrt((double)result->length));
    bl = ROUND(temp) + 1;

    sample_count = (int)ceil((double)result->length/(double)bl);

    index = 0;
    //max = the maximal acceptable index for getting a block of length bl.
    //Plus 1 because of random generated values are in the interval [-0.5, max + 0.5) so every
    //possible index has the same probability after rounding the random value
    max = (double)(result->length - bl);

    for(i = 0; i < sample_count; i++)
    {
        //get a random number out of interval [-0.5, max + 0.5)
        random = RANDOM(-0.5, max + 0.5);

        //round the random number, clamping is only for array-access safety
        random = ROUND(random);
        random = CLAMP(random, 0.0, max);

        //Get a maximum of bl values from the original values starting at index random. If index == #elements the resampling is complete
        for(j = 0; j < bl && index < result->length; j++, index++)
        {
            //printf("INDEX: %d; IND: %d, ORIGINALS: %d\n", index, (int)random + j, originals->values[(int)random + j]);
            result->values[index] = originals->values[(int)random + j];
        }
    }
}

/*
 * Calculates the final three results (binarized_vector, threshold, p-value).
 */
void calc_final_results(final_result* result, int_array* v, dbl_array* vect, dbl_array* vect_sorted, double tau, int numberofsamples)
{
    int i;
    int_array* samples;
    double nom, t_zero, t_star, mdm;

    //calculate the threshold and the binarized vector:
    i = (int)floor(median(v));
    *(result->threshold) = (vect_sorted->values[i] + vect_sorted->values[i-1]) * 0.5;
    for(i = 0; i < result->binarized_vector->length; i++)
    {
        result->binarized_vector->values[i] = (int)(vect->values[i] > *(result->threshold));
    }

    //calculate the statistics:
    samples = init_int_array(0, v->length, 0);

    //nom = MDM(v')
    nom = normabsmedian(v, vect_sorted);

    t_zero = tau - nom;

    *(result->p) = 1.0;

    if(v->length < 3)
    {
        warning("Too few members in the vector of strongest discontinuities of the optimal step functions. The computed p-value may not be reliable.");
    }

    #if !DEBUG_MODE
    GetRNGstate();
    #endif
    for(i = 0; i < numberofsamples; i++)
    {
        //resample the values and calc t(v*) = MDM(v') - MDM(v*)
        blockwiseboot(samples, v);
        mdm = normabsmedian(samples, vect_sorted);
        t_star = nom - mdm;

        //sum up the number of t(v*) >= t0
        *(result->p) += (double)(t_star >= t_zero);
    }
    #if !DEBUG_MODE
    PutRNGstate();
    #endif

    //divide p by the number of samples, which is the maximal possible value for p
    //so p is in interval [0,1]
    *(result->p) /= ((double)numberofsamples + 1);

    destroy_int_array(samples);
}
