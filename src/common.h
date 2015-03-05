/* 
 * File:   common.h
 * Author: stefan
 *
 * Created on 7. Juli 2010, 15:29
 */

#ifndef _COMMON_H
#define	_COMMON_H

#ifdef	__cplusplus
extern "C" {
#endif
    #define DEBUG_MODE 0
    #define SAVE_BESSEL_VALUES 0
    #define BESSEL_FILE "/home/stefan/Praktikum_Neuro/bessel.txt"
    #define MIN(a,b) ((a)<(b)?(a):(b))
    #define MAX(a,b) ((a)>(b)?(a):(b))
    #define CLAMP(x,a,b) (MAX(MIN((x),(b)),(a)))
    #define ROUND(x) (((x)-floor((x))>=0.5)?ceil((x)):floor((x)))
    #if DEBUG_MODE
        #define RANDOM(low,high) (genrand_real2()*((high)-(low))+(low))
        #define PRINT printf
    #else
        #define RANDOM(low,high) (runif((low),(high)))
        #define PRINT Rprintf
    #endif
    #define ROUND_DIGITS(x,digits) (ROUND((x)*pow(10.0,(digits)))/pow(10.0,(digits)))
    //#define RANDOM(low,high) (is_called_from_R?runif((low),(high)):genrand_real2()*((high)-(low))+(low))
    //#define RANDOM(low,high) (runif((low),(high)))

    typedef struct dbl_array{
        double* values;
        int length;
        int allocated_values;
    }
    dbl_array;

    typedef struct int_array{
        int* values;
        int length;
        int allocated_values;
    }
    int_array;

    typedef struct dbl_matrix{
        double** values;
        int rows;
        int cols;
        int allocated_values;
    }
    dbl_matrix;

    typedef struct int_matrix{
        int** values;
        int rows;
        int cols;
        int allocated_values;
    }
    int_matrix;

    typedef struct final_result{
        int_array* binarized_vector;
        double* p;
        double* threshold;
    }
    final_result;

    //extern int is_called_from_R;

    dbl_array* init_dbl_array(double* values, int length, int keep);
    //dbl_array* init_dbl_array_keep(double* values, int length);
    int_array* init_int_array(int* values, int length, int keep);
    //int_array* init_int_array_keep(int* values, int length);
    dbl_matrix* init_dbl_matrix(double* values, int rows, int cols, int keep);
    //dbl_matrix* init_dbl_matrix_keep(double* values, int rows, int cols);
    int_matrix* init_int_matrix(int* values, int rows, int cols, int keep);
    //int_matrix* init_int_matrix_keep(int* values, int rows, int cols);
    void cut_dbl_matrix(dbl_matrix* mat, double* values, int row_begin, int row_end, int col_begin, int col_end);
    void cut_int_matrix(int_matrix* mat, int* values, int row_begin, int row_end, int col_begin, int col_end);
    void cut_dbl_array(dbl_array* arr, double* values, int begin, int end);
    void cut_int_array(int_array* arr, int* values, int begin, int end);
    void destroy_dbl_matrix(dbl_matrix* mat);
    void destroy_int_matrix(int_matrix* mat);
    void destroy_dbl_array(dbl_array* arr);
    void destroy_int_array(int_array* arr);
    void print_dbl_array(dbl_array* arr, const char* name);
    void print_int_array(int_array* arr, const char* name);
    void print_dbl_matrix(dbl_matrix* mat, const char* name);
    void print_int_matrix(int_matrix* mat, const char* name);

    int comp(const void* a, const void* b);
    int comp_int(const void* a, const void* b);
    int comp_desc(const void* a, const void* b);
    int comp_desc_int(const void* a, const void* b);

    int indexOf(double value, dbl_array* array);
    int indexOf_int(int value, int_array* array);

    double mean(double* values, int a, int b);
    double median(int_array* vect);
    double normabsmedian(int_array* vect, dbl_array* vector_sorted);
    void blockwiseboot(int_array* result, int_array* originals);
    void calc_final_results(final_result* result, int_array* v, dbl_array* vect, dbl_array* vect_sorted, double tau, int numberofsamples);



#ifdef	__cplusplus
}
#endif

#endif	/* _COMMON_H */

