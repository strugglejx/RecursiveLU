
#ifndef __TEST_H_
#define __TEST_H_

#include <sys/time.h>
#include <time.h>

struct algorithm_struct
{
	char name[128];
	int (*test)(double *src, double *dst, int n0);
	int cnt;
	struct timeval start, end;
	// struct sum;
	long sum;
	double error;
};

int order_arr[6] = {32,64,128,256,512,1024};

/****************************** ./cal/LU_C.c **************************************/

extern int __TEST3_org_inverse_matrix_n(double *src, double *dst, int n0);
extern int __TEST3_opt1_inverse_matrix_n(double *src, double *dst, int n0);
extern int __TEST3_opt20_inverse_matrix_n(double *src, double *dst, int n0);
extern int __TEST3_opt21_inverse_matrix_n(double *src, double *dst, int n0);
extern int __TEST3_opt22_inverse_matrix_n(double *src, double *dst, int n0);
extern int __TEST3_opt31_inverse_matrix_n(double *src, double *dst, int n0);
extern int __TEST3_opt3_inverse_matrix_n(double *src, double *dst, int n0);
extern int __TEST3_opt40_inverse_matrix_n(double *src, double *dst, int n0);
extern int __TEST3_opt41_inverse_matrix_n(double *src, double *dst, int n0);
extern int __TEST3_opt42_inverse_matrix_n(double *src, double *dst, int n0);
extern int __TEST3_opt5_f_inverse_matrix_n(double *src, double *dst, int n0);
extern int __TEST3_opt6_inverse_matrix_n(double *src, double *dst, int n0);
extern int __TEST6_org_inverse_matrix_n(double *src, double *dst, int n0);

/****************************** ./cal/LU_C_Optimize.c **************************************/

extern int __TEST5_org_inverse_matrix_n(double *src, double *dst, int n);
extern int __TEST5_org_inverse_matrix_opt1_n(double *src, double *dst, int n);
extern int __TEST5_org_inverse_matrix_n_Square(double *src, double *dst, int n);
extern int __TEST5_org_inverse_matrix_n_Rectangle1(double *src, double *dst, int n);
extern int __TEST5_org_inverse_matrix_n_Rectangle2(double *src, double *dst, int n);

extern int LU_Decomposition_BR(double *tri, int n, int m, int d);
extern int LU_L21_mul_U12_Square(double *dst, double *tri, int n, int m, int d);
extern int LU_Decomposition_TR_Square(double *dst, double *tri, int n, int m, int d);
extern int LU_Decomposition_L(double *dst, double *tri, int n, int m, int d);

extern int LU_L21_mul_U12(double *dst, double *tri, int n, int m, int d);
extern int LU_Decomposition_TR(double *dst, double *tri, int n, int m, int d);

extern int LU_L21_mul_U12_Rectangle1(double *dst, double *tri, int n, int m, int d);
extern int LU_L21_mul_U12_Rectangle2(double *dst, double *tri, int n, int m, int d);
extern int LU_Decomposition_TR_Rectangle(double *dst, double *tri, int n, int m, int d);

extern int LU_L21_mul_U12_opt1(double *dst, double *tri, int n, int m, int d);

/****************************** ./cal/LU_C_ASM_1.c **************************************/

extern int __TEST5_org_inverse_matrix_n_Rectangle_8(double *src, double *dst, int n);

extern int LU_Decomposition_BR_8(double *tri, int n, int m, int d);
extern int LU_L21_mul_U12_Rectangle_8(double *tri, int n, int m, int d);
extern int LU_Decomposition_TR_Square_8(double *tri, int n, int m, int d);
extern int LU_Decomposition_L_8(double *tri, int n, int m, int d);

/****************************** ./cal/LU_C_ASM_2.c **************************************/

extern int __TEST5_org_inverse_matrix_n_Rectangle_9(double *src, double *dst, int n);

extern int LU_Decomposition_BR_9(double *tri, int n, int m, int d);
extern int LU_L21_mul_U12_Rectangle_9(double *dst, double *tri, int n, int m, int d);
extern int LU_Decomposition_TR_Square_9(double *dst, double *tri, int n, int m, int d);
extern int LU_Decomposition_L_9(double *dst, double *tri, int n, int m, int d);

/****************************** ./cal/LU_C_ASM_3.c **************************************/

extern int __TEST5_org_inverse_matrix_n_Rectangle_10(double *src, double *dst, int n);

extern int LU_Decomposition_BR_10(double *tri, int n, int m, int d);
extern int LU_L21_mul_U12_Rectangle_10(double *dst, double *tri, int n, int m, int d);
extern int LU_Decomposition_TR_Square_10(double *dst, double *tri, int n, int m, int d);
extern int LU_Decomposition_L_10(double *dst, double *tri, int n, int m, int d);

#define INDEX 32
struct algorithm_struct test[INDEX] = {

	{
		.name = "__TEST3_opt42_",
		.test = __TEST3_opt42_inverse_matrix_n,
		.cnt = 0,
		.sum = 0,
		.error = 0,
	},
	{
		.name = "__TEST_5_ORIGIN_4_",
		.test = __TEST5_org_inverse_matrix_n,
		.cnt = 0,
		.sum = 0,
		.error = 0,
	},
	{
		.name = "__TEST_5_ORIGIN_4_OPT1_",
		.test = __TEST5_org_inverse_matrix_opt1_n,
		.cnt = 0,
		.sum = 0,
		.error = 0,
	},
	{
		.name = "__TEST_5_RECTANGLE1_4_",
		.test = __TEST5_org_inverse_matrix_n_Rectangle1,
		.cnt = 0,
		.sum = 0,
		.error = 0,
	},
	{
		.name = "__TEST_5_RECTANGLE_9_BL_MUL_TR_(SERIAL)",
		.test = __TEST5_org_inverse_matrix_n_Rectangle_9,
		.cnt = 0,
		.sum = 0,
		.error = 0,
	},
	{
		.name = "__TEST_5_RECTANGLE_10_BL_MUL_TR_(SINGLE_BLOCK)",
		.test = __TEST5_org_inverse_matrix_n_Rectangle_10,
		.cnt = 0,
		.sum = 0,
		.error = 0,
	},

};

#define for_each_test(pt) for (pt = test; pt->test != NULL; pt++)

int ORDER;
int COUNT;

#endif /*__TEST_H_*/
