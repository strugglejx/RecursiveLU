#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define BLOCK_SIZE 8
#define BLOCK_SIZE_ROW 4
#define BLOCK_SIZE_COL 64

/*
        BL_MUL_TR做汇编优化（只做一重分块）
*/

/*
        使用LU分解进行矩阵求逆（8*8正方块改进为 dn * dm（如8*16,8*32））
*/

/*
    测试不同大分块大小（d）对性能的影响
    #define BLOCK_SIZE 8
    #define BLOCK_SIZE_ROW 4
    #define BLOCK_SIZE_COL 16

*/

int __TEST5_org_inverse_matrix_n_Square_10(double *src, double *dst, int n0);

/*********************SQUARE BLOCKING***************************/

inline int LU_Decomposition_BR_10(double *tri, int n, int m, int d)
{
  /*
        对U22进行Crout分解，并保存在tri中（U22的大小小于等于d）
  */

  // int i, j, k, id;
  // for (i = 0; i < m; i++) {
  //   for (j = 0; j < m; j++) {
  //     // int in = i * n;
  //     if (j <= i) { // 下三角
  //       for (k = 0; k < j; k++) {
  //         tri[i * n + j] -= tri[i * n + k] * tri[k * n + j];
  //       }
  //     } else { // 上三角
  //       for (k = 0; k < i; k++) {
  //         tri[i * n + j] -= tri[i * n + k] * tri[k * n + j];
  //       }
  //       tri[i * n + j] /= tri[i * n + i];
  //     }
  //   }
  // }

  // 指针优化
  int i, j, k, id;
  double *tri_in_j;
  double *tri_in;
  double *tri_kn;
  for (j = 0; j < m; j++)
  {
    tri_in_j = tri + j;
    tri_in = tri;
    for (i = 0; i < m; i++, tri_in_j += n, tri_in += n)
    {
      tri_kn = tri + j;
      if (j <= i)
      { // 下三角
        for (k = 0; k < j; k++, tri_kn += n)
        {
          *tri_in_j -= *(tri_in + k) * *(tri_kn);
        }
      }
      else
      { // 上三角
        for (k = 0; k < i; k++, tri_kn += n)
        {
          *tri_in_j -= *(tri_in + k) * *(tri_kn);
        }
        *tri_in_j /= *(tri_in + i);
      }
    }
  }
}

inline int LU_Decomposition_TR_Square_10(double *dst, double *tri, int n, int m,
                                         int d)
{
  /*
      计算出U12/TR（上三角），并保存在tri中（按列一定比按行快）
  */

  // 分块；指针化
  int i, j, jj, k, id;
  double *tri_in;
  double *tri_in_k;
  double *tri_kn;
  for (j = d; j < m; j += BLOCK_SIZE)
  {
    tri_in = tri;
    for (i = 0; i < d; i++, tri_in += n)
    {
      int min_j = (BLOCK_SIZE + j < m) ? (BLOCK_SIZE + j) : (m);
      tri_kn = tri;
      tri_in_k = tri_in;
      for (k = 0; k < i; k++, tri_kn += n, tri_in_k++)
      {
        for (jj = j; jj < min_j; jj++)
        {
          *(tri_in + jj) -= *tri_in_k * *(tri_kn + jj);
        }
      }
      for (jj = j; jj < min_j; jj++)
      {
        *(tri_in + jj) /= *tri_in_k;
      }
    }
  }
}

inline int LU_L21_mul_U12_Rectangle_10(double *dst, double *tri, int n, int m,
                                       int d)
{
  // 长方块

  // 尽量避免流水线中的数据冲突

  // 分块内联汇编（8*8）（默认d为8）

  // (1*8) * (8*8) = (1*8) --- 按列循环
  // tri[i, 8:16] -= tri[i, 0:8] * tri[0:8, 0:8]
  /*
    v24   v25   v26   v27   v28     v29     v30     v31
    */

  /*
    v8   v9   v10  v11
    v12  v13  v14  v15
    v16  v17  v18  v19
    v20  v21  v22  v23
    v8   v9   v10  v11
    v12  v13  v14  v15
    v16  v17  v18  v19
    v20  v21  v22  v23
    */

  /*
    v4   v5   v6   v7
    */

  int i, j, jj;
  double *tri_in = tri + d * n;
  for (i = d; i < m; i += BLOCK_SIZE_ROW, tri_in += n * BLOCK_SIZE_ROW)
  {
    for (j = d; j < m; j += BLOCK_SIZE_COL)
    {
      int min_j = (BLOCK_SIZE_COL + j < m) ? (BLOCK_SIZE_COL + j) : (m);
      int clip_col = (m - j < BLOCK_SIZE_COL) ? (m - j) : (BLOCK_SIZE_COL);
      int tmp = 0;

      __asm__ __volatile(
          "mov x9, %[tri]\n"    // 指针
          "lsl x15, %[j], #3\n" // 行偏移量

          "mov x11, %[tri_in]\n"             // BL指针 tri_in
          "lsl x10, %[cnt], #3\n"            // 每行偏移量为n * 8，共偏移8行
          "add x13, x9, x15\n"               // TR指针:tri + j * 8
          "add x12, x11, x15\n"              // BR指针: tri_in + j * 8
          "mov x14, %[block_size_row]\n"     // 行计数器
          "add x16, x11, #32\n"              // BL指针+4


          "1:"
          // load tri[i, 0:8]
          "ld4r {v24.2d-v27.2d}, [x11]\n" // 读取连续的8个double数据（tri[i,
                                          // 0:8]）
          // load tri[i, 0:8]
          "ld4r {v28.2d-v31.2d}, [x16]\n" // 读取连续的8个double数据（tri[i,
                                          // 0:8]）

          "lsr x17, %[block_size_col], #3\n" // 列计数器
          "2:"
          // load tri[i, 8:16]
          "ld1 {v4.2d-v7.2d}, [x12]\n" // 读取连续的8个double数据（tri[i,
                                       // 8:16]）

          // // TR第1-4行
          "ld1 {v8.2d-v11.2d}, [x13]\n"  // 读取连续的8个double数据
          "add x13, x13, x10\n"          // 指针偏移
          "ld1 {v12.2d-v15.2d}, [x13]\n" // 读取连续的8个double数据
          "add x13, x13, x10\n"          // 指针偏移
          "ld1 {v16.2d-v19.2d}, [x13]\n" // 读取连续的8个double数据
          "add x13, x13, x10\n"          // 指针偏移
          "ld1 {v20.2d-v23.2d}, [x13]\n" // 读取连续的8个double数据
          "add x13, x13, x10\n"          // 指针偏移

          "fmls v4.2d, v24.2d, v8.2d\n"  // 乘减
          "fmls v5.2d, v24.2d, v9.2d\n"  // 乘减
          "fmls v6.2d, v24.2d, v10.2d\n" // 乘减
          "fmls v7.2d, v24.2d, v11.2d\n" // 乘减

          "fmls v4.2d, v25.2d, v12.2d\n" // 乘减
          "fmls v5.2d, v25.2d, v13.2d\n" // 乘减
          "fmls v6.2d, v25.2d, v14.2d\n" // 乘减
          "fmls v7.2d, v25.2d, v15.2d\n" // 乘减

          "fmls v4.2d, v26.2d, v16.2d\n" // 乘减
          "fmls v5.2d, v26.2d, v17.2d\n" // 乘减
          "fmls v6.2d, v26.2d, v18.2d\n" // 乘减
          "fmls v7.2d, v26.2d, v19.2d\n" // 乘减

          "fmls v4.2d, v27.2d, v20.2d\n" // 乘减
          "fmls v5.2d, v27.2d, v21.2d\n" // 乘减
          "fmls v6.2d, v27.2d, v22.2d\n" // 乘减
          "fmls v7.2d, v27.2d, v23.2d\n" // 乘减

          "ld1 {v8.2d-v11.2d}, [x13]\n"  // 读取连续的8个double数据
          "add x13, x13, x10\n"          // 指针偏移
          "ld1 {v12.2d-v15.2d}, [x13]\n" // 读取连续的8个double数据
          "add x13, x13, x10\n"          // 指针偏移
          "ld1 {v16.2d-v19.2d}, [x13]\n" // 读取连续的8个double数据
          "add x13, x13, x10\n"          // 指针偏移
          "ld1 {v20.2d-v23.2d}, [x13]\n" // 读取连续的8个double数据
          "add x13, x13, x10\n"          // 指针偏移

          "fmls v4.2d, v28.2d, v8.2d\n"  // 乘减
          "fmls v5.2d, v28.2d, v9.2d\n"  // 乘减
          "fmls v6.2d, v28.2d, v10.2d\n" // 乘减
          "fmls v7.2d, v28.2d, v11.2d\n" // 乘减

          "fmls v4.2d, v29.2d, v12.2d\n" // 乘减
          "fmls v5.2d, v29.2d, v13.2d\n" // 乘减
          "fmls v6.2d, v29.2d, v14.2d\n" // 乘减
          "fmls v7.2d, v29.2d, v15.2d\n" // 乘减

          "fmls v4.2d, v30.2d, v16.2d\n" // 乘减
          "fmls v5.2d, v30.2d, v17.2d\n" // 乘减
          "fmls v6.2d, v30.2d, v18.2d\n" // 乘减
          "fmls v7.2d, v30.2d, v19.2d\n" // 乘减

          "fmls v4.2d, v31.2d, v20.2d\n" // 乘减
          "fmls v5.2d, v31.2d, v21.2d\n" // 乘减
          "fmls v6.2d, v31.2d, v22.2d\n" // 乘减
          "fmls v7.2d, v31.2d, v23.2d\n" // 乘减

          "st1 {v4.2d-v7.2d}, [x12]\n" // 保存

          "add x15, x15, #64\n" // j = j + 8
          "add x12, x12, #64\n" // BR指针向右平移8个单位
          "add x13, x9, x15\n"  // TR指针向右平移8个单位
          "sub x17, x17, #1\n"  // 列计数器（每次跳8个单位）
          "cbnz x17, 2b\n"      // 跳转按行遍历

          "lsl x15, %[j], #3\n" // 行偏移量
          "sub x14, x14, #1\n"  // 行计数器
          "add x11, x11, x10\n" // BL指针向下偏移1行
          "add x16, x16, x10\n" // BL+4指针向下偏移1行
          "add x12, x11, x15\n" // BR指针: tri_in + j * 8
          "add x13, x9, x15\n"  // TR指针:tri + j * 8

          "cbnz x14, 1b\n" // 跳转按列遍历

          : 
          : [tri] "r"(tri), [tri_in] "r"(tri_in), [cnt] "r"(n), [j] "r"(j), [block_size_row] "r"(BLOCK_SIZE_ROW), [block_size_col] "r"(clip_col)
          : "memory");

      // printf("tmp:%d\n", tmp);
    }
  }
}

inline int LU_Decomposition_L_10(double *dst, double *tri, int n, int m, int d)
{
  /*
        Crout分解（单位上三角矩阵）计算出L11,L21,U11，并保存在tri中
  */

  int i, j, k;
  double *tri_in = tri;
  double *tri_in_i = tri;
  for (i = 0; i < m; i++, tri_in += n, tri_in_i += (n + 1))
  {
    double *tri_in_j = tri_in;
    for (j = 0; j < d; j++, tri_in_j++)
    {
      double *tri_kn_j = tri + j;
      double *tri_in_k = tri_in;
      if (j <= i)
      { // 下三角
        for (k = 0; k < j; k++, tri_kn_j += n, tri_in_k++)
        {
          *(tri_in_j) -= *(tri_in_k) * *tri_kn_j; // tri[k * n + j];
        }
      }
      else
      { // 上三角
        for (k = 0; k < i; k++, tri_kn_j += n, tri_in_k++)
        {
          *(tri_in_j) -= *(tri_in_k) * *tri_kn_j;
        }
        *(tri_in_j) /= *(tri_in_i);
      }
    }
  }
}

int __TEST5_org_inverse_matrix_n_Rectangle_10(double *src, double *dst, int n)
{
  /*
        LU分解的递归分块算法
  */

  int i, j, k, id;
  // int d = 2;
  // int loop_count = 3;
  // int d = 4;
  // int loop_count = 2;
  int d = 8;
  int loop_count = n/d-1; //n / d - 1;
  // int d = 16;
  // int loop_count = 31;
  // int d = 32;
  // int loop_count = 15;

  memset(dst, 0, sizeof(double) * n * n);

  // 上三角矩阵初始化
  // dst初始化为单位矩阵
  double *dst_ii = dst;

  for (int i = 0; i < n * n; i++)
  {
    dst[i] = src[i];
  }

  for (int i = 0; i < loop_count; i++)
  {
    int temp1 = (n + 1) * d * i;
    int temp2 = n - d * i;
    LU_Decomposition_L_10(dst + temp1, dst + temp1, n, temp2, d);
    LU_Decomposition_TR_Square_10(dst + temp1, dst + temp1, n, temp2, d);
    LU_L21_mul_U12_Rectangle_10(dst + temp1, dst + temp1, n, temp2, d);
  }

  LU_Decomposition_BR_10(dst + loop_count * d * (n + 1), n, n - loop_count * d, d);
}
