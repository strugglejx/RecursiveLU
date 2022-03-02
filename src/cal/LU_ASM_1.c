#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define BLOCK_SIZE 8
#define BLOCK_SIZE_ROW 4
#define BLOCK_SIZE_COL 64
/*
        BL_MUL_TR、TR、BL三个模块都做汇编优化
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

int __TEST5_org_inverse_matrix_n_Rectangle_8(double *src, double *dst, int n0);

/*********************SQUARE BLOCKING***************************/

inline int LU_Decomposition_BR_8(double *tri, int n, int m, int d)
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

inline int LU_L21_mul_U12_Rectangle_8(double *tri, int n, int m,
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
      for (jj = j; jj < min_j; jj += 8)
      {
        __asm__ __volatile(

            "mov x9, %[tri]\n"         // 指针
            "lsl x15, %[j], #3\n"      // 行偏移量
            "mov x11, %[tri_in]\n"     // BL指针 tri_in
            "lsl x10, %[cnt], #3\n"    // 每行偏移量为n * 8，共偏移8行
            "add x13, x9, x15\n"       // TR指针:tri + j * 8
            "add x12, x11, x15\n"      // BR指针: tri_in + j * 8
            "mov x14, %[block_size]\n" // 计数器
            "add x19, x11, #32\n"      // BL指针+4

            // // TR第1-4行
            "ld1 {v8.2d-v11.2d}, [x13]\n"  // 读取连续的8个double数据
            "add x13, x13, x10\n"          // 指针偏移
            "ld1 {v12.2d-v15.2d}, [x13]\n" // 读取连续的8个double数据
            "add x13, x13, x10\n"          // 指针偏移
            "ld1 {v16.2d-v19.2d}, [x13]\n" // 读取连续的8个double数据
            "add x13, x13, x10\n"          // 指针偏移
            "ld1 {v20.2d-v23.2d}, [x13]\n" // 读取连续的8个double数据
            "add x13, x13, x10\n"          // 指针偏移

            "1:"
            // load tri[i, 0:8]
            "ld4r {v24.2d-v27.2d}, [x11]\n" // 读取连续的8个double数据（tri[i,
                                            // 0:8]）
            // load tri[i, 8:16]
            "ld1 {v4.2d-v7.2d}, [x12]\n" // 读取连续的8个double数据（tri[i,
                                         // 8:16]）

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

            "st1 {v4.2d-v7.2d}, [x12]\n" // 保存

            "sub x14, x14, #1\n"
            "add x12, x12, x10\n" // BR指针向下偏移1行
            "add x11, x11, x10\n" // BL指针向下偏移1行
            "cbnz x14, 1b\n"      // 跳转按列遍历

            /**********************************************/

            "add x12, %[tri_in], x15\n" // BR指针: tri_in + j * 8
            "mov x14, %[block_size]\n"  // 计数器

            // TR第5-8行
            "ld1 {v8.2d-v11.2d}, [x13]\n"  // 读取连续的8个double数据
            "add x13, x13, x10\n"          // 指针偏移
            "ld1 {v12.2d-v15.2d}, [x13]\n" // 读取连续的8个double数据
            "add x13, x13, x10\n"          // 指针偏移
            "ld1 {v16.2d-v19.2d}, [x13]\n" // 读取连续的8个double数据
            "add x13, x13, x10\n"          // 指针偏移
            "ld1 {v20.2d-v23.2d}, [x13]\n" // 读取连续的8个double数据
            "add x13, x13, x10\n"          // 指针偏移

            "2:"
            // load tri[i, 0:8]
            "ld4r {v28.2d-v31.2d}, [x19]\n" // 读取连续的8个double数据（tri[i,
                                            // 0:8]）

            // load tri[i, 8:16]
            "ld1 {v4.2d-v7.2d}, [x12]\n" // 读取连续的8个double数据（tri[i,
                                         // 8:16]）

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

            "sub x14, x14, #1\n"
            "add x12, x12, x10\n" // BR指针向下偏移1行
            "add x19, x19, x10\n" // BL+4指针向下偏移1行
            "cbnz x14, 2b\n"      // 跳转按列遍历

            :
            : [tri] "r"(tri), [tri_in] "r"(tri_in), [cnt] "r"(n), [j] "r"(jj), [block_size] "r"(BLOCK_SIZE_ROW)
            : "memory");
      }
    }
  }
}

inline int LU_Decomposition_TR_Square_8(double *tri, int n, int m,
                                        int d)
{
  /*
      计算出U12/TR（上三角），并保存在tri中（按列一定比按行快）
  */

  int m_d = m - d;
  if (m_d <= 0)
  {
    return 0;
  }

  // 分块；指针化
  // int i, j, jj, k, id;
  // double *tri_in;
  // double *tri_in_k;
  // double *tri_kn;
  // for (j = d; j < m; j += BLOCK_SIZE)
  // {
  //   tri_in = tri;
  //   for (i = 0; i < d; i++, tri_in += n)
  //   {
  //     int min_j = (BLOCK_SIZE + j < m) ? (BLOCK_SIZE + j) : (m);
  //     tri_kn = tri;
  //     tri_in_k = tri_in;
  //     for (k = 0; k < i; k++, tri_kn += n, tri_in_k++)
  //     {
  //       for (jj = j; jj < min_j; jj++)
  //       {
  //         *(tri_in + jj) -= *tri_in_k * *(tri_kn + jj);
  //       }
  //     }
  //     for (jj = j; jj < min_j; jj++)
  //     {
  //       *(tri_in + jj) /= *tri_in_k;
  //     }
  //   }
  // }
  
  int j;
  for (j = d; j < d+m; j += 8)
  {
    double *tri_j = tri + j;
    __asm__ __volatile(

  "mov x9, %[tri]\n"      // 指针
  "mov x12, %[tri]\n"     // TL指针
  "lsl x10, %[cnt], #3\n" // 每行偏移量为n * 8

  "1:"
  "mov x13, %[tri_j]\n" // TR指针

  // 第一行
  "ld1r {v0.2d}, [x12]\n"       // 读取连续的2个double数据（tri[0,
                                // 0:1]）
  "ld1 {v8.2d-v11.2d}, [x13]\n" // 读取连续的8个double数据

  "add x12, x12, x10\n"

  "fdiv v8.2d, v8.2d, v0.2d\n"   // 除
  "fdiv v9.2d, v9.2d, v0.2d\n"   // 除
  "fdiv v10.2d, v10.2d, v0.2d\n" // 除
  "fdiv v11.2d, v11.2d, v0.2d\n" // 除

  "st1 {v8.2d-v11.2d}, [x13]\n" // 保存在内存中
  "add x13, x13, x10\n"         // 指针偏移

  // 第二行
  "ld2r {v0.2d-v1.2d}, [x12]\n"  // 读取连续的2个double数据（tri[1,
                                 // 0:1]）
  "ld1 {v12.2d-v15.2d}, [x13]\n" // 读取连续的8个double数据

  "add x12, x12, x10\n"

  "fmls v12.2d, v0.2d, v8.2d\n"  // 乘减
  "fmls v13.2d, v0.2d, v9.2d\n"  // 乘减
  "fmls v14.2d, v0.2d, v10.2d\n" // 乘减
  "fmls v15.2d, v0.2d, v11.2d\n" // 乘减

  "fdiv v12.2d, v12.2d, v1.2d\n" // 除
  "fdiv v13.2d, v13.2d, v1.2d\n" // 除
  "fdiv v14.2d, v14.2d, v1.2d\n" // 除
  "fdiv v15.2d, v15.2d, v1.2d\n" // 除

  "st1 {v12.2d-v15.2d}, [x13]\n" // 保存在内存中
  "add x13, x13, x10\n"          // 指针偏移

  // 第三行
  "ld3r {v0.2d-v2.2d}, [x12]\n"  // 读取连续的4个double数据（tri[1,
                                 // 0:3]）
  "ld1 {v16.2d-v19.2d}, [x13]\n" // 读取连续的8个double数据

  "add x12, x12, x10\n"

  "fmls v16.2d, v0.2d, v8.2d\n"  // 乘减
  "fmls v17.2d, v0.2d, v9.2d\n"  // 乘减
  "fmls v18.2d, v0.2d, v10.2d\n" // 乘减
  "fmls v19.2d, v0.2d, v11.2d\n" // 乘减

  "fmls v16.2d, v1.2d, v12.2d\n" // 乘减
  "fmls v17.2d, v1.2d, v13.2d\n" // 乘减
  "fmls v18.2d, v1.2d, v14.2d\n" // 乘减
  "fmls v19.2d, v1.2d, v15.2d\n" // 乘减

  "fdiv v16.2d, v16.2d, v2.2d\n" // 除
  "fdiv v17.2d, v17.2d, v2.2d\n" // 除
  "fdiv v18.2d, v18.2d, v2.2d\n" // 除
  "fdiv v19.2d, v19.2d, v2.2d\n" // 除

  "st1 {v16.2d-v19.2d}, [x13]\n" // 保存在内存中
  "add x13, x13, x10\n"          // 指针偏移

  // 第四行
  "ld4r {v0.2d-v3.2d}, [x12]\n"  // 读取连续的4个double数据（tri[1,
                                 // 0:3]）
  "ld1 {v20.2d-v23.2d}, [x13]\n" // 读取连续的8个double数据

  "add x12, x12, x10\n"

  "fmls v20.2d, v0.2d, v8.2d\n"  // 乘减
  "fmls v21.2d, v0.2d, v9.2d\n"  // 乘减
  "fmls v22.2d, v0.2d, v10.2d\n" // 乘减
  "fmls v23.2d, v0.2d, v11.2d\n" // 乘减

  "fmls v20.2d, v1.2d, v12.2d\n" // 乘减
  "fmls v21.2d, v1.2d, v13.2d\n" // 乘减
  "fmls v22.2d, v1.2d, v14.2d\n" // 乘减
  "fmls v23.2d, v1.2d, v15.2d\n" // 乘减

  "fmls v20.2d, v2.2d, v16.2d\n" // 乘减
  "fmls v21.2d, v2.2d, v17.2d\n" // 乘减
  "fmls v22.2d, v2.2d, v18.2d\n" // 乘减
  "fmls v23.2d, v2.2d, v19.2d\n" // 乘减

  "fdiv v20.2d, v20.2d, v3.2d\n" // 除
  "fdiv v21.2d, v21.2d, v3.2d\n" // 除
  "fdiv v22.2d, v22.2d, v3.2d\n" // 除
  "fdiv v23.2d, v23.2d, v3.2d\n" // 除

  "st1 {v20.2d-v23.2d}, [x13]\n" // 保存在内存中
  "add x13, x13, x10\n"          // 指针偏移

  // 第五行
  "add x14, x12, #32\n"
  "ld4r {v0.2d-v3.2d}, [x12]\n"  // 读取连续的6个double数据（tri[1,
                                 // 0:3]）
  "ld1 {v24.2d-v27.2d}, [x13]\n" // 读取连续的8个double数据

  "fmls v24.2d, v0.2d, v8.2d\n"  // 乘减
  "fmls v25.2d, v0.2d, v9.2d\n"  // 乘减
  "fmls v26.2d, v0.2d, v10.2d\n" // 乘减
  "fmls v27.2d, v0.2d, v11.2d\n" // 乘减

  "ld1r {v4.2d}, [x14]\n"

  "fmls v24.2d, v1.2d, v12.2d\n" // 乘减
  "fmls v25.2d, v1.2d, v13.2d\n" // 乘减
  "fmls v26.2d, v1.2d, v14.2d\n" // 乘减
  "fmls v27.2d, v1.2d, v15.2d\n" // 乘减

  "fmls v24.2d, v2.2d, v16.2d\n" // 乘减
  "fmls v25.2d, v2.2d, v17.2d\n" // 乘减
  "fmls v26.2d, v2.2d, v18.2d\n" // 乘减
  "fmls v27.2d, v2.2d, v19.2d\n" // 乘减

  "add x12, x12, x10\n"
  "add x14, x14, x10\n"

  "fmls v24.2d, v3.2d, v20.2d\n" // 乘减
  "fmls v25.2d, v3.2d, v21.2d\n" // 乘减
  "fmls v26.2d, v3.2d, v22.2d\n" // 乘减
  "fmls v27.2d, v3.2d, v23.2d\n" // 乘减

  "fdiv v24.2d, v24.2d, v4.2d\n" // 除
  "fdiv v25.2d, v25.2d, v4.2d\n" // 除
  "fdiv v26.2d, v26.2d, v4.2d\n" // 除
  "fdiv v27.2d, v27.2d, v4.2d\n" // 除

  "st1 {v24.2d-v27.2d}, [x13]\n" // 保存在内存中
  "add x13, x13, x10\n"          // 指针偏移

  // 第六行
  "ld4r {v0.2d-v3.2d}, [x12]\n"  // 读取连续的6个double数据（tri[1,
                                 // 0:3]）
  "ld1 {v28.2d-v31.2d}, [x13]\n" // 读取连续的8个double数据

  "fmls v28.2d, v0.2d, v8.2d\n"  // 乘减
  "fmls v29.2d, v0.2d, v9.2d\n"  // 乘减
  "fmls v30.2d, v0.2d, v10.2d\n" // 乘减
  "fmls v31.2d, v0.2d, v11.2d\n" // 乘减

  "fmls v28.2d, v1.2d, v12.2d\n" // 乘减
  "fmls v29.2d, v1.2d, v13.2d\n" // 乘减
  "fmls v30.2d, v1.2d, v14.2d\n" // 乘减
  "fmls v31.2d, v1.2d, v15.2d\n" // 乘减

  "ld2r {v4.2d-v5.2d}, [x14]\n"

  "fmls v28.2d, v2.2d, v16.2d\n" // 乘减
  "fmls v29.2d, v2.2d, v17.2d\n" // 乘减
  "fmls v30.2d, v2.2d, v18.2d\n" // 乘减
  "fmls v31.2d, v2.2d, v19.2d\n" // 乘减

  "fmls v28.2d, v3.2d, v20.2d\n" // 乘减
  "fmls v29.2d, v3.2d, v21.2d\n" // 乘减
  "fmls v30.2d, v3.2d, v22.2d\n" // 乘减
  "fmls v31.2d, v3.2d, v23.2d\n" // 乘减

  "add x12, x12, x10\n"
  "add x14, x14, x10\n"

  "fmls v28.2d, v4.2d, v24.2d\n" // 乘减
  "fmls v29.2d, v4.2d, v25.2d\n" // 乘减
  "fmls v30.2d, v4.2d, v26.2d\n" // 乘减
  "fmls v31.2d, v4.2d, v27.2d\n" // 乘减

  "fdiv v28.2d, v28.2d, v5.2d\n" // 除
  "fdiv v29.2d, v29.2d, v5.2d\n" // 除
  "fdiv v30.2d, v30.2d, v5.2d\n" // 除
  "fdiv v31.2d, v31.2d, v5.2d\n" // 除

  "st1 {v28.2d-v31.2d}, [x13]\n" // 保存在内存中
  "add x13, x13, x10\n"

  // 第七行
  "ld4r {v0.2d-v3.2d}, [x12]\n" // 读取连续的6个double数据（tri[1,
                                // 0:3]）
  "ld3r {v4.2d-v6.2d}, [x14]\n"
  "ld1 {v28.2d-v31.2d}, [x13]\n" // 读取连续的8个double数据

  "sub x13, x13, x10\n"
  "add x12, x12, x10\n"
  "add x14, x14, x10\n"

  "fmls v28.2d, v0.2d, v8.2d\n"  // 乘减
  "fmls v29.2d, v0.2d, v9.2d\n"  // 乘减
  "fmls v30.2d, v0.2d, v10.2d\n" // 乘减
  "fmls v31.2d, v0.2d, v11.2d\n" // 乘减

  "fmls v28.2d, v1.2d, v12.2d\n" // 乘减
  "fmls v29.2d, v1.2d, v13.2d\n" // 乘减
  "fmls v30.2d, v1.2d, v14.2d\n" // 乘减
  "fmls v31.2d, v1.2d, v15.2d\n" // 乘减

  "ld1 {v8.2d-v11.2d}, [x13]\n" // 读取TR第6行连续的8个double数据
  "add x13, x13, x10\n"         // 指针偏移

  "fmls v28.2d, v2.2d, v16.2d\n" // 乘减
  "fmls v29.2d, v2.2d, v17.2d\n" // 乘减
  "fmls v30.2d, v2.2d, v18.2d\n" // 乘减
  "fmls v31.2d, v2.2d, v19.2d\n" // 乘减

  "fmls v28.2d, v3.2d, v20.2d\n" // 乘减
  "fmls v29.2d, v3.2d, v21.2d\n" // 乘减
  "fmls v30.2d, v3.2d, v22.2d\n" // 乘减
  "fmls v31.2d, v3.2d, v23.2d\n" // 乘减

  "fmls v28.2d, v4.2d, v24.2d\n" // 乘减
  "fmls v29.2d, v4.2d, v25.2d\n" // 乘减
  "fmls v30.2d, v4.2d, v26.2d\n" // 乘减
  "fmls v31.2d, v4.2d, v27.2d\n" // 乘减

  "fmls v28.2d, v5.2d, v8.2d\n"  // 乘减
  "fmls v29.2d, v5.2d, v9.2d\n"  // 乘减
  "fmls v30.2d, v5.2d, v10.2d\n" // 乘减
  "fmls v31.2d, v5.2d, v11.2d\n" // 乘减

  "fdiv v28.2d, v28.2d, v6.2d\n" // 除
  "fdiv v29.2d, v29.2d, v6.2d\n" // 除
  "fdiv v30.2d, v30.2d, v6.2d\n" // 除
  "fdiv v31.2d, v31.2d, v6.2d\n" // 除

  "st1 {v28.2d-v31.2d}, [x13]\n" // 保存在内存中
  "add x13, x13, x10\n"          // 指针偏移

  // 第八行
  "ld4r {v0.2d-v3.2d}, [x12]\n"  // 读取连续的6个double数据（tri[1,
                                 // 0:3]）
  "ld1 {v28.2d-v31.2d}, [x13]\n" // 读取TR第8行连续的8个double数据
  "sub x13, x13, x10\n"          // 指针偏移

  "fmls v28.2d, v1.2d, v12.2d\n" // 乘减
  "fmls v29.2d, v1.2d, v13.2d\n" // 乘减
  "fmls v30.2d, v1.2d, v14.2d\n" // 乘减
  "fmls v31.2d, v1.2d, v15.2d\n" // 乘减

  "ld4r {v4.2d-v7.2d}, [x14]\n"

  "fmls v28.2d, v2.2d, v16.2d\n" // 乘减
  "fmls v29.2d, v2.2d, v17.2d\n" // 乘减
  "fmls v30.2d, v2.2d, v18.2d\n" // 乘减
  "fmls v31.2d, v2.2d, v19.2d\n" // 乘减

  // 第6行
  "fmls v28.2d, v5.2d, v8.2d\n"  // 乘减
  "fmls v29.2d, v5.2d, v9.2d\n"  // 乘减
  "fmls v30.2d, v5.2d, v10.2d\n" // 乘减
  "fmls v31.2d, v5.2d, v11.2d\n" // 乘减

  "fmls v28.2d, v3.2d, v20.2d\n" // 乘减
  "fmls v29.2d, v3.2d, v21.2d\n" // 乘减
  "fmls v30.2d, v3.2d, v22.2d\n" // 乘减
  "fmls v31.2d, v3.2d, v23.2d\n" // 乘减

  "ld1 {v8.2d-v11.2d}, [%[tri_j]]\n" // 读取TR第1行连续的8个double数据
  "ld1 {v12.2d-v15.2d}, [x13]\n"     // 读取TR第7行连续的8个double数据
  "add x13, x13, x10\n"              // 指针偏移

  "fmls v28.2d, v4.2d, v24.2d\n" // 乘减
  "fmls v29.2d, v4.2d, v25.2d\n" // 乘减
  "fmls v30.2d, v4.2d, v26.2d\n" // 乘减
  "fmls v31.2d, v4.2d, v27.2d\n" // 乘减

  "fmls v28.2d, v0.2d, v8.2d\n"  // 乘减
  "fmls v29.2d, v0.2d, v9.2d\n"  // 乘减
  "fmls v30.2d, v0.2d, v10.2d\n" // 乘减
  "fmls v31.2d, v0.2d, v11.2d\n" // 乘减

  "fmls v28.2d, v6.2d, v12.2d\n" // 乘减
  "fmls v29.2d, v6.2d, v13.2d\n" // 乘减
  "fmls v30.2d, v6.2d, v14.2d\n" // 乘减
  "fmls v31.2d, v6.2d, v15.2d\n" // 乘减

  "fdiv v28.2d, v28.2d, v7.2d\n" // 除
  "fdiv v29.2d, v29.2d, v7.2d\n" // 除
  "fdiv v30.2d, v30.2d, v7.2d\n" // 除
  "fdiv v31.2d, v31.2d, v7.2d\n" // 除

  "st1 {v28.2d-v31.2d}, [x13]\n" // 保存在内存中

        :
        : [tri_j] "r"(tri_j), [tri] "r"(tri), [cnt] "r"(n)
        : "memory");
  }
}

inline int LU_Decomposition_L_8(double *tri, int n, int m, int d)
{
  /*
        Crout分解（单位上三角矩阵）计算出L11,L21,U11，并保存在tri中
  */


  __asm__ __volatile(

      "mov x9, %[tri]\n" // 指针
      "mov x10, %[cnt]\n"
      "lsl x10, x10, #3\n" // 每行偏移量为n * 8

      // 第1行
      "ld1 {v4.2d-v7.2d}, [x9]\n" // 读取连续的8个double数据
      "dup v0.2d, v4.d[0]\n"      // 临时变量（除数tri[0][0]）
      "mov d1, v4.d[1]\n"
      "fdiv d1, d1, d4\n"
      "mov v4.d[1], v1.d[0]\n"
      "fdiv v5.2d, v5.2d, v0.2d\n"
      "fdiv v6.2d, v6.2d, v0.2d\n"
      "fdiv v7.2d, v7.2d, v0.2d\n"

      // 保存第1行
      "st1 {v4.2d-v7.2d}, [x9]\n"

      // 第2行前2个元素
      "add x9, x9, x10\n" // 指针增加n * 8
      "ld1 {v8.2d-v11.2d}, [x9]\n"

      "mov d0, v8.d[1]\n"
      "fmls d0, d8, v4.d[1]\n"
      "mov v8.d[1], v0.d[0]\n"

      // 第2行后6个元素
      "dup v0.2d, v8.d[0]\n"         // 临时变量（被乘数）
      "dup v1.2d, v8.d[1]\n"         // 临时变量（除数tri[1][1]）
      "fmls v9.2d, v0.2d, v5.2d\n"   // 乘减
      "fdiv v9.2d, v9.2d, v1.2d\n"   // 除
      "fmls v10.2d, v0.2d, v6.2d\n"  // 乘减
      "fdiv v10.2d, v10.2d, v1.2d\n" // 除
      "fmls v11.2d, v0.2d, v7.2d\n"  // 乘减
      "fdiv v11.2d, v11.2d, v1.2d\n" // 除

      // 保存第2行
      "st1 {v8.2d-v11.2d}, [x9]\n"

      // 第3行第2个元素
      "add x9, x9, x10\n"
      "ld1 {v12.2d-v15.2d}, [x9]\n"
      "mov d0, v12.d[1]\n"
      "fmls d0, d12, v4.d[1]\n" // 1
      "mov v12.d[1], v0.d[0]\n"

      // 第3行后6个元素（减去第1行）
      "dup v0.2d, v12.d[0]\n"       // 临时变量（被乘数tri[2][0]）
      "fmls v13.2d, v0.2d, v5.2d\n" // 乘减
      "fmls v14.2d, v0.2d, v6.2d\n" // 乘减
      "fmls v15.2d, v0.2d, v7.2d\n" // 乘减

      // 第3行后6个元素（减去第2行）
      "dup v0.2d, v12.d[1]\n"        // 临时变量（被乘数tri[2][1]
      "fmls v13.2d, v0.2d, v9.2d\n"  // 乘减
      "fmls v14.2d, v0.2d, v10.2d\n" // 乘减
      "fmls v15.2d, v0.2d, v11.2d\n" // 乘减

      // 第3行后6个元素（除以tri[2][2]）
      "dup v1.2d, v13.d[0]\n" // 临时变量（除数tri[2][2]）
      "mov d0, v13.d[1]\n"
      "fdiv d0, d0, d1\n" // 除（v13.d[0]无需除）
      "mov v13.d[1], v0.d[0]\n"
      "fdiv v14.2d, v14.2d, v1.2d\n" // 除
      "fdiv v15.2d, v15.2d, v1.2d\n" // 除

      // 保存第3行
      "st1 {v12.2d-v15.2d}, [x9]\n"

      // 第4行前2个元素
      "add x9, x9, x10\n"
      "ld1 {v16.2d-v19.2d}, [x9]\n"
      "mov d2, v16.d[1]\n"
      "fmls d2, d16, v4.d[1]\n" // 1
      "mov v16.d[1], v2.d[0]\n"

      // 第4行后6个元素（减去第1行）
      "dup v0.2d, v16.d[0]\n"       // 临时变量（被乘数tri[3][0]）
      "fmls v17.2d, v0.2d, v5.2d\n" // 乘减
      "fmls v18.2d, v0.2d, v6.2d\n" // 乘减
      "fmls v19.2d, v0.2d, v7.2d\n" // 乘减

      // 第4行后6个元素（减去第2行）
      "dup v0.2d, v16.d[1]\n"        // 临时变量（被乘数tri[3][1]）
      "fmls v17.2d, v0.2d, v9.2d\n"  // 乘减
      "fmls v18.2d, v0.2d, v10.2d\n" // 乘减
      "fmls v19.2d, v0.2d, v11.2d\n" // 乘减

      // 第4行后6个元素（减去第3行）
      "dup v0.2d, v17.d[0]\n" // 临时变量（被乘数tri[3][2]）
      "mov d2, v17.d[1]\n"
      "fmls d2, d0, v13.d[1]\n" // 乘减（v17.d[0]无需减）
      "mov v17.d[1], v2.d[0]\n"
      "fmls v18.2d, v0.2d, v14.2d\n" // 乘减
      "fmls v19.2d, v0.2d, v15.2d\n" // 乘减

      // 第4行后6个元素（除以tri[3][3]）
      "dup v1.2d, v17.d[1]\n"        // 临时变量（除数tri[3][3]）
      "fdiv v18.2d, v18.2d, v1.2d\n" // 除
      "fdiv v19.2d, v19.2d, v1.2d\n" // 除

      // 保存第4行
      "st1 {v16.2d-v19.2d}, [x9]\n"

      // 第5行前2个元素
      "add x9, x9, x10\n"
      "ld1 {v20.2d-v23.2d}, [x9]\n"
      "mov d2, v20.d[1]\n"
      "fmls d2, d20, v4.d[1]\n" // 1
      "mov v20.d[1], v2.d[0]\n"

      // 第5行后6个元素（减去第1行）
      "dup v0.2d, v20.d[0]\n"       // 临时变量（被乘数tri[4][0]）
      "fmls v21.2d, v0.2d, v5.2d\n" // 乘减
      "fmls v22.2d, v0.2d, v6.2d\n" // 乘减
      "fmls v23.2d, v0.2d, v7.2d\n" // 乘减

      // 第5行后6个元素（减去第2行）
      "dup v0.2d, v20.d[1]\n"        // 临时变量（被乘数tri[4][1]）
      "fmls v21.2d, v0.2d, v9.2d\n"  // 乘减
      "fmls v22.2d, v0.2d, v10.2d\n" // 乘减
      "fmls v23.2d, v0.2d, v11.2d\n" // 乘减

      // 第5行后6个元素（减去第3行）
      "dup v0.2d, v21.d[0]\n" // 临时变量（被乘数tri[4][2]）
      "mov d2, v21.d[1]\n"
      "fmls d2, d0, v13.d[1]\n" // 乘减（v21.d[0]无需减）
      "mov v21.d[1], v2.d[0]\n"
      "fmls v22.2d, v0.2d, v14.2d\n" // 乘减
      "fmls v23.2d, v0.2d, v15.2d\n" // 乘减

      // 第5行后6个元素（减去第4行）
      "dup v0.2d, v21.d[1]\n"        // 临时变量（被乘数tri[4][3]）
      "fmls v22.2d, v0.2d, v18.2d\n" // 乘减
      "fmls v23.2d, v0.2d, v19.2d\n" // 乘减

      // 第5行后6个元素（除以tri[4][4]）
      "dup v1.2d, v22.d[0]\n" // 临时变量（除数tri[4][4]）
      "mov d2, v22.d[1]\n"
      "fdiv d2, d2, d1\n" // 除
      "mov v22.d[1], v2.d[0]\n"
      "fdiv v23.2d, v23.2d, v1.2d\n" // 除

      // 保存第5行
      "st1 {v20.2d-v23.2d}, [x9]\n"

      // 第6行前2个元素
      "add x9, x9, x10\n"
      "ld1 {v24.2d-v27.2d}, [x9]\n"
      "mov d2, v24.d[1]\n"
      "fmls d2, d24, v4.d[1]\n" // 1
      "mov v24.d[1], v2.d[0]\n"

      // 第6行后6个元素（减去第1行）
      "dup v0.2d, v24.d[0]\n"       // 临时变量（被乘数tri[5][0]）
      "fmls v25.2d, v0.2d, v5.2d\n" // 乘减
      "fmls v26.2d, v0.2d, v6.2d\n" // 乘减
      "fmls v27.2d, v0.2d, v7.2d\n" // 乘减

      // 第6行后6个元素（减去第2行）
      "dup v0.2d, v24.d[1]\n"        // 临时变量（被乘数tri[5][1]）
      "fmls v25.2d, v0.2d, v9.2d\n"  // 乘减
      "fmls v26.2d, v0.2d, v10.2d\n" // 乘减
      "fmls v27.2d, v0.2d, v11.2d\n" // 乘减

      // 第6行后6个元素（减去第3行）
      "dup v0.2d, v25.d[0]\n" // 临时变量（被乘数tri[5][2]）
      "mov d2, v25.d[1]\n"
      "fmls d2, d0, v13.d[1]\n" // 乘减（v25.d[0]无需减）
      "mov v25.d[1], v2.d[0]\n"
      "fmls v26.2d, v0.2d, v14.2d\n" // 乘减
      "fmls v27.2d, v0.2d, v15.2d\n" // 乘减

      // 第6行后6个元素（减去第4行）
      "dup v0.2d, v25.d[1]\n"        // 临时变量（被乘数tri[5][3]）
      "fmls v26.2d, v0.2d, v18.2d\n" // 乘减
      "fmls v27.2d, v0.2d, v19.2d\n" // 乘减

      // 第6行后6个元素（减去第5行）
      "dup v0.2d, v26.d[0]\n" // 临时变量（被乘数tri[5][4]）
      "mov d2, v26.d[1]\n"
      "fmls d2, d0, v22.d[1]\n" // 乘减
      "mov v26.d[1], v2.d[0]\n"
      "fmls v27.2d, v0.2d, v23.2d\n" // 乘减

      // 第6行后6个元素（除以tri[5][5]）
      "dup v1.2d, v26.d[1]\n"        // 临时变量（除数tri[5][5]）
      "fdiv v27.2d, v27.2d, v1.2d\n" // 除

      // 保存第6行
      "st1 {v24.2d-v27.2d}, [x9]\n"

      // 第7行前2个元素
      "add x9, x9, x10\n"
      "ld1 {v28.2d-v31.2d}, [x9]\n"
      "mov d2, v28.d[1]\n"
      "fmls d2, d28, v4.d[1]\n" // 1
      "mov v28.d[1], v2.d[0]\n"

      // 第7行后6个元素（减去第1行）
      "dup v0.2d, v28.d[0]\n"       // 临时变量（被乘数tri[6][0]）
      "fmls v29.2d, v0.2d, v5.2d\n" // 乘减
      "fmls v30.2d, v0.2d, v6.2d\n" // 乘减
      "fmls v31.2d, v0.2d, v7.2d\n" // 乘减

      // 第7行后6个元素（减去第2行）
      "dup v0.2d, v28.d[1]\n"        // 临时变量（被乘数tri[6][1]）
      "fmls v29.2d, v0.2d, v9.2d\n"  // 乘减
      "fmls v30.2d, v0.2d, v10.2d\n" // 乘减
      "fmls v31.2d, v0.2d, v11.2d\n" // 乘减

      // 第7行后6个元素（减去第3行）
      "dup v0.2d, v29.d[0]\n" // 临时变量（被乘数tri[6][2]）
      "mov d2, v29.d[1]\n"
      "fmls d2, d0, v13.d[1]\n" // 乘减（v29.d[0]无需减）
      "mov v29.d[1], v2.d[0]\n"
      "fmls v30.2d, v0.2d, v14.2d\n" // 乘减
      "fmls v31.2d, v0.2d, v15.2d\n" // 乘减

      // 第7行后6个元素（减去第4行）
      "dup v0.2d, v29.d[1]\n"        // 临时变量（被乘数tri[6][3]）
      "fmls v30.2d, v0.2d, v18.2d\n" // 乘减
      "fmls v31.2d, v0.2d, v19.2d\n" // 乘减

      // 第7行后6个元素（减去第5行）
      "dup v0.2d, v30.d[0]\n" // 临时变量（被乘数tri[6][4]）
      "mov d2, v30.d[1]\n"
      "fmls d2, d0, v22.d[1]\n"      // 乘减（v30.d[0]无需减）
      "fmls v31.2d, v0.2d, v23.2d\n" // 乘减

      // 第7行后6个元素（减去第6行）
      "dup v0.2d, v2.d[0]\n"         // 临时变量（被乘数tri[6][5]）
      "fmls v31.2d, v0.2d, v27.2d\n" // 乘减

      // 第7行后6个元素（除以tri[6][6]）
      "mov v3.d[0], v31.d[1]\n"
      "fdiv d3, d3, d31\n" // 除
      "mov v30.d[1], v2.d[0]\n"
      "mov v31.d[1], v3.d[0]\n"

      // 保存第7行
      "st1 {v28.2d-v31.2d}, [x9]\n"

      // 第8行前2个元素
      "add x9, x9, x10\n"
      "ld1 {v12.2d}, [x9]\n"
      "add x9, x9, #16\n"
      "ld1 {v16.2d}, [x9]\n"
      "add x9, x9, #16\n"
      "ld1 {v20.2d}, [x9]\n"
      "add x9, x9, #16\n"
      "ld1 {v24.2d}, [x9]\n"

      "mov d2, v12.d[1]\n"
      "fmls d2, d12, v4.d[1]\n" // 1
      "mov v12.d[1], v2.d[0]\n"

      // 第8行后6个元素（减去第1行）
      "dup v0.2d, v12.d[0]\n"       // 临时变量（被乘数tri[7][0]）
      "fmls v16.2d, v0.2d, v5.2d\n" // 乘减
      "fmls v20.2d, v0.2d, v6.2d\n" // 乘减
      "fmls v24.2d, v0.2d, v7.2d\n" // 乘减

      // 第8行后6个元素（减去第2行）
      "dup v0.2d, v12.d[1]\n"        // 临时变量（被乘数tri[7][1]）
      "fmls v16.2d, v0.2d, v9.2d\n"  // 乘减
      "fmls v20.2d, v0.2d, v10.2d\n" // 乘减
      "fmls v24.2d, v0.2d, v11.2d\n" // 乘减

      // 第8行后6个元素（减去第3行）
      "dup v0.2d, v16.d[0]\n" // 临时变量（被乘数tri[7][2]）
      "mov d2, v16.d[1]\n"
      "fmls d2, d0, v13.d[1]\n" // 乘减（v16.d[0]无需减）
      "mov v16.d[1], v2.d[0]\n"
      "fmls v20.2d, v0.2d, v14.2d\n" // 乘减
      "fmls v24.2d, v0.2d, v15.2d\n" // 乘减

      // 第8行后6个元素（减去第4行）
      "dup v0.2d, v16.d[1]\n"        // 临时变量（被乘数tri[6][3]）
      "fmls v20.2d, v0.2d, v18.2d\n" // 乘减
      "fmls v24.2d, v0.2d, v19.2d\n" // 乘减

      // 第8行后6个元素（减去第5行）
      "dup v0.2d, v20.d[0]\n" // 临时变量（被乘数tri[6][4]）
      "mov d2, v20.d[1]\n"
      "fmls d2, d0, v22.d[1]\n"      // 乘减（v20.d[0]无需减）
      "fmls v24.2d, v0.2d, v23.2d\n" // 乘减

      // 第8行后6个元素（减去第6/7行）
      "dup v0.2d, v2.d[0]\n"
      "fmls v24.2d, v0.2d, v27.2d\n" // 乘减
      "mov d3, v24.d[1]\n"
      "fmls d3, d24, v31.d[1]\n" // 乘减
      "mov v20.d[1], v2.d[0]\n"
      "mov v24.d[1], v3.d[0]\n"

      // 保存第8行
      // "st1 {v12.2d,v16.2d,v20.2d,v24.2d}, [x9], #64\n"
      "st1 {v24.2d}, [x9]\n"
      "sub x9, x9, #16\n"
      "st1 {v20.2d}, [x9]\n"
      "sub x9, x9, #16\n"
      "st1 {v16.2d}, [x9]\n"
      "sub x9, x9, #16\n"
      "st1 {v12.2d}, [x9]\n"

      :
      : [tri] "r"(tri), [cnt] "r"(n)
      : "memory");

  int m_d = m - d;
  if (m_d <= 0)
  {
    return;
  }

  __asm__ __volatile(
      "mov x9, %[tri]\n"      // 指针
      "lsl x10, %[cnt], #3\n" // 每行偏移量为n * 8

      "lsl x11, %[cnt], #6\n" // n * 8 * 8
      "add x11, x9, x11\n"
      "add x12, x10, #16\n"

      // load
      "ld1 {v0.2d-v3.2d}, [x9]\n" // 读取连续的8个double数据

      "add x9, x9, x12\n"
      "ld1 {v4.2d-v6.2d}, [x9]\n" // 读取连续的6个double数据
      "add x9, x9, x10\n"
      "ld1 {v7.2d-v9.2d}, [x9]\n" // 读取连续的6个double数据

      "add x9, x9, x12\n"
      "ld1 {v10.2d-v11.2d}, [x9]\n" // 读取连续的4个double数据
      "add x9, x9, x10\n"
      "ld1 {v12.2d-v13.2d}, [x9]\n" // 读取连续的4个double数据

      "add x9, x9, x12\n"
      "ld1 {v14.2d}, [x9]\n" // 读取连续的2个double数据
      "add x9, x9, x10\n"
      "ld1 {v15.2d}, [x9]\n" // 读取连续的2个double数据

      "1:"
      "ld1 {v20.2d-v23.2d}, [x11]\n" // 读取连续的8个double数据

      // 计算
      // 第一行
      "dup v24.2d, v20.d[0]\n" // v24为临时变量
      "mov d16, v20.d[1]\n"
      "fmls d16, d24, v0.d[1]\n" // 乘减
      "mov v20.d[1], v16.d[0]\n"
      "dup v25.2d, v20.d[1]\n"       // 第二行临时变量
      "fmls v21.2d, v1.2d, v24.2d\n" // 乘减
      "fmls v22.2d, v2.2d, v24.2d\n" // 乘减
      "fmls v23.2d, v3.2d, v24.2d\n" // 乘减

      // 第二行
      "fmls v21.2d, v4.2d, v25.2d\n" // 乘减
      "fmls v22.2d, v5.2d, v25.2d\n" // 乘减
      "fmls v23.2d, v6.2d, v25.2d\n" // 乘减

      // 第三行
      "mov d16, v21.d[1]\n"
      "dup v24.2d, v21.d[0]\n"   // v24为临时变量
      "fmls d16, d21, v7.d[1]\n" // 乘减
      "mov v21.d[1], v16.d[0]\n"
      "fmls v22.2d, v8.2d, v24.2d\n" // 乘减
      "fmls v23.2d, v9.2d, v24.2d\n" // 乘减

      // 第四行
      "dup v25.2d, v21.d[1]\n"        // v24为临时变量
      "fmls v22.2d, v10.2d, v25.2d\n" // 乘减
      "fmls v23.2d, v11.2d, v25.2d\n" // 乘减

      // 第五行
      "mov d16, v22.d[1]\n"
      "dup v24.2d, v22.d[0]\n"    // v24为临时变量
      "fmls d16, d22, v12.d[1]\n" // 乘减
      "mov v22.d[1], v16.d[0]\n"
      "fmls v23.2d, v13.2d, v24.2d\n" // 乘减

      // 第六行
      "dup v24.2d, v22.d[1]\n"        // v24为临时变量
      "fmls v23.2d, v14.2d, v24.2d\n" // 乘减

      // 第七行
      "mov d16, v23.d[1]\n"
      "fmls d16, d23, v15.d[1]\n" // 乘减
      "mov v23.d[1], v16.d[0]\n"

      "st1 {v20.2d-v23.2d}, [x11]\n" // 读取连续的8个double数据

      "add x11, x11, x10\n" // BL指针向下偏移1行
      "sub %[m], %[m], #1\n"
      "cbnz %[m], 1b\n" // 跳转按列遍历

      :
      : [tri] "r"(tri), [cnt] "r"(n), [m] "r"(m_d)
      : "memory");
}

int __TEST5_org_inverse_matrix_n_Rectangle_8(double *src, double *dst, int n)
{
  /*
        LU分解的递归分块算法
  */

  int i, j, k, id;
  // int d = 2;
  // int loop_count = 3;
  // int d = 8;
  // int loop_count = 2;
  int d = 8;
  int loop_count = n/d -1;
  // int d = 16;
  // int loop_count = 31;
  // int d = 32;
  // int loop_count = 15;

  for (int i = 0; i < n * n; i++)
  {
    dst[i] = src[i];
  }

  for (int i = 0; i < loop_count; i++)
  {
    int temp1 = (n + 1) * d * i;
    int temp2 = n - d * i;
    // printf("temp2:%d\n",temp2);
    LU_Decomposition_L_8(dst + temp1, n, temp2, d);
    LU_Decomposition_TR_Square_8(dst + temp1, n, temp2, d);
    LU_L21_mul_U12_Rectangle_8(dst + temp1, n, temp2, d);
  }

  LU_Decomposition_BR_8(dst + loop_count * d * (n + 1), n, n - loop_count * d, d);
}
