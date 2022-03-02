#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define BLOCK_SIZE 64
#define BLOCK_SIZE_ROW 4
#define BLOCK_SIZE_COL 64
/*
        使用分块LU分解进行矩阵求逆（8*8正方块改进为 dn * dm（如8*16,8*32））
*/

/*
    测试不同大分块大小（d）对性能的影响
    #define BLOCK_SIZE 8
    #define BLOCK_SIZE_ROW 4
    #define BLOCK_SIZE_COL 16

*/

int __TEST5_org_inverse_matrix_n_Square(double *src, double *dst, int n0);
void print_mat(double *tri, int n, char *str);

/*********************TMP***************************/

int LU_Decomposition_8_mul_8(double *tri, int n)
{

  // 8*8
  /*
    v4   v5   v6   v7
    v8   v9   v10  v11
    v12  v13  v14  v15
    v16  v17  v18  v19
    v20  v21  v22  v23
    v24  v25  v26  v27
    v28  v29  v30  v31
    v12  v16  v20  v24
    */
}

/*********************PRINT***************************/

void print_mat(double *tri, int n, char *str)
{

  printf("tri_%s\n", str);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      printf("%f\t", tri[i * n + j]);
    }
    printf("\n");
  }
  printf("\n");
}

void print_L_mul_R(double *tri, int n)
{

  printf("L * R\n");
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {

      if (j <= i)
      {
        double res = tri[i * n + j];
        for (int k = 0; k < j; k++)
        {
          res += tri[i * n + k] * tri[k * n + j];
        }
        printf("%f\t", res);
      }
      else
      {
        double res = 0;
        for (int k = 0; k <= i; k++)
        {
          res += tri[i * n + k] * tri[k * n + j];
        }
        printf("%f\t", res);
      }
    }
    printf("\n");
  }
  printf("\n");
}

void check(double *src, double *tri, int n)
{

  printf("check...\n");
  double error = 0.0;
  double res = 0;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {

      if (j <= i)
      {
        res = tri[i * n + j];
        for (int k = 0; k < j; k++)
        {
          res += tri[i * n + k] * tri[k * n + j];
        }
      }
      else
      {
        res = 0.0;
        for (int k = 0; k <= i; k++)
        {
          res += tri[i * n + k] * tri[k * n + j];
        }
      }
      error += fabs(res - src[i * n + j]);
    }
  }
  printf("Accumulated error: %f\n", error);
}

/*********************ORIGIN***************************/

inline int LU_L21_mul_U12(double *dst, double *tri, int n, int m,
                                 int d)
{
  /*
        对L21,U12进行矩阵乘法，并与tri相减并保存在tri中
  */
  // 原代码
  int i, j, k;
  double *tri_in = tri + d * n;
  for (i = d; i < m; i++, tri_in += n)
  {
    for (j = d; j < m; j++)
    {
      double *tri_kn = tri;
      for (k = 0; k < d; k++, tri_kn += n)
      {
        *(tri_in + j) -= *(tri_in + k) * *(tri_kn + j);
      }
    }
  }
}


/*********************SQUARE BLOCKING***************************/

inline int LU_Decomposition_BR(double *tri, int n, int m, int d)
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

inline int LU_L21_mul_U12_Square(double *dst, double *tri, int n, int m,
                                 int d)
{
  /*
        对L21,U12进行矩阵乘法，并与tri相减并保存在tri中
  */
  // 原代码
  // int i, j, k;
  // double *tri_in = tri + d * n;
  // for (i = d; i < m; i++, tri_in += n) {
  //   for (j = d; j < m; j++) {
  //     double *tri_kn = tri;
  //     for (k = 0; k < d; k++, tri_kn += n) {
  //       *(tri_in + j) -= *(tri_in + k) * *(tri_kn + j);
  //     }
  //   }
  // }

  // 分块（错误）相比于原来的有提升，但提升不多
  // int i, j, jj, k;
  // double *tri_in = tri + d * n;
  // for (i = d; i < m; i++, tri_in += n) {
  //   for (j = d; j < m; j += BLOCK_SIZE) {
  //     double *tri_kn = tri;
  //     double *tri_in_k = tri_in;
  //     int min_j = (BLOCK_SIZE + j < m) ? (BLOCK_SIZE + j) : (m);
  //     for (k = 0; k < d; k++, tri_in_k++, tri_kn += n) {
  //       for (jj = j; jj < min_j; jj++) {
  //         *(tri_in + jj) -= *(tri_in_k) * *(tri_kn + jj);
  //       }
  //     }
  //   }
  // }

  // 按列分块（正确）（方块）
  int i, j, ii, jj, k;
  double *tri_in = tri + d * n;
  for (i = d; i < m; i += BLOCK_SIZE, tri_in += n * BLOCK_SIZE)
  {
    int min_i = (BLOCK_SIZE + i < m) ? (BLOCK_SIZE + i) : (m);
    for (j = d; j < m; j += BLOCK_SIZE)
    {
      int min_j = (BLOCK_SIZE + j < m) ? (BLOCK_SIZE + j) : (m);
      double *tri_iin = tri_in;
      for (ii = i; ii < min_i; ii++, tri_iin += n)
      {
        double *tri_kn = tri;
        double *tri_in_k = tri_iin;
        for (k = 0; k < d; k++, tri_in_k++, tri_kn += n)
        {
          for (jj = j; jj < min_j; jj++)
          {
            *(tri_iin + jj) -= *(tri_in_k) * *(tri_kn + jj);
          }
        }
      }
    }
  }
}

inline int LU_Decomposition_TR_Square(double *dst, double *tri, int n, int m,
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

inline int LU_Decomposition_L(double *dst, double *tri, int n, int m, int d)
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

int __TEST5_org_inverse_matrix_n_Square(double *src, double *dst, int n)
{
  /*
        LU分解的递归分块算法
  */

  int i, j, k, id;
  // int d = 2;
  // int loop_count = 3;
  // int d = 4;
  // int loop_count = 127;
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
    LU_Decomposition_L(dst + temp1, dst + temp1, n, temp2, d);
    LU_Decomposition_TR_Square(dst + temp1, dst + temp1, n, temp2, d);
    LU_L21_mul_U12_Square(dst + temp1, dst + temp1, n, temp2, d);
  }

  LU_Decomposition_BR(dst + loop_count * d * (n + 1), n, n - loop_count * d, d);
}

int __TEST5_org_inverse_matrix_n(double *src, double *dst, int n)
{
  /*
        LU分解的递归分块算法
  */

  int i, j, k, id;
  // int d = 2;
  // int loop_count = 3;
  // int d = 4;
  // int loop_count = 127;
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
    LU_Decomposition_L(dst + temp1, dst + temp1, n, temp2, d);
    LU_Decomposition_TR_Square(dst + temp1, dst + temp1, n, temp2, d);
    LU_L21_mul_U12(dst + temp1, dst + temp1, n, temp2, d);
  }

  LU_Decomposition_BR(dst + loop_count * d * (n + 1), n, n - loop_count * d, d);
}


/*********************RECTANGLE BLOCKING***************************/

inline int LU_L21_mul_U12_Rectangle1(double *dst, double *tri, int n, int m,
                                     int d)
{
  /*
        对L21,U12进行矩阵乘法，并与tri相减并保存在tri中
  */

  // 按列分块（正确）（长方块）
  int i, j, ii, jj, k;
  double *tri_in = tri + d * n;
  for (i = d; i < m; i += BLOCK_SIZE_ROW, tri_in += n * BLOCK_SIZE_ROW)
  {
    int min_i = (BLOCK_SIZE_ROW + i < m) ? (BLOCK_SIZE_ROW + i) : (m);
    for (j = d; j < m; j += BLOCK_SIZE_COL)
    {
      int min_j = (BLOCK_SIZE_COL + j < m) ? (BLOCK_SIZE_COL + j) : (m);
      double *tri_iin = tri_in;
      for (ii = i; ii < min_i; ii++, tri_iin += n)
      {
        double *tri_kn = tri;
        double *tri_in_k = tri_iin;
        for (k = 0; k < d; k++, tri_in_k++, tri_kn += n)
        {
          for (jj = j; jj < min_j; jj++)
          {
            *(tri_iin + jj) -= *(tri_in_k) * *(tri_kn + jj);
          }
        }
      }
    }
  }
}

inline int LU_L21_mul_U12_Rectangle2(double *dst, double *tri, int n, int m,
                                     int d)
{
  /*
        对L21,U12进行矩阵乘法，并与tri相减并保存在tri中
        对分块矩阵内部进一步进行分块（与汇编分块具有相同的计算逻辑）
  */

  int i, j, ii, jj, k, jjj;

  double *tri_in = tri + d * n;
  for (i = d; i < m; i += BLOCK_SIZE_ROW, tri_in += BLOCK_SIZE_ROW * n)
  {
    int min_i = (BLOCK_SIZE_ROW + i < m) ? (BLOCK_SIZE_ROW + i) : (m);
    for (j = d; j < m; j += BLOCK_SIZE_COL)
    {
      int min_j = (BLOCK_SIZE_COL + j < m) ? (BLOCK_SIZE_COL + j) : (m);
      double *tri_iin;
      for (jj = j; jj < min_j; jj += 8)
      {
        tri_iin = tri_in;
        for (ii = i; ii < min_i; ii++, tri_iin += n)
        {
          double *tri_iin_k = tri_iin;
          double *tri_kn_jj = tri + jj;
          for (k = 0; k < d / 2; k++, tri_iin_k++, tri_kn_jj += n)
          {
            *(tri_iin + jj) -= *tri_iin_k * *(tri_kn_jj);
            *(tri_iin + jj + 1) -= *tri_iin_k * *(tri_kn_jj + 1);
            *(tri_iin + jj + 2) -= *tri_iin_k * *(tri_kn_jj + 2);
            *(tri_iin + jj + 3) -= *tri_iin_k * *(tri_kn_jj + 3);
            *(tri_iin + jj + 4) -= *tri_iin_k * *(tri_kn_jj + 4);
            *(tri_iin + jj + 5) -= *tri_iin_k * *(tri_kn_jj + 5);
            *(tri_iin + jj + 6) -= *tri_iin_k * *(tri_kn_jj + 6);
            *(tri_iin + jj + 7) -= *tri_iin_k * *(tri_kn_jj + 7);
          }
        }
      }
      for (jj = j; jj < min_j; jj += 8)
      {
        tri_iin = tri_in;
        for (ii = i; ii < min_i; ii++, tri_iin += n)
        {
          double *tri_iin_k = tri_iin + d / 2;
          double *tri_kn_jj = tri + d / 2 * n + jj;
          for (k = d / 2; k < d; k++, tri_iin_k++, tri_kn_jj += n)
          {
            *(tri_iin + jj) -= *tri_iin_k * *(tri_kn_jj);
            *(tri_iin + jj + 1) -= *tri_iin_k * *(tri_kn_jj + 1);
            *(tri_iin + jj + 2) -= *tri_iin_k * *(tri_kn_jj + 2);
            *(tri_iin + jj + 3) -= *tri_iin_k * *(tri_kn_jj + 3);
            *(tri_iin + jj + 4) -= *tri_iin_k * *(tri_kn_jj + 4);
            *(tri_iin + jj + 5) -= *tri_iin_k * *(tri_kn_jj + 5);
            *(tri_iin + jj + 6) -= *tri_iin_k * *(tri_kn_jj + 6);
            *(tri_iin + jj + 7) -= *tri_iin_k * *(tri_kn_jj + 7);
          }
        }
      }
    }
  }
}

inline int LU_Decomposition_TR_Rectangle(double *dst, double *tri, int n, int m,
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

int __TEST5_org_inverse_matrix_n_Rectangle1(double *src, double *dst, int n)
{
  /*
        LU分解的递归分块算法
  */

  int i, j, k, id;
  // int d = 2;
  // int loop_count = 2;
  // int d = 4;
  // int loop_count = 127;
  int d = 8;
  int loop_count = n / d - 1;
  // int d = 16;
  // int loop_count = 31;
  // int d = 32;
  // int loop_count = 15;

  // 上三角矩阵初始化
  // dst初始化为单位矩阵
  for (int i = 0; i < n * n; i++)
  {
    dst[i] = src[i];
  }

  for (int i = 0; i < loop_count; i++)
  {
    int temp1 = (n + 1) * d * i;
    int temp2 = n - d * i;
    LU_Decomposition_L(dst + temp1, dst + temp1, n, temp2, d);
    LU_Decomposition_TR_Rectangle(dst + temp1, dst + temp1, n, temp2, d);
    LU_L21_mul_U12_Rectangle1(dst + temp1, dst + temp1, n, temp2, d);
  }

  LU_Decomposition_BR(dst + loop_count * d * (n + 1), n, n - loop_count * d, d);
}

int __TEST5_org_inverse_matrix_n_Rectangle2(double *src, double *dst, int n)
{
  /*
        LU分解的递归分块算法
  */

  int i, j, k, id;
  // int d = 2;
  // int loop_count = 2;
  // int d = 4;
  // int loop_count = 127;
  int d = 8;
  int loop_count = n / d - 1;
  // int d = 16;
  // int loop_count = 31;
  // int d = 32;
  // int loop_count = 15;

  // 上三角矩阵初始化
  // dst初始化为单位矩阵
  for (int i = 0; i < n * n; i++)
  {
    dst[i] = src[i];
  }

  for (int i = 0; i < loop_count; i++)
  {
    int temp1 = (n + 1) * d * i;
    int temp2 = n - d * i;
    LU_Decomposition_L(dst + temp1, dst + temp1, n, temp2, d);
    LU_Decomposition_TR_Rectangle(dst + temp1, dst + temp1, n, temp2, d);
    LU_L21_mul_U12_Rectangle2(dst + temp1, dst + temp1, n, temp2, d);
  }

  LU_Decomposition_BR(dst + loop_count * d * (n + 1), n, n - loop_count * d, d);
}


/*********************ORIGIN---JKI***************************/

inline int LU_L21_mul_U12_opt1(double *dst, double *tri, int n, int m,
                                 int d)
{
  /*
        对L21,U12进行矩阵乘法，并与tri相减并保存在tri中
  */
  // 原代码
  int i, j, k;
  // double *tri_in = tri + d * n;
  // for (i = d; i < m; i++, tri_in += n)
  // {
  //   for (j = d; j < m; j++)
  //   {
  //     double *tri_kn = tri;
  //     for (k = 0; k < d; k++, tri_kn += n)
  //     {
  //       *(tri_in + j) -= *(tri_in + k) * *(tri_kn + j);
  //     }
  //   }
  // }

  // jki:3,432,818 us
  // for(j=d; j<m; j++){
  //   for(k=0;k<d;k++){
  //     for(i=d;i<m;i++){
  //       tri[i *n+j] -= tri[i*n+k] * tri[k*n+j];
  //     }
  //   }
  // }

  // ijk:412,541 us
  // for(i=d;i<m;i++){
  //   for(j=d; j<m; j++){
  //     for(k=0;k<d;k++){
  //       tri[i *n+j] -= tri[i*n+k] * tri[k*n+j];
  //     }
  //   }
  // }

  // ikj:416,698 us
  // for(i=d;i<m;i++){
  //   for(k=0;k<d;k++){
  //     for(j=d; j<m; j++){
  //       tri[i *n+j] -= tri[i*n+k] * tri[k*n+j];
  //     }
  //   }
  // }
  
  // ikj:425,750 us
  double *tri_in = tri + d * n;
  for(i=d;i<m;i++, tri_in+=n){
    double *tri_kn = tri;
    for(k=0;k<d;k++, tri_kn += n){
      for(j=d; j<m; j++){
        *(tri_in + j) -= *(tri_in + k) * *(tri_kn+j);
      }
    }
  }

  // kij:435,198 us
  // for(k=0;k<d;k++){
  //   for(i=d;i<m;i++){
  //     for(j=d; j<m; j++){
  //       tri[i *n+j] -= tri[i*n+k] * tri[k*n+j];
  //     }
  //   }
  // }

}

int __TEST5_org_inverse_matrix_opt1_n(double *src, double *dst, int n)
{
  /*
        LU分解的递归分块算法
  */

  int i, j, k, id;
  // int d = 2;
  // int loop_count = 3;
  // int d = 4;
  // int loop_count = 127;
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
    LU_Decomposition_L(dst + temp1, dst + temp1, n, temp2, d);
    LU_Decomposition_TR_Square(dst + temp1, dst + temp1, n, temp2, d);
    LU_L21_mul_U12_opt1(dst + temp1, dst + temp1, n, temp2, d);
  }

  LU_Decomposition_BR(dst + loop_count * d * (n + 1), n, n - loop_count * d, d);
}


// int main() {

//   // double src[36] = {2, 3, 1,  5, 1, 6, 3, 7, 9, 1, 3, 4, 6, 1, 7,  8,  6,
//   2,
//   //                   4, 6, 10, 2, 5, 7, 9, 1, 2, 3, 4, 5, 8, 9, 10, 23, 12,
//   //                   7};
//   // double dst[36] = {0};
//   // double tri[36] = {0};

//   int d = 1;
//   int n = 5;

//   double *src = malloc(sizeof(double) * n * n);
//   double *dst = malloc(sizeof(double) * n * n);
//   double *tri = malloc(sizeof(double) * n * n);

//   srand((unsigned)time(NULL));
//   for (int i = 0; i < n; i++) {
//     for (int j = 0; j < n; j++) {
//       src[i * n + j] = rand() % 10000;
//       // printf("%f", src[i * n + j]);
//     }
//   }

//   __TEST5_org_inverse_matrix_n(src, dst, n);

//   return 0;
// }
