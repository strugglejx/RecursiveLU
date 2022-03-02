#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*
        使用LU分解进行矩阵求逆
*/

int __TEST3_opt6_inverse_matrix_n(double *src, double *dst, int n0);
int __TEST3_opt5_f_inverse_matrix_n(double *src, double *dst, int n0);
int __TEST3_opt40_inverse_matrix_n(double *src, double *dst, int n0);
int __TEST3_opt41_inverse_matrix_n(double *src, double *dst, int n0);
int __TEST3_opt42_inverse_matrix_n(double *src, double *dst, int n0);
int __TEST3_opt31_inverse_matrix_n(double *src, double *dst, int n0);
int __TEST3_opt3_inverse_matrix_n(double *src, double *dst, int n0);
int __TEST3_opt22_inverse_matrix_n(double *src, double *dst, int n0);
int __TEST3_opt21_inverse_matrix_n(double *src, double *dst, int n0);
int __TEST3_opt20_inverse_matrix_n(double *src, double *dst, int n0);
int __TEST3_opt1_inverse_matrix_n(double *src, double *dst, int n0);
int __TEST3_org_inverse_matrix_n(double *src, double *dst, int n0);

int __TEST3_opt6_inverse_matrix_n(double *src, double *dst, int n) {
  /*
          LU/Crout分解（其中，L为下三角矩阵，U为单位上三角矩阵）
          优化1：LU保存在同一个矩阵中
          优化2：所有二维数组改为一维数组；
                顺序遍历时设置一个变量id进行索引，以代替乘法索引得到性能的提升
                in和knj提前计算
          优化3：按列遍历
          优化4：数组改为指针
  */

  return 0;
}

int __TEST3_opt5_f_inverse_matrix_n(double *src, double *dst, int n) {
  /*
          LU/Crout分解（其中，L为下三角矩阵，U为单位上三角矩阵）
          优化1：LU保存在同一个矩阵中
          优化2：所有二维数组改为一维数组；
                顺序遍历时设置一个变量id进行索引，以代替乘法索引得到性能的提升
                in和knj提前计算
          优化3：按列遍历
          优化4：数组改为指针
          优化5：动态数组改静态数组
  */

  int i, j, k, id;
  // 用于保存单位上三角和下三角矩阵
  // double *tri = (double *)malloc(n * n * sizeof(double));
  // memset(tri, 0, sizeof(tri));
  memset(dst, 0, sizeof(dst));

  // 上三角矩阵初始化
  // dst初始化为单位矩阵
  double *dst_ii = dst;
  int np1 = n + 1;
  for (i = 0; i < n; i++,dst_ii += np1) {
    *dst_ii = 1;
  }

  // Crout分解——计算出L和U
  // 优化2：
  double *dst_id = dst;
  double *src_id = src;
  for (id = 0, j = 0; j < n; j++, id = j, dst_id = dst + j, src_id = src + j) {
    double *dst_in = dst;
    for (i = 0; i < n; i++, dst_in += n, id += n, dst_id += n, src_id += n) {
      if (j <= i) { // 下三角
        *dst_id = *src_id;
        double *dst_kn_j = dst + j;
        double *dst_in_j = dst_in + j;
        double *dst_in_k = dst_in;
        for (; dst_in_k < dst_in_j; dst_in_k++, dst_kn_j += n) {
          *dst_id -= *dst_in_k * *dst_kn_j;
        }
      } else { // 上三角
        *dst_id = *src_id;
        double *dst_kn_j = dst + j;
        double *dst_in_i = dst_in + i;
        double *dst_in_k = dst_in;
        for (; dst_in_k < dst_in_i; dst_in_k++, dst_kn_j += n) {
          *dst_id -= *dst_in_k * *dst_kn_j;
        }
        *dst_id /= *dst_in_i;
      }
    }
  }
}

int __TEST3_opt42_inverse_matrix_n(double *src, double *dst, int n) {
  /*
          LU/Crout分解（其中，L为下三角矩阵，U为单位上三角矩阵）
          优化1：LU保存在同一个矩阵中
          优化2：所有二维数组改为一维数组；
                顺序遍历时设置一个变量id进行索引，以代替乘法索引得到性能的提升
                in和knj提前计算
          优化3：按列遍历
          优化4：数组改为指针
          优化42：在31的基础上将数组改为指针
  */
  int i, j, k, id;
  // 用于保存单位上三角和下三角矩阵
  memset(dst, 0, sizeof(double) * n * n);

  // 上三角矩阵初始化
  // dst初始化为单位矩阵
  double *dst_ii = dst;
  int np1 = n + 1;
  for (i = 0; i < n; i++, dst_ii += np1) {
    *dst_ii = 1;
  }

  // Crout分解——计算出L和U
  // 优化2：
  double *dst_id = dst;
  double *src_id = src;
  for (j = 0; j < n; j++, id = j, dst_id = dst + j, src_id = src + j) {
    double *dst_in = dst;
    for (i = 0; i < n; i++, dst_in += n, id += n, dst_id += n, src_id += n) {
      if (j <= i) { // 下三角
        *dst_id = *src_id;
        double *dst_kn_j = dst + j;
        double *dst_in_j = dst_in + j;
        double *dst_in_k = dst_in;
        for (; dst_in_k < dst_in_j; dst_in_k++, dst_kn_j += n) {
          *dst_id -= *dst_in_k * *dst_kn_j;
        }
      } else { // 上三角
        *dst_id = *src_id;
        double *dst_kn_j = dst + j;
        double *dst_in_i = dst_in + i;
        double *dst_in_k = dst_in;
        for (; dst_in_k < dst_in_i; dst_in_k++, dst_kn_j += n) {
          *dst_id -= *dst_in_k * *dst_kn_j;
        }
        *dst_id /= *dst_in_i;
      }
    }
  }
}

int __TEST3_opt41_inverse_matrix_n(double *src, double *dst, int n) {
  /*
          LU/Crout分解（其中，L为下三角矩阵，U为单位上三角矩阵）
          优化1：LU保存在同一个矩阵中
          优化2：所有二维数组改为一维数组；
                顺序遍历时设置一个变量id进行索引，以代替乘法索引得到性能的提升
                in和knj提前计算
          优化3：按列遍历
          优化4：数组改为指针
  */

  int i, j, k, id;
  // 用于保存单位上三角和下三角矩阵
  double *tri = (double *)malloc(n * n * sizeof(double));

  memset(tri, 0, sizeof(tri));
  memset(dst, 0, sizeof(dst));

  // 上三角矩阵初始化
  // dst初始化为单位矩阵
  double *tri_ii = tri;
  double *dst_ii = dst;
  int np1 = n + 1;
  for (i = 0; i < n; i++, tri_ii += np1, dst_ii += np1) {
    *tri_ii = 1; // 在循环中*(tri + i * n + i)与tri[i * n + i]等价
    *dst_ii = 1;
  }

  // Crout分解——计算出L和U
  // 优化2：
  double *tri_id = tri;
  double *src_id = src;
  for (j = 0; j < n; j++, id = j, tri_id = tri + j, src_id = src + j) {
    double *tri_in = tri;
    for (i = 0; i < n; i++, tri_in += n, id += n, tri_id += n, src_id += n) {
      if (j <= i) { // 下三角
        *tri_id = *src_id;
        double *tri_kn_j = tri + j;
        double *tri_in_j = tri_in + j;
        double *tri_in_k = tri_in;
        for (; tri_in_k < tri_in_j; tri_in_k++, tri_kn_j += n) {
          *tri_id -= *tri_in_k * *tri_kn_j;
        }
      } else { // 上三角
        *tri_id = *src_id;
        double *tri_kn_j = tri + j;
        double *tri_in_i = tri_in + i;
        double *tri_in_k = tri_in;
        for (; tri_in_k < tri_in_i; tri_in_k++, tri_kn_j += n) {
          *tri_id -= *tri_in_k * *tri_kn_j;
        }
        *tri_id /= *tri_in_i;
      }
    }
  }

  // 求解Ly=b（顺序）
  double *tri_in_j = tri;
  double *dst_in_k = dst;
  double *dst_k = dst;
  for (k = 0; k < n; k++, dst_k++, dst_in_k = dst_k) {
    tri_in_j = tri;
    for (int in = 0, i = 0; i < n;
         in += n, tri_in_j += (n - i), i++, dst_in_k += n) {
      double *dst_jn_k = dst_k;
      for (j = 0; j < i; j++, tri_in_j++, dst_jn_k += n) {
        *dst_in_k -= *tri_in_j * *dst_jn_k;
      }
      *dst_in_k /= *tri_in_j; // *(dst + in + k) /= *(tri + in + i);
    }
  }

  // 求解Ux=y（逆序）
  dst_in_k = dst;
  dst_k = dst;
  int temp_n_n1 = n * (n - 1);
  for (k = 0; k < n; k++, dst_k++) {
    dst_in_k = dst_k + temp_n_n1;
    double *tri_in_j = tri + temp_n_n1 + n - 1;
    for (
        i = n - 1; i >= 0;
        i--, dst_in_k -= n,
       tri_in_j -=
        (i +
         2)) { // 内循环结束后j等于i，由于前面i--了，且in-=n，另外内循环j初始化为n-1。因此tri_in_j
               // -= (i+2)
      double *dst_jn_k = dst_k + temp_n_n1;
      for (j = n - 1; j > i; j--, dst_jn_k -= n, tri_in_j--) {
        *dst_in_k -=
            *(tri_in_j) * *dst_jn_k; // tri_in_j通常调试并找到相应规律后得到的
      }
    }
  }

  free(tri);

  return 0;
}

int __TEST3_opt40_inverse_matrix_n(double *src, double *dst, int n) {
  /*
          LU/Crout分解（其中，L为下三角矩阵，U为单位上三角矩阵）
          优化1：LU保存在同一个矩阵中
          优化2：所有二维数组改为一维数组；
                顺序遍历时设置一个变量id进行索引，以代替乘法索引得到性能的提升
                in和knj提前计算
          优化3：按列遍历
          优化4：数组改为指针
  */

  int i, j, k, id;
  // 用于保存单位上三角和下三角矩阵
  double *tri = (double *)malloc(n * n * sizeof(double));

  memset(tri, 0, sizeof(tri));
  memset(dst, 0, sizeof(dst));

  // 上三角矩阵初始化
  // dst初始化为单位矩阵
  for (i = 0; i < n; i++) {
    *(tri + i * n + i) = 1; // 在循环中*(tri + i * n + i)与tri[i * n + i]等价
    *(dst + i * n + i) = 1;
  }

  // Crout分解——计算出L和U
  // 优化2：
  for (id = 0, j = 0; j < n; j++, id = j) {
    for (int in = 0, i = 0; i < n; i++, in += n, id += n) {
      if (j <= i) { // 下三角
        *(tri + id) = *(src + id);
        for (int kn_j = j, k = 0; k < j; k++, kn_j += n) {
          *(tri + id) -= *(tri + in + k) * *(tri + kn_j);
        }
      } else { // 上三角
        *(tri + id) = *(src + id);
        for (int kn_j = j, k = 0; k < i; k++, kn_j += n) {
          *(tri + id) -= *(tri + in + k) * *(tri + kn_j);
        }
        *(tri + id) /= *(tri + in + i);
      }
    }
  }

  // 求解Ly=b（顺序）
  for (k = 0; k < n; k++) {
    for (int in = 0, i = 0; i < n; i++, in += n) {
      for (j = 0; j < i; j++) {
        *(dst + in + k) -= *(tri + in + j) * *(dst + j * n + k);
      }
      *(dst + in + k) /= *(tri + in + i);
    }
  }

  // 求解Ux=y（逆序）
  for (k = 0; k < n; k++) {
    for (int in = (n - 1) * n, i = n - 1; i >= 0; i--, in -= n) {
      for (j = n - 1; j > i; j--) {
        *(dst + in + k) -= *(tri + in + j) * *(dst + j * n + k);
      }
    }
  }

  free(tri);

  return 0;
}

int __TEST3_opt31_inverse_matrix_n(double *src, double *dst, int n) {
  /*
          LU/Crout分解（其中，L为下三角矩阵，U为单位上三角矩阵）
          优化1：LU保存在同一个矩阵中
          优化2：所有二维数组改为一维数组；
                顺序遍历时设置一个变量id进行索引，以代替乘法索引得到性能的提升
                in和knj提前计算
          优化3：按列遍历
          优化31：回代部分按行遍历
  */

  int i, j, k, id;
  // 用于保存单位上三角和下三角矩阵
  double *tri = (double *)malloc(n * n * sizeof(double));

  memset(tri, 0, sizeof(tri));
  memset(dst, 0, sizeof(dst));

  // 上三角矩阵初始化
  // dst初始化为单位矩阵
  for (i = 0; i < n; i++) {
    tri[i * n + i] = 1;
    dst[i * n + i] = 1;
  }

  // Crout分解——计算出L和U
  // 优化2：
  for (id = 0, j = 0; j < n; j++, id = j) {
    for (int in = 0, i = 0; i < n; i++, in += n, id += n) {
      if (j <= i) { // 下三角
        tri[id] = src[id];
        for (int kn_j = j, k = 0; k < j; k++, kn_j += n) {
          tri[id] -= tri[in + k] * tri[kn_j];
        }
      } else { // 上三角
        tri[id] = src[id];
        for (int kn_j = j, k = 0; k < i; k++, kn_j += n) {
          tri[id] -= tri[in + k] * tri[kn_j];
        }
        tri[id] /= tri[in + i];
      }
    }
  }

  // 求解Ly=b（顺序）
  for (int in = 0, i = 0; i < n; i++, in += n) {
    for (j = 0; j < i; j++) {
      for (k = 0; k < n; k++) {
        dst[in + k] -= tri[in + j] * dst[j * n + k];
      }
    }
    for (k = 0; k < n; k++) {
      dst[in + k] /= tri[in + i];
    }
  }

  // 求解Ux=y（逆序）
  for (int in = (n - 1) * n, i = n - 1; i >= 0; i--, in -= n) {
    for (j = n - 1; j > i; j--) {
      for (k = 0; k < n; k++) {
        dst[in + k] -= tri[in + j] * dst[j * n + k];
      }
    }
  }

  free(tri);
}

int __TEST3_opt3_inverse_matrix_n(double *src, double *dst, int n) {
  /*
          LU/Crout分解（其中，L为下三角矩阵，U为单位上三角矩阵）
          优化1：LU保存在同一个矩阵中
          优化2：所有二维数组改为一维数组；
                顺序遍历时设置一个变量id进行索引，以代替乘法索引得到性能的提升
                in和knj提前计算
          优化3：按列遍历
  */

  int i, j, k, id;
  // 用于保存单位上三角和下三角矩阵
  double *tri = (double *)malloc(n * n * sizeof(double));

  memset(tri, 0, sizeof(tri));
  memset(dst, 0, sizeof(dst));

  // 上三角矩阵初始化
  // dst初始化为单位矩阵
  for (i = 0; i < n; i++) {
    tri[i * n + i] = 1;
    dst[i * n + i] = 1;
  }

  // Crout分解——计算出L和U
  // 优化2：
  for (id = 0, j = 0; j < n; j++, id = j) {
    for (int in = 0, i = 0; i < n; i++, in += n, id += n) {
      if (j <= i) { // 下三角
        tri[id] = src[id];
        for (int kn_j = j, k = 0; k < j; k++, kn_j += n) {
          tri[id] -= tri[in + k] * tri[kn_j];
        }
      } else { // 上三角
        tri[id] = src[id];
        for (int kn_j = j, k = 0; k < i; k++, kn_j += n) {
          tri[id] -= tri[in + k] * tri[kn_j];
        }
        tri[id] /= tri[in + i];
      }
    }
  }

  // 求解Ly=b（顺序）
  for (k = 0; k < n; k++) {
    for (int in = 0, i = 0; i < n; i++, in += n) {
      for (j = 0; j < i; j++) {
        dst[in + k] -= tri[in + j] * dst[j * n + k];
      }
      dst[in + k] /= tri[in + i];
    }
  }

  // 求解Ux=y（逆序）
  for (k = 0; k < n; k++) {
    for (int in = (n - 1) * n, i = n - 1; i >= 0; i--, in -= n) {
      for (j = n - 1; j > i; j--) {
        dst[in + k] -= tri[in + j] * dst[j * n + k];
      }
    }
  }

  free(tri);
}

int __TEST3_opt22_inverse_matrix_n(double *src, double *dst, int n) {
  /*
          LU/Crout分解（其中，L为下三角矩阵，U为单位上三角矩阵）
          优化1：LU保存在同一个矩阵中
          优化2：所有二维数组改为一维数组；
                顺序遍历时设置一个变量id进行索引，以代替乘法索引得到性能的提升
                in和knj提前计算
  */

  int i, j, k, id;
  // 用于保存单位上三角和下三角矩阵
  double *tri = (double *)malloc(n * n * sizeof(double));

  memset(tri, 0, sizeof(tri));
  memset(dst, 0, sizeof(dst));

  // 上三角矩阵初始化
  // dst初始化为单位矩阵
  for (int in = 0, i = 0; i < n; i++, in += n) {
    tri[in + i] = 1;
    dst[in + i] = 1;
  }

  // Crout分解——计算出L和U
  // 优化2：
  for (int in = 0, i = 0, id = 0; i < n; i++, in += n) {
    for (j = 0; j < n; j++, id++) {

      if (j <= i) { // 下三角
        tri[id] = src[id];
        for (int kn_j = j, k = 0; k < j; k++, kn_j += n) {
          tri[id] -= tri[in + k] * tri[kn_j];
        }
      } else { // 上三角
        tri[id] = src[id];
        for (int kn_j = j, k = 0; k < i; k++, kn_j += n) {
          tri[id] -= tri[in + k] * tri[kn_j];
        }
        tri[id] /= tri[in + i];
      }
    }
  }

  // 求解Ly=b（顺序）
  for (k = 0; k < n; k++) {
    for (int in = 0, i = 0; i < n; i++, in += n) {
      for (j = 0; j < i; j++) {
        dst[in + k] -= tri[in + j] * dst[j * n + k];
      }
      dst[in + k] /= tri[in + i];
    }
  }

  // 求解Ux=y（逆序）
  for (k = 0; k < n; k++) {
    for (int in = (n - 1) * n, i = n - 1; i >= 0; i--, in -= n) {
      for (j = n - 1; j > i; j--) {
        dst[in + k] -= tri[in + j] * dst[j * n + k];
      }
    }
  }

  free(tri);
}

int __TEST3_opt21_inverse_matrix_n(double *src, double *dst, int n) {
  /*
          LU/Crout分解（其中，L为下三角矩阵，U为单位上三角矩阵）
          优化1：LU保存在同一个矩阵中
          优化2：所有二维数组改为一维数组；
                顺序遍历时设置一个变量id进行索引，以代替乘法索引得到性能的提升
  */

  int i, j, k, id;
  // 用于保存单位上三角和下三角矩阵
  double *tri = (double *)malloc(n * n * sizeof(double));

  memset(tri, 0, sizeof(tri));
  memset(dst, 0, sizeof(dst));

  // 上三角矩阵初始化
  // dst初始化为单位矩阵
  for (i = 0; i < n; i++) {
    tri[i * n + i] = 1;
    dst[i * n + i] = 1;
  }

  // Crout分解——计算出L和U
  // 优化2：
  for (i = 0, id = 0; i < n; i++) {
    for (j = 0; j < n; j++, id++) {
      // int in = i * n;

      if (j <= i) { // 下三角
        tri[id] = src[id];
        for (k = 0; k < j; k++) {
          tri[id] -= tri[i * n + k] * tri[k * n + j];
        }
      } else { // 上三角
        tri[id] = src[id];
        for (k = 0; k < i; k++) {
          tri[id] -= tri[i * n + k] * tri[k * n + j];
        }
        tri[id] /= tri[i * n + i];
      }
    }
  }

  // 求解Ly=b（顺序）
  for (k = 0; k < n; k++) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < i; j++) {
        dst[i * n + k] -= tri[i * n + j] * dst[j * n + k];
      }
      dst[i * n + k] /= tri[i * n + i];
    }
  }

  // 求解Ux=y（逆序）
  for (k = 0; k < n; k++) {
    for (i = n - 1; i >= 0; i--) {
      for (j = n - 1; j > i; j--) {
        dst[i * n + k] -= tri[i * n + j] * dst[j * n + k];
      }
    }
  }

  free(tri);
}

int __TEST3_opt20_inverse_matrix_n(double *src, double *dst, int n) {
  /*
          LU/Crout分解（其中，L为下三角矩阵，U为单位上三角矩阵）
          优化1：LU保存在同一个矩阵中
          优化20：所有二维数组改为一维数组；
  */

  int i, j, k;

  // Crout分解——计算出L和U
  // 优化2：
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (j <= i) { // 下三角
        dst[i * n + j] = src[i * n + j];
        for (k = 0; k < j; k++) {
          dst[i * n + j] -= dst[i * n + k] * dst[k * n + j];
        }
      } else { // 上三角
        dst[i * n + j] = src[i * n + j];
        for (k = 0; k < i; k++) {
          dst[i * n + j] -= dst[i * n + k] * dst[k * n + j];
        }
        dst[i * n + j] /= dst[i * n + i];
      }
    }
  }

}

int __TEST3_opt1_inverse_matrix_n(double *src, double *dst, int n) {
  /*
          LU/Crout分解（其中，L为下三角矩阵，U为单位上三角矩阵）
          优化1：LU保存在同一个矩阵中
  */

  int i, j, k;
  // 用于保存单位上三角和下三角矩阵
  double **tri = (double **)malloc(n * sizeof(double *));
  for (i = 0; i < n; i++) {
    tri[i] = (double *)malloc(n * sizeof(double));
  }

  // L、U和dst初始化
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      tri[i][j] = 0;
      dst[i * n + j] = 0;
    }
  }

  // 上三角矩阵初始化
  // dst初始化为单位矩阵
  for (i = 0; i < n; i++) {
    tri[i][i] = 1;
    dst[i * n + i] = 1;
  }

  // Crout分解——计算出L和U
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (j <= i) {
        tri[i][j] = src[i * n + j];
        for (k = 0; k < j; k++) {
          tri[i][j] -= tri[i][k] * tri[k][j];
        }
      } else {
        tri[i][j] = src[i * n + j];
        for (k = 0; k < i; k++) {
          tri[i][j] -= tri[i][k] * tri[k][j];
        }
        tri[i][j] /= tri[i][i];
      }
    }
  }

  // 求解Ly=b（顺序）
  for (k = 0; k < n; k++) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < i; j++) {
        dst[i * n + k] -= tri[i][j] * dst[j * n + k];
      }
      dst[i * n + k] /= tri[i][i];
    }
  }

  // 求解Ux=y（逆序）
  for (k = 0; k < n; k++) {
    for (i = n - 1; i >= 0; i--) {
      for (j = n - 1; j > i; j--) {
        dst[i * n + k] -= tri[i][j] * dst[j * n + k];
      }
    }
  }

  for (i = 0; i < n; i++) {
    free(tri[i]);
  }
  free(tri);
}

int __TEST3_org_inverse_matrix_n(double *src, double *dst, int n) {
  /*
          LU分解（其中，L为单位下三角矩阵，U为上三角矩阵）
  */

  int i, j, k;
  // 下三角
  double **L = (double **)malloc(n * sizeof(double *));
  for (i = 0; i < n; i++) {
    L[i] = (double *)malloc(n * sizeof(double));
  }

  // 单位上三角
  double **U = (double **)malloc(n * sizeof(double *));
  for (i = 0; i < n; i++) {
    U[i] = (double *)malloc(n * sizeof(double));
  }

  // L、U和dst初始化
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      L[i][j] = 0;
      U[i][j] = 0;
      dst[i * n + j] = 0;
    }
  }

  // 上三角矩阵初始化
  // dst初始化为单位矩阵
  for (i = 0; i < n; i++) {
    U[i][i] = 1;
    dst[i * n + i] = 1;
  }

  //
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (j <= i) {
        L[i][j] = src[i * n + j];
        for (k = 0; k < j; k++) {
          L[i][j] -= L[i][k] * U[k][j];
        }
      } else {
        U[i][j] = src[i * n + j];
        for (k = 0; k < i; k++) {
          U[i][j] -= L[i][k] * U[k][j];
        }
        U[i][j] /= L[i][i];
      }
    }
  }

  // 求解Ly=b（顺序）
  for (k = 0; k < n; k++) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < i; j++) {
        dst[i * n + k] -= L[i][j] * dst[j * n + k];
      }
      dst[i * n + k] /= L[i][i];
    }
  }

  // 求解Ux=y（逆序）
  for (k = 0; k < n; k++) {
    for (i = n - 1; i >= 0; i--) {
      for (j = n - 1; j > i; j--) {
        dst[i * n + k] -= U[i][j] * dst[j * n + k];
      }
    }
  }
}

// int main() {

//   double src[9] = {1, 1, 1, 2, 3, 2, 3, 8, 2};
//   double dst[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

//   __TEST3_opt3_inverse_matrix_n(src, dst, 3);

//   printf("输出结果\n");
//   for (int i = 0; i < 9; i++) {
//     printf("%f ", dst[i]);
//     if (i % 3 == 2) {
//       printf("\n");
//     }
//   }
//   printf("*******************\n");

//   return 0;
// }
