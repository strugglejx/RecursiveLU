
#define __USE_GNU
#define _GNU_SOURCE

#include "cal/LU_C.c"
// #include "TEST5.c"
#include "cal/LU_C_Optimize.c"
#include "cal/LU_ASM_1.c"
#include "cal/LU_ASM_2.c"
#include "cal/LU_ASM_3.c"
#include "test.h"
#include <getopt.h>
#include <malloc.h>
#include <pthread.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> // POSIX操作系统

#define WARM_UP_COUNT 5

#define app_panic(format, args...) \
  do                               \
  {                                \
    printf(format, ##args);        \
    abort();                       \
  } while (0)

static void usage(void)
{
  printf("使用参数说明:./test -o 4 -n 1\n"
         "-o\t输入测试矩阵的阶\n"
         "-n\t输入测试矩阵的数量\n"
         "-h\t帮助信息\n");
}
static const char short_options[] = "o:n:h";
static const struct option long_options[] = {
    {"order", required_argument, NULL, 'o'},
    {"count", required_argument, NULL, 'n'},
    {"help", no_argument, NULL, 'h'},
    {0, 0, 0, 0},
};

static void parameters_parser(int argc, char *argv[])
{

  if (argc < 2)
  {
    printf("[Error]: too few input parameters\n");
    usage();
    exit(-1);
  }
  ORDER = 4;
  COUNT = 1;
  for (;;)
  {
    int index, c;
    c = getopt_long(argc, argv, short_options, long_options, &index);
    if (-1 == c)
      break;

    switch (c)
    {
    case 'o':
      ORDER = atoi(optarg);
      break;
    case 'n':
      COUNT = atoi(optarg);
      break;
    case 'h':
    default:
      usage();
      exit(-1);
    }
  }
}

static void run(struct algorithm_struct *pt, double *src, int n, int warm_up)
{
  double *dst = (double *)malloc(n * n * sizeof(double));
  if (dst == NULL)
  {
    printf("[ERROR]: dst malloc error\n");
    exit(-1);
  }

  gettimeofday(&(pt->start), NULL);
  __asm__ __volatile__("" ::
                           : "memory");
  pt->test(src, dst, n);
  __asm__ __volatile__("" ::
                           : "memory");
  gettimeofday(&(pt->end), NULL);

  if (warm_up == 0)
  {
    pt->sum += (pt->end.tv_sec - pt->start.tv_sec) * 1000000 +
                       pt->end.tv_usec - pt->start.tv_usec;
    pt->cnt++;

    double *mul = (double *)malloc(n * n * sizeof(double));
    double error = 0;
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {

        if (j <= i)
        {
          double res = dst[i * n + j];
          for (int k = 0; k < j; k++)
          {
            res += dst[i * n + k] * dst[k * n + j];
          }
          error += fabs(src[i * n + j] - res);
        }
        else
        {
          double res = 0;
          for (int k = 0; k <= i; k++)
          {
            res += dst[i * n + k] * dst[k * n + j];
          }
          error += fabs(src[i * n + j] - res);
        }
      }
    }
    pt->error += error / (COUNT * 1.0) / ORDER / ORDER;
  }
  else if (warm_up == 1)
  {
    return;
  }

  //检测LU分解运算是否正确

  if (COUNT == 1)
  {
  double *mul = (double *)malloc(n * n * sizeof(double));
  double error = 0;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {

      if (j <= i)
      {
        double res = dst[i * n + j];
        for (int k = 0; k < j; k++)
        {
          res += dst[i * n + k] * dst[k * n + j];
        }
        error += fabs(src[i * n + j] - res);
      }
      else
      {
        double res = 0;
        for (int k = 0; k <= i; k++)
        {
          res += dst[i * n + k] * dst[k * n + j];
        }
        error += fabs(src[i * n + j] - res);
      }
    }
  }
  pt->error += error / (COUNT * 1.0) / n / n;

    printf("Accumulated error of %s: %f\n", pt->name, fabs(error));
    printf("**************************\n\n");
  }

  free(dst);
}

static inline void get_cpu_mask(pid_t pid, cpu_set_t *mask)
{
  if (sched_getaffinity(pid, sizeof(cpu_set_t), mask) == -1)
  {
    app_panic("Get cpu affinity failed.\n");
  }
}

static inline void set_cpu_mask(pid_t pid, cpu_set_t *mask)
{
  if (sched_setaffinity(pid, sizeof(cpu_set_t), mask) == -1)
  {
    app_panic("Set cpu affinity failed.\n");
  }
}

static void thread_bind(int cpu)
{
  cpu_set_t cpu_mask;
  get_cpu_mask(0, &cpu_mask);
  CPU_ZERO(&cpu_mask);
  CPU_SET(cpu, &cpu_mask);
  set_cpu_mask(0, &cpu_mask);
}

int main(int argc, char *argv[])
{
  struct algorithm_struct *pt;
  double *src, *psrc;
  parameters_parser(argc, argv);

  thread_bind(0);

  src = (double *)malloc(COUNT * ORDER * ORDER * sizeof(double));
  if (src == NULL)
  {
    printf("[ERROR]: src malloc error\n");
    exit(-1);
  }
  srand((unsigned)time(NULL));
  for (int i = 0; i < COUNT; i++)
  {
    psrc = src + i * ORDER * ORDER;
    for (int j = 0; j < ORDER * ORDER; j++)
      psrc[j] = 1.2 * (j + 1) + i * 5.6 + (rand() % 20);
  }

  for (int i = 0; i < COUNT + WARM_UP_COUNT; i++)
  {
    //  warm up
    if (i < WARM_UP_COUNT)
    {
      for_each_test(pt)
      {
        psrc = src + (rand() % COUNT) * ORDER * ORDER;
        run(pt, psrc, ORDER, 1);
      }
    }
    else
    {
      for_each_test(pt)
      {
        psrc = src + (rand() % COUNT) * ORDER * ORDER;
        run(pt, psrc, ORDER, 0);
      }
    }
  }

  // for_each_test(pt) printf("%s:\t%ld(us)\n", pt->name, pt->sum);

  for_each_test(pt) printf("%s:\t%ld(us);error:\t%.16lf\n", pt->name, pt->sum, pt->error);

  free(src);
  return 0;
}
