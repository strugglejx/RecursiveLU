# RecursiveLU
## 项目介绍
>>RecursiveLU是基于ARMv8处理器架构的SIMD汇编优化程序，实现了现有递归LU 分解算法的SIMD汇编优化，并达到了1.5 倍的加速效果。
## 目录结构
├── src <br>
│     ├── main.c <br>
│     ├── test.h <br>
└     └── cal <br>
        ├── LU_ASM_1.c <br>
        ├── LU_ASM_2.c <br>
        ├── LU_ASM_3.c <br>
        ├── LU_C_Optimize.c <br>
        └── LU_C.c <br>
## 编译
''' Shell
gcc -O3 main.c -o main
'''
## 实验结果
>>以下结果在以下参数下得到： 矩阵大小为512 × 512，分块矩阵的大小为 4 × 64。
   算法  | 时间（微秒/100 轮）
 ---- | ----- | ------  
 递归分块的 LU 分解（C 语言）  |  3839746
 递归分块的 LU 分解（内联汇编）|  3266809
 递归分块的 LU 分解+双重分块策略优化（C 语言）| 6464961
 递归分块的 LU 分解+双重分块策略优化（内联汇编）| 2598728
