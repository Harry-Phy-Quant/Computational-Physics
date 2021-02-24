"""
幂次法求矩阵最大模的本征值和本征矢
"""
import numpy as np
from math import *


# 输入参数为符合幂次法要求的矩阵 A，初始化向量 q，最大迭代次数 maxN
def powerIteration(A, q, maxN):
    for i in range(maxN):
        z = np.dot(A, q)
        q = z / sqrt(np.dot(z, z))
    v = np.dot(q, np.dot(A, q))
    return v, q


if __name__ == '__main__':
    # 构造矩阵 A
    n = 10
    A = np.eye(n, n)
    for i in range(n):
        A[i, i] = 2
    for i in range(n-1):
        A[i, i + 1] = -1
        A[i + 1, i] = -1
    # 添加边界条件
    A[0, n-1] = -1
    A[n-1, 0] = -1

    # 随机生成一个 q 作为初始化向量
    q0 = np.random.rand(n)
    v, q = powerIteration(A, q0, 1000)
    print("迭代1000次，最大本征值为：\n" + str(v) + "\n")
    print("迭代1000次，最大本征值对应的本征矢量为：\n" + str(q))

