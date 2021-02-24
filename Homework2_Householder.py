"""
Householder 变换实现 QR 分解 （对矩阵 A 不停地左乘以适当的 Householder 矩阵即可）
"""
import numpy as np
from math import *
import copy
import time


# 输入参数为一个任意的实矩阵 B
def Householder(B):
    A = copy.deepcopy(B)
    n = len(A)
    Q = np.eye(n, n)  # 初始化 Q 矩阵为单位阵
    # 对矩阵 A 的第 k 列的第 k 到 n 行进行 Householder 变换，i = k - 1
    for i in range(n-1):
        # 取出需要做变换的矢量（注意要做深拷贝），RK_x 是 n-k+1 维的，并将其变换为 Householder 矩阵中的 v 矢量
        x = copy.deepcopy(A[i:, i])
        x[0] -= sqrt(np.dot(x, x))
        # 计算 PA，计算过程中由于 P 的左上角是一个 k-1 维的单位阵，所以 A 只有右下角 n-k+1 维的元素发生变化
        # 将 A 这个 n-k+1 维子矩阵提取出来，并用 R_n-k+1 作用
        squareNorm = np.dot(x, x)
        for j in range(i, n):
            A[i:, j] = A[i:, j] - 2 * np.dot(x, np.dot(x, A[i:, j])) / squareNorm
        # 计算 QP，其中 P 是由一个 k-1 维单位阵和 n-k+1 维 Householder 矩阵构成
        # 在原来的 RK_x 前面添加 k-1 个 0，得到 Pk 的矢量表达形式
        for _ in range(i):
            x = np.insert(x, 0, 0)
        # 将 Q 的每一行与 Pk（矢量形式）进行乘积
        for m in range(n):
            Q[m, :] = Q[m, :] - 2 * np.dot(np.dot(Q[m, :], x), x) / squareNorm
    return Q, A


if __name__ == '__main__':
    A = np.array([[1.0, 2.0, 3.0, 4.0], [2.0, 3.0, 4.0, 5.0], [3.0, 4.0, 5.0, 6.0], [4.0, 5.0, 6.0, 7.0]])
    # A = np.random.rand(1200, 1200)
    begin = time.time()
    Q, R = Householder(A)
    end = time.time()
    print("正交矩阵 Q：\n" + str(Q) + "\n")
    print("上三角矩阵 R：\n" + str(R) + "\n")
    print("验证 Q 正交 Q.T*Q =：\n" + str(np.dot(Q, Q.T)) + "\n")
    print("验证 A = QR：\n" + str(np.dot(Q, R)) + "\n")
    print("用时：" + str(end-begin) + "s")


