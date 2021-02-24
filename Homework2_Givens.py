"""
Givens 变换实现 QR 分解 （对矩阵 A 不停地左乘以适当的 Givens 矩阵即可）
"""
import numpy as np
from math import *
import copy
import time


# 输入参数为一个任意的实矩阵 B
def Givens(B):
    A = copy.deepcopy(B)
    n = len(A)
    Q = np.eye(n, n)  # 初始化 Q 矩阵为单位阵
    # 对矩阵 A_pp, A_pq, A_qp, A_qq 组成的 2×2 子矩阵进行 Givens 旋转，使得 A_qp = 0
    for p in range(n-1):
        for q in range(p+1, n):
            # 计算 Givens 矩阵的变换元而不构造矩阵本身
            c = A[p, p] / sqrt(A[p, p] ** 2 + A[q, p] ** 2)
            s = A[q, p] / sqrt(A[p, p] ** 2 + A[q, p] ** 2)
            # 将 Givens 矩阵作用在 A 左边
            # 取出 A 的第 p 行和第 q 行的第 p 列至第 n 列
            A_p = copy.deepcopy(A[p, p:])
            A_q = copy.deepcopy(A[q, p:])
            A[p, p:] = c * A_p + s * A_q
            A[q, p:] = -s * A_p + c * A_q
            # 将 Givens 矩阵作用在 Q 右边
            # 取出 Q 的第 p 列和第 q 列
            Q_p = copy.deepcopy(Q[:, p])
            Q_q = copy.deepcopy(Q[:, q])
            Q[:, p] = c * Q_p + s * Q_q
            Q[:, q] = -s * Q_p + c * Q_q
    return Q, A


if __name__ == '__main__':
    A = np.array([[1.0, 2.0, 3.0, 4.0], [2.0, 3.0, 4.0, 5.0], [3.0, 4.0, 5.0, 6.0], [4.0, 5.0, 6.0, 7.0]])
    # A = np.random.rand(1200, 1200)
    begin = time.time()
    Q, R = Givens(A)
    end = time.time()
    print("正交矩阵 Q：\n" + str(Q) + "\n")
    print("上三角矩阵 R：\n" + str(R) + "\n")
    print("验证 Q 正交 Q.T*Q =：\n" + str(np.dot(Q, Q.T)) + "\n")
    print("验证 A = QR：\n" + str(np.dot(Q, R)) + "\n")
    print("用时：" + str(end-begin) + "s")


