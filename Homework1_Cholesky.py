"""
Cholesky分解法解希尔伯特矩阵线性方程组
"""
from math import *


def cholesky(n):  # 矩阵的维数为 n
    # 通过二维列表的形式创建 n 维 Hilbert 矩阵
    H = [[1 / (i + j - 1) for i in range(1, n + 1)] for j in range(1, n + 1)]
    # b 为方程右端的常矢量
    b = [1 for _ in range(n)]
    # H = A\dagger A，通过数学归纳法求下三角矩阵 A\dagger
    for i in range(1, n + 1):  # 对 1 阶 ~ n 阶主子矩阵进行遍历，已知第 i-1 个主子阵，现在求第 i 个主子阵
        for j in range(1, i):  # 用反代法求新补充的一行的元素
            for k in range(1, j):
                H[i - 1][j - 1] -= H[i - 1][k - 1] * H[j - 1][k - 1]
            H[i - 1][j - 1] /= H[j - 1][j - 1]
            H[j - 1][i - 1] = 0  # 将右上角部分填充为 0
            H[i - 1][i - 1] -= (H[i - 1][j - 1]) ** 2  # 计算右下角元素（\beta)
        H[i - 1][i - 1] = sqrt(H[i - 1][i - 1])
    # 接下来分两步求解 RK_x
    # 第一步：求解 RK_y 使得 A\dagger RK_y = b
    y = [0 for _ in range(n)]  # 先初始化 RK_y
    for i in range(1, n + 1):  # 需要从上往下正代 n 次
        y[i - 1] = b[i - 1]
        for k in range(1, i):
            y[i - 1] -= H[i - 1][k - 1] * y[k - 1]
        y[i - 1] /= H[i - 1][i - 1]
    # 第二步：求解 RK_x 使得 Ax = RK_y，其中 A 为上三角矩阵
    H = [[H[j][i] for j in range(n)] for i in range(n)]  # 转置
    x = [0 for _ in range(n)]  # 先初始化 RK_x
    for i in range(n, 0, -1):  # 需要从下往上反代 n 次
        x[i - 1] = y[i - 1]
        for k in range(n, i, -1):
            x[i - 1] -= H[i - 1][k - 1] * x[k - 1]
        x[i - 1] /= H[i - 1][i - 1]
    return x


if __name__ == '__main__':
    # 矩阵的维数为 n，从 1 到 10 进行循环
    for i in range(1, 11):
        print(cholesky(i))
