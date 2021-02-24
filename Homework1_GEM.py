"""
高斯消元法解希尔伯特矩阵线性方程组
"""


def GEM(n):  # 矩阵的维数为 n
    # 通过二维列表的形式创建 n 维 Hilbert 矩阵
    H = [[1 / (i + j - 1) for i in range(1, n + 1)] for j in range(1, n + 1)]
    # b 为方程右端的常矢量
    b = [1 for _ in range(n)]
    # 对 Hilbert 矩阵消元，得到上三角矩阵
    for k in range(1, n):  # 最外层循环，需要进行 n-1 次消元过程
        for i in range(k + 1, n + 1):  # 对需要消元的每行进行遍历
            l = -H[i - 1][k - 1] / H[k - 1][k - 1]  # 计算该行的倍数因子
            for j in range(k, n + 1):  # 对每列进行遍历
                H[i - 1][j - 1] += l * H[k - 1][j - 1]  # 更新矩阵元素
            b[i - 1] += l * b[k - 1]  # 更新右端项
    # 反向迭代，得到待求矢量 RK_x
    x = [0 for _ in range(n)]  # 先初始化结果矢量
    for i in range(n, 0, -1):  # 需要从下往上反代 n 次
        x[i - 1] = b[i - 1]
        for k in range(n, i, -1):
            x[i - 1] -= H[i - 1][k - 1] * x[k - 1]
        x[i - 1] /= H[i - 1][i - 1]
    return x


if __name__ == '__main__':
    # 矩阵的维数为 n，从 1 到 10 进行循环
    for i in range(1, 11):
        print(GEM(i))
