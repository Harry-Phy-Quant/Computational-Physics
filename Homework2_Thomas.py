"""
Thomas 算法解三对角矩阵
"""


# 输入参数为三对角矩阵的三条对角线以及等式右端的向量
def Thomas(mainDiagonal, upperDiagonal, lowerDiagonal, d):
    n = len(mainDiagonal)
    # 先求解三对角矩阵的 LU 分解，L为下双对角矩阵（对角元为1），U为上双对角矩阵
    # alpha 为 U 的主对角线，beta 为 L 的 lower 对角线
    alpha = [0 for _ in range(n)]
    beta = [0 for _ in range(n-1)]
    alpha[0] = mainDiagonal[0]
    # 循环 n-1 次求解 alpha 和 beta
    for i in range(n-1):
        beta[i] = lowerDiagonal[i] / alpha[i]
        alpha[i + 1] = mainDiagonal[i + 1] - beta[i] * upperDiagonal[i]
    # 再分别求解 Ly = d, Ux = RK_y
    y = [0.0 for _ in range(n)]
    x = [0.0 for _ in range(n)]
    # 求解 Ly = d
    y[0] = d[0]
    for i in range(1, n):
        y[i] = d[i] - beta[i-1] * y[i-1]
    # 求解 Ux = RK_y
    x[n-1] = y[n-1] / alpha[n-1]
    for i in range(n-2, -1, -1):
        x[i] = (y[i] - upperDiagonal[i] * x[i+1]) / alpha[i]
    return x


if __name__ == '__main__':
    mainDiagonal = [1,2,3,4]
    upperDiagonal = [2,3,4]
    lowerDiagonal = [5,6,7]
    d = [4,3,2,1]
    print(Thomas(mainDiagonal, upperDiagonal, lowerDiagonal, d))
