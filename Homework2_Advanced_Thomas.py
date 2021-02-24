"""
推广的 Thomas 算法解含边界条件的三对角矩阵
"""
from Homework2_Thomas import Thomas


# 输入参数为三对角矩阵的三条对角线，右上角矩阵元，左下角矩阵元，以及等式右端的向量
def Advanced_Thomas(mainDiagonal, upperDiagonal, lowerDiagonal, upperRight, lowerLeft, d):
    # 先划去第一行第一列得到三对角子矩阵 A_tilda 的三角对角线
    A_tilda_mainDiagonal = mainDiagonal[1:]
    A_tilda_upperDiagonal = upperDiagonal[1:]
    A_tilda_lowerDiagonal = lowerDiagonal[1:]
    # 去除 d 的第一个元素，得到的向量记为 d_tilda
    d_tilda = d[1:]
    # 先用 Thomas 算法求解 A_tilda · u = d_tilda
    u = Thomas(A_tilda_mainDiagonal, A_tilda_upperDiagonal, A_tilda_lowerDiagonal, d_tilda)
    # 再用 Thomas 算法求解 A_tilda · v = psai
    psai = [0 for _ in range(len(d_tilda))]
    psai[0] = lowerDiagonal[0]
    psai[-1] = lowerLeft
    v = Thomas(A_tilda_mainDiagonal, A_tilda_upperDiagonal, A_tilda_lowerDiagonal, [-i for i in psai])
    # 下面计算 x1
    x = [0 for _ in range(len(d))]
    x[0] = (d[0] - upperDiagonal[0] * u[0] - upperRight * u[-1]) / (mainDiagonal[0] + upperDiagonal[0] * v[0] + upperRight * v[-1])
    for i in range(1, len(x)):
        x[i] = u[i-1] + x[0] * v[i-1]
    return x


if __name__ == '__main__':
    mainDiagonal = [1,2,3,4]
    upperDiagonal = [2,3,4]
    lowerDiagonal = [5,6,7]
    d = [4,3,2,1]
    print(Advanced_Thomas(mainDiagonal, upperDiagonal, lowerDiagonal, 1, 2, d))
