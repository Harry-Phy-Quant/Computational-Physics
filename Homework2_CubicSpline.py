"""
等间距情况下的三次样条函数插值
"""
from Homework2_Thomas import Thomas
from prettytable import PrettyTable
import matplotlib.pyplot as plt


# 输入五个参数，x_list, y_list 为一一对应的支撑点列表，m0, mn 为两个端点处的二阶导数值（边界条件），RK_x 为自变量的值
def CubicSpline(x_list, y_list, m0, mn, x):
    N = len(x_list) - 1  # 区间数
    h = x_list[1] - x_list[0]  # 两点间隔
    d = [0 for _ in range(N + 1)]
    # 由边界条件确定 d[0]，d[N]
    d[0] = 2 * m0
    d[N] = 2 * mn
    # 计算 d[1] ~ d[N-1]
    for i in range(1, N):
        d[i] = 3 * (y_list[i-1] + y_list[i+1] - 2 * y_list[i]) / (h**2)
    # 给出 λ 和 μ 向量
    lambdaList = [0.5 for _ in range(N)]
    lambdaList[0] = 0
    miuList = [0.5 for _ in range(N)]
    miuList[N - 1] = 0
    mainDiagonal = [2 for _ in range(N + 1)]
    # 用 Thomas 算法解三对角矩阵得到 m 向量
    m = Thomas(mainDiagonal, lambdaList, miuList, d)
    # 用得到的 m 向量来计算任意给定 RK_x 值的预测值
    i = int((x - x_list[0]) / h)  # 确定 RK_x 所属的区间范围 RK_x[i] ~ RK_x[i+1]
    if i == N:
        return y_list[-1]
    else:
        A = ((y_list[i + 1] - y_list[i]) / h) - (h * (m[i + 1] - m[i]) / 6)
        B = y_list[i] - m[i] * (h ** 2) / 6
        S = - m[i] * ((x - x_list[i + 1]) ** 3) / (6 * h) + m[i + 1] * ((x - x_list[i]) ** 3) / (6 * h) + A * (
                    x - x_list[i]) + B
        return S


if __name__ == '__main__':
    # 以 Runge 现象为例
    def Runge(x):
        return 1 / (1 + 25 * (x ** 2))

    support_x_list = [-1 + 0.1 * i for i in range(21)]
    support_y_list = [Runge(i) for i in support_x_list]
    x_list = [-1 + 0.05 * i for i in range(41)]
    y_list = [Runge(i) for i in x_list]
    CubicSpline_list = [CubicSpline(support_x_list, support_y_list, 925/4394, 925/4394, i) for i in x_list]
    difference = [abs(CubicSpline_list[i] - y_list[i]) for i in range(len(y_list))]

    tb = PrettyTable()
    tb.add_column('No.', [i for i in range(1, 42)])
    tb.add_column('RK_x', x_list)
    tb.add_column('f(RK_x)', y_list)
    tb.add_column('S(RK_x)', CubicSpline_list)
    tb.add_column('|S(RK_x) - f(RK_x)|', difference)
    print(tb)

    plt.plot(x_list, y_list, 'b.-', x_list, CubicSpline_list, 'r.:')
    plt.ylim((0, 1.2))
    plt.title('Cubic Spline')
    plt.legend(['  f(RK_x)', r'$S(RK_x)$'], loc=1)
    plt.show()

    '''
    # 方便将输出结果直接存入 Latex 表格
    for i in range(len(x_list)):
        print(str(i+1)+"&"+format(x_list[i], '0.4f')+"&"+format(y_list[i], '0.4f')+"&"+format(CubicSpline_list[i], '0.4f')+"&"+format(difference[i], '0.4f')+"\\"+"\\")
    '''

