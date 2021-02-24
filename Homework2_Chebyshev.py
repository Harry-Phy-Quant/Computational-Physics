"""
Chebyshev 近似
"""
from math import *
from prettytable import PrettyTable
import matplotlib.pyplot as plt


# 输入两个参数： y_list 为 N 阶 Chebyshev 多项式在零点的函数值，RK_x 为自变量
def Chebyshev(y_list, x):
    N = len(y_list)  # 阶数
    c = [0 for _ in range(N)]  # 初始化系数向量
    c[0] = (1 / N) * sum(y_list)
    for m in range(1, N):  # 循环求系数
        for k in range(N):
            c[m] += (2 / N) * cos(m*pi*(k+0.5)/N) * y_list[k]
    # 下面根据系数与对应的 Chebyshev 多项式求近似值
    # 先用迭代关系求 Tm(RK_x) 的值, m = 0 ~ N-1, 其中 T0 = 1，T1 = RK_x
    T = [1, x]
    for i in range(N-2):
        T.append(2 * x * T[-1] - T[-2])
    # 再与系数组合得到结果
    S = 0
    for m in range(N):
        S += c[m] * T[m]
    return S


if __name__ == '__main__':
    # 以 Runge 现象为例
    def Runge(x):
        return 1 / (1 + 25 * (x**2))

    support_x_list = [cos(pi*(k+0.5)/20) for k in range(20)]
    support_y_list = [Runge(i) for i in support_x_list]
    # 选取的测试点是切比雪夫节点和它们的中点
    x_list = []
    for i in range(len(support_x_list)-1):
        insert = (support_x_list[i] + support_x_list[i+1]) / 2
        x_list.append(support_x_list[i])
        x_list.append(insert)
    x_list.append(support_x_list[-1])
    y_list = [Runge(i) for i in x_list]
    Chebyshev_list = [Chebyshev(support_y_list, i) for i in x_list]
    difference = [abs(Chebyshev_list[i] - y_list[i]) for i in range(len(y_list))]

    tb = PrettyTable()
    tb.add_column('No.', [i for i in range(1, 40)])
    tb.add_column('RK_x', x_list)
    tb.add_column('f(RK_x)', y_list)
    tb.add_column('S20(RK_x)', Chebyshev_list)
    tb.add_column('|S20(RK_x) - f(RK_x)|', difference)
    print(tb)

    plt.plot(x_list, y_list, 'b.-', x_list, Chebyshev_list, 'r.:')
    plt.ylim((0, 1.2))
    plt.title('Chebyshev Polynomial')
    plt.legend(['  f(RK_x)', r'$S_{20}(RK_x)$'], loc=1)
    plt.show()

    """
    # 方便将输出结果直接存入 Latex 表格
    for i in range(len(x_list)):
        print(str(i+1)+"&"+format(x_list[i], '0.4f')+"&"+format(y_list[i], '0.4f')+"&"+format(Chebyshev_list[i], '0.4f')+"&"+format(difference[i], '0.4f')+"\\"+"\\")
    """
