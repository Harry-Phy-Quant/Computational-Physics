"""
4(c)所用代码，优化版
"""

from math import *
import scipy.integrate as integrate
import time

start = time.time()

# 计算第一项
first_result = []
for Lambda in range(10, 50, 10):
    first = 0
    for k in range(-Lambda, Lambda + 1):
        for m in range(-Lambda, Lambda + 1):
            for n in range(-Lambda, Lambda + 1):
                    first += exp(-(k ** 2 + m ** 2 + n ** 2 - 0.5)) / (k ** 2 + m ** 2 + n ** 2 - 0.5)
    first_result.append(first)
    # 若计算结果精度已经到达 10^-10 则终止循环
    if len(first_result) >= 2 and abs(first_result[-1] - first_result[-2]) <= 10 ** (-10):
        print("第一项的结果已达到 10^-10 精度！")
        break
# 第一项的结果
first = first_result[-1]


# 计算第二项
second_result = []
for Lambda in range(10, 50, 10):
    second = 0
    for k in range(-Lambda, Lambda + 1):
        for m in range(-Lambda, Lambda + 1):
            for n in range(-Lambda, Lambda + 1):
                if (k ** 2 + m ** 2 + n ** 2) != 0:
                    # 计算积分，代码可能会产生高亮，应该是版本冲突问题，不影响计算结果！！！
                    second += integrate.quad(lambda t: exp(0.5 * t) * ((pi / t) ** 1.5) * exp(-pi * pi * (k ** 2 + m ** 2 + n ** 2) / t), 0, 1)[0]
    second_result.append(second)
    # 若计算结果精度已经到达 10^-10 则终止循环
    if len(second_result) >= 2 and abs(second_result[-1] - second_result[-2]) <= 10 ** (-10):
        print("第二项的结果已达到 10^-10 精度！")
        break
# 第二项的结果
second = second_result[-1]


# 计算第三项和第四项
# 计算积分，代码可能会产生高亮，应该是版本冲突问题，不影响计算结果！！！
third = integrate.quad(lambda t: (exp(0.5 * t) - 1) * ((pi / t) ** 1.5), 0, 1)[0]
fourth = - 2 * (pi ** 1.5)

end = time.time()

print("最终结果为： " + str(first + second + third + fourth))
print("计算耗时为： " + str(end - start) + "s")
