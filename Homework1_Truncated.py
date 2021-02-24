"""
4(a)、4(b)所用代码，暴力法
"""

from math import *
import matplotlib.pyplot as plt

Lambda_list = []
result_list = []

# Lambda = 10~90，以 10 为间隔
for Lambda in range(10, 100, 10):
    Lambda_list.append(Lambda)
    # 计算级数和
    sum = 0
    maxK = Lambda
    for k in range(-maxK, maxK + 1):
        # 给定 k 以后，m 能取的值实际上可以限定在 -maxM 到 maxM中间
        maxM = int(sqrt(Lambda ** 2 - k ** 2))
        for m in range(-maxM, maxM + 1):
            # 给定 k, m 以后，N 能取的值实际上可以限定在 -maxN 到 maxN中间
            maxN = int(sqrt(Lambda ** 2 - k ** 2 - m ** 2))
            for n in range(-maxN, maxN + 1):
                sum += 1 / (k ** 2 + m ** 2 + n ** 2 - 0.5)
    # 计算积分
    integral = pow(2, 0.5) * pi * log(abs((pow(2, 0.5) - 2 * Lambda) / (pow(2, 0.5) + 2 * Lambda))) + 4 * pi * Lambda
    # 计算级数与积分的差
    result = sum - integral
    result_list.append(result)
    print("Lambda= " + str(Lambda) + " completed! Result= " + str(result))
    # 如果已经达到精度要求，则停止循环
    if Lambda != 10 and (abs(result_list[-1] - result_list[-2])) <= 0.00001:
        break

# Lambda = 100~370，以 30 为间隔
for Lambda in range(100, 400, 30):
    Lambda_list.append(Lambda)
    # 计算级数和
    sum = 0
    maxK = Lambda
    for k in range(-maxK, maxK + 1):
        maxM = int(sqrt(Lambda ** 2 - k ** 2))
        for m in range(-maxM, maxM + 1):
            maxN = int(sqrt(Lambda ** 2 - k ** 2 - m ** 2))
            for n in range(-maxN, maxN + 1):
                sum += 1 / (k ** 2 + m ** 2 + n ** 2 - 0.5)
    # 计算积分
    integral = pow(2, 0.5) * pi * log(abs((pow(2, 0.5) - 2 * Lambda) / (pow(2, 0.5) + 2 * Lambda))) + 4 * pi * Lambda
    # 计算级数与积分的差
    result = sum - integral
    result_list.append(result)
    print("Lambda= " + str(Lambda) + " completed! Result= " + str(result))
    if abs(result_list[-1] - result_list[-2]) <= 0.00001:
        break

# Lambda = 400~1200，以 80 为间隔
for Lambda in range(400, 1280, 80):
    Lambda_list.append(Lambda)
    # 计算级数和
    sum = 0
    maxK = Lambda
    for k in range(-maxK, maxK + 1):
        maxM = int(sqrt(Lambda ** 2 - k ** 2))
        for m in range(-maxM, maxM + 1):
            maxN = int(sqrt(Lambda ** 2 - k ** 2 - m ** 2))
            for n in range(-maxN, maxN + 1):
                sum += 1 / (k ** 2 + m ** 2 + n ** 2 - 0.5)
    # 计算积分
    integral = pow(2, 0.5) * pi * log(abs((pow(2, 0.5) - 2 * Lambda) / (pow(2, 0.5) + 2 * Lambda))) + 4 * pi * Lambda
    # 计算级数与积分的差
    result = sum - integral
    result_list.append(result)
    print("Lambda= " + str(Lambda) + " completed! Result= " + str(result))
    if abs(result_list[-1] - result_list[-2]) <= 0.00001:
        break

# Lambda = 1300~1500，以 100 为间隔
for Lambda in range(1300, 1600, 100):
    Lambda_list.append(Lambda)
    # 计算级数和
    sum = 0
    maxK = Lambda
    for k in range(-maxK, maxK + 1):
        maxM = int(sqrt(Lambda ** 2 - k ** 2))
        for m in range(-maxM, maxM + 1):
            maxN = int(sqrt(Lambda ** 2 - k ** 2 - m ** 2))
            for n in range(-maxN, maxN + 1):
                sum += 1 / (k ** 2 + m ** 2 + n ** 2 - 0.5)
    # 计算积分
    integral = pow(2, 0.5) * pi * log(abs((pow(2, 0.5) - 2 * Lambda) / (pow(2, 0.5) + 2 * Lambda))) + 4 * pi * Lambda
    # 计算级数与积分的差
    result = sum - integral
    result_list.append(result)
    print("Lambda= " + str(Lambda) + " completed! Result= " + str(result))
    if abs(result_list[-1] - result_list[-2]) <= 0.00001:
        break

# 作图
x_data = Lambda_list
y_data = result_list
plt.plot(x_data, y_data)
plt.show()
