import pandas as pd
from sklearn.linear_model import LinearRegression
import numpy as np

df = pd.read_excel(r"D:\Wei Zhijian\Wheat_prod_report\两列数据找相关关系.xlsx")

column1 = df['ref'].values.reshape(-1, 1)
column2 = df['asd'].values.reshape(-1, 1)

# 创建线性回归模型
model = LinearRegression()
model.fit(column1, column2)

# 获取斜率a和截距b
slope = model.coef_[0][0]
intercept = model.intercept_[0]
print(f'Slope (a): {slope}')
print(f'Intercept (b): {intercept}')

# 构建数学变换公式(如果第一列数据经过线性变换可以表示为：y = ax + b)
print(f'The mathematical transformation is: Column2 = {slope} * Column1 + {intercept}')

