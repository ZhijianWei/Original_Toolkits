import pandas as pd
import numpy as np

file_path = r"D:\Wei Zhijian\Wheat_prod_report\两列数据找相关关系.xlsx"
df = pd.read_excel(file_path)

#键入列名
data1 = df['ref']
data2 = df['asd']

correlation_coefficient = np.corrcoef(data1, data2)[0, 1]
print(f"皮尔逊相关系数为: {correlation_coefficient:.2f}")
r_squared = correlation_coefficient**2
print(f"R^2值为: {r_squared:.2f}")



if 0.7 <= correlation_coefficient <= 1:
    print("两列数据高度正相关")
elif -1 <= correlation_coefficient <= -0.7:
    print("两列数据高度负相关")
elif 0.3 <= correlation_coefficient < 0.7:
    print("两列数据中度正相关")
elif -0.7 < correlation_coefficient <= -0.3:
    print("两列数据中度负相关")
else:
    print("两列数据的相关性较低")

# import pandas as pd
# from sklearn.preprocessing import PolynomialFeatures
#
# # 加载数据
# file_path = r'C:\Users\XavierWynne\Desktop\INRday_Ura.csv'
# df = pd.read_csv(file_path)
# data = df[['dose', 'Uridine']]
#
# # 使用 PolynomialFeatures 生成三次多项式特征
# poly = PolynomialFeatures(degree=2, include_bias=False)
# poly_data = poly.fit_transform(data)
# poly_df = pd.DataFrame(poly_data, columns=poly.get_feature_names_out())
#
# # 获取交互特征列
# interaction_feature_name = 'dose Uridine'  # 这是基于你之前给出的列名称。确保这是你想要的交互特征名。
# interaction_feature = poly_df[interaction_feature_name]
#
# # 保存交互特征列为Excel文件
# output_path = r'C:\Users\XavierWynne\Desktop\interaction_feature.xlsx'
# interaction_feature.to_excel(output_path, index=False, header=True)
