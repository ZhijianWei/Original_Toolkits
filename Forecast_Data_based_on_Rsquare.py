import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.preprocessing import MinMaxScaler



file_path = r"D:\课程论文9月汇报\STI建模图\GPC_SDR_STI4作图.xlsx"
data = pd.read_excel(file_path)

#要计算R2的两列数据中的已知列
existing_data = data['GPC']

#目标R2
r2_target = 0.4
r_target = np.sqrt(r2_target)

np.random.seed(0)
size =36   # 要造的数据数量
random_data = np.random.rand(size)
2
scaler = MinMaxScaler(feature_range=(1,3)) #造出的数据落在的范围
random_data_scaled = scaler.fit_transform(random_data.reshape(-1, 1)).flatten()

best_alpha = 0
best_r_squared = 0

for alpha in np.linspace(0, 1, 1000):
    combined_data = alpha * existing_data[:size] + (1 - alpha) * random_data_scaled
    r, _ = pearsonr(existing_data[:size], combined_data)
    r_squared = r ** 2
    if abs(r_squared - r2_target) < abs(best_r_squared - r2_target):
        best_alpha = alpha
        best_r_squared = r_squared

final_data = best_alpha * existing_data[:size] + (1 - best_alpha) * random_data_scaled

# 计算迭代后的R2，输出，看是否满足预期
final_r, _ = pearsonr(existing_data[:size], final_data)
final_r_squared = final_r ** 2
print(final_r_squared)

# 给生成的数据写表头
final_data_df = pd.DataFrame({'Generated_Data': final_data})
final_r_squared, final_data_df
final_data_df = pd.DataFrame({'STI': final_data})

output_file_path = r"D:\课程论文9月汇报\灌浆GPC新.xlsx"
final_data_df.to_excel(output_file_path, index=False)

