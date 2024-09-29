import pandas as pd
import numpy as np


class hyperspectralTomultispectral:
    def __init__(self, input_hyperspectral, SRF, output_multispectral):
        """
        :param input_hyperspectral:  输入的高光谱文件 .xlsx
        :param SRF:  传感器光谱响应函数文件 .xlsx
        :param output_multispectral: 输入的多光谱文件 .xlsx
        @ created by Sunhaoran, 20240507, From NJAU
        """

        self.input_hyperspectral = input_hyperspectral
        self.SRF = SRF
        self.output_multispectral = output_multispectral
        self.read_hyperspectral()

    def read_hyperspectral(self):
        self.Hyperspectral_df = pd.read_csv(self.input_hyperspectral)
        self.station_name = self.Hyperspectral_df.columns[1:]
        self.Hyperspectral_data = self.Hyperspectral_df.values[:, 1:]
        self.srf_df = pd.read_csv(self.SRF)
        self.srf_data = self.srf_df.values[:, 1:]
        self.band_name = self.srf_df.index[1:]
        self.c1 = self.Hyperspectral_data.shape[1]
        self.c2 = self.srf_data.shape[1]

        self.srf_data[np.isnan(self.srf_data)] = 0
        self.srf_data[self.srf_data < 0] = 0




    def cal_erf(self):
        self.erf = np.zeros((self.c2, self.c1) )# 创建存储等效Rrs的矩阵
        for i in range(self.c2):
            for j in range(self.c1):
                arr = np.where(self.srf_data>0)                    # 积分下限
                l1 = arr[0]
                l2 = arr[-1]         # 积分上限
                self.Spectral_integration = np.trapz(self.Hyperspectral_data[l1:l2, i]*self.SRF[l1:l2, i])
                self.SRF_integration = np.trapz(SRF[l1:l2, i])
                self.erf = self.Spectral_integration/self.SRF_integration
        print(self.erf)



if __name__ == '__main__':
    H = r'D:\Lidong_data\HtoM\ASD.csv'
    SRF = r'D:\Lidong_data\HtoM\Sentinel2A_SRF.csv'
    Out = r"D:\Lidong_data\HtoM\out.csv"
    obj = hyperspectralTomultispectral(H, SRF, Out)
    obj.cal_erf()


