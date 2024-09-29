from get_data import GRID
import os
import glob
from osgeo import gdal, ogr
import geopandas as gpd
import numpy as np
import rasterio as rio
import rasterio.mask


class cal_NDVI():
    def __init__(self, input_path, out_path):
        self.input_path = input_path
        self.out_path = out_path
        self.grids = GRID()

    # def add_mean_value_field(self, raster_file, shp_file, out_shp_file):
    #     """
    #     计算矢量图斑范围内的栅格均值并增加字段
    #     :param raster_file: 待计算均值的栅格底图
    #     :param shp_file: 与栅格数据叠加的矢量掩膜数据
    #     :param out_shp_file: 增加均值字段后输出的矢量数据，可与shp_file相同，即重写shp_file文件
    #     :return
    #     """
    #
    #     shp_data = gpd.GeoDataFrame.from_file(shp_file)
    #     raster_data = rio.open(raster_file)
    #     profile = raster_data.profile
    #     # 保存投影信息一致v
    #     shp_data = shp_data.to_crs(raster_data.crs)


    def run(self):
        file_name = glob.glob(os.path.join(self.input_path, '*.tif'))
        for name in file_name:
            code_name = name.split('\\')[-1].split('.')[0]
            new_name = os.path.join(self.out_path, code_name+'.tif')
            im_geotrans, im_proj, im_data, im_width, im_height = self.grids.read_img(name)
            NIR_band = im_data[3, :, :]
            Red_band = im_data[2, :, :]
            NDVI = (NIR_band-Red_band)/(NIR_band+Red_band+0.0000000001)#防止计算结果为0
            NDVI[NDVI > 1] = 1
            NDVI[NDVI < -1] = -1
            self.grids.write_img(new_name, im_proj, im_geotrans, NDVI)

if __name__ == '__main__':
    input_path = r"D:\DeepLearning_RS_data\KDXF_VHR\JINLIN-1\CGDZ8"
    out_path = r"D:\DeepLearning_RS_data\KDXF_VHR\JINLIN-1\CGDZ8\NDVI"
    obj = cal_NDVI(input_path, out_path)
    obj.run()



