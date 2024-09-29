from get_data import GRID
import os
import glob
from osgeo import gdal, ogr
import geopandas as gpd
import numpy as np
import rasterio as rio
import rasterio.mask


def add_mean_value_field(raster_file, shp_file, out_shp_file):
    """
    计算矢量图斑范围内的栅格均值并增加字段
    :param raster_file: 待计算均值的栅格底图
    :param shp_file: 与栅格数据叠加的矢量掩膜数据
    :param out_shp_file: 增加均值字段后输出的矢量数据，可与shp_file相同，即重写shp_file文件
    :return
    """

    shp_data = gpd.GeoDataFrame.from_file(shp_file)
    raster_data = rio.open(raster_file)
    profile = raster_data.profile
    # 保存投影信息一致
    shp_data = shp_data.to_crs(raster_data.crs)

    out_shp_data = shp_data.copy()
    mean_value_field = []
    for i in range(0, len(shp_data)):
        # 获取矢量数据的features
        geo = shp_data.geometry[i]
        feature = [geo.__geo_interface__]
        # 通过feature裁剪栅格影像
        out_image, out_transform = rasterio.mask.mask(raster_data, feature, all_touched=True, crop=True,
                                                      nodata=raster_data.nodata)
        # 获取影像Value值，并转化为list
        out_list = out_image.data.tolist()
        # 除去list中的Nodata值
        out_list = out_list[0]
        out_data = []
        for k in range(len(out_list)):
            for j in range(len(out_list[k])):
                if out_list[k][j] > 0:
                    out_data.append(out_list[k][j])
        # 求数据的平均数
        if len(out_data):
            mean_data = np.mean(out_data)
        else:
            mean_data = None
        mean_value_field.append(mean_data)

    # 增加属性字段，并将GeodataFrame导出为shp文件
    out_shp_data.insert(out_shp_data.shape[1], 'des', mean_value_field)
    out_shp_data.to_file(out_shp_file)

raster_file = r"D:\SR\test0917\output\Des_0\temp\step1.tif"
shp_file = r"D:\SR\test0917\output\Des_0\temp\res2.shp"
out_shp_file = r"D:\SR\test0917\output\Des_0\temp\res3.shp"
add_mean_value_field(raster_file, shp_file, out_shp_file)