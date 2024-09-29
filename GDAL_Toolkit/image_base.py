from osgeo import gdal, ogr, osr
import numpy as np
import cv2
import tempfile
import numpy.ma as ma
import os

# =====栅格处理基本功能函数（类）=====
# 获取二值图像
def get_binary_pic(img_array, threshold=0, mask=255):
    """获取二值图像"""
    # img_array = GRID.get_data(img,1)
    # gray = cv2.cvtColor(np.array(img_array), cv2.COLOR_BGR2GRAY)
    gray = img_array
    _, thresh = cv2.threshold(gray, threshold, mask, cv2.THRESH_BINARY)
    return thresh
# 求取栅格最大面积矩形
def largestRectangleArea(heights):
    """直方图获取最大矩形"""
    n = len(heights)
    left, right = [0] * n, [n] * n

    mono_stack = list()
    for i in range(n):
        while mono_stack and heights[mono_stack[-1]] >= heights[i]:
            right[mono_stack[-1]] = i
            mono_stack.pop()
        left[i] = mono_stack[-1] if mono_stack else -1
        mono_stack.append(i)
    ans = 0
    y = 0
    w = 0
    h = 0
    for i in range(n):
        tmp = (right[i] - left[i] - 1) * heights[i]
        if tmp > ans:
            ans = tmp
            y = left[i] + 1
            w = int(heights[i])
            h = right[i] - left[i] - 1
    return ans, y, w, h
def maximalRectangle(matrix, mask=255):
    """二值矩阵获取最大矩形"""
    n, m = matrix.shape
    left = np.zeros([n, m])
    for i in range(n):
        for j in range(m):
            if matrix[i, j] == mask:
                if j == 0:
                    left[i, j] = 1
                else:
                    left[i, j] = left[i, j - 1] + 1
    ret = 0
    ans_x = 0
    ans_y = 0
    ans_w = 0
    ans_h = 0
    for i in range(m):
        tmp_ret, y, w, h = largestRectangleArea(left[:, i])
        if tmp_ret > ret:
            ret = tmp_ret
            ans_x = i - w + 1
            ans_y = y
            ans_w = w
            ans_h = h
    return ans_x, ans_y, ans_w, ans_h

def stretch_percentclip(data, percent=None, nodata=None, gamma=None):
    """
    param data: np 数组
    param percent: 百分比
    nodata: 空值
    gamma: 1.5(这里是越小越亮)
    """
    if percent is None: percent = 0
    if nodata is None: nodata = -9999999
    if gamma is None: gamma = 0.5
    data = ma.masked_equal(data, nodata)
    height, width = data.shape
    hist, bin_edges = np.histogram(data, bins=1000, density=False)
    cumhist = np.cumsum(hist) / (height * width)
    q_min = np.min(data)
    q_max = np.max(data)
    data = (1.0 / (q_max - q_min) * (data - q_min))
    data_gamma = (np.power(data, gamma) * 255)

    return data_gamma.astype('uint8')

class GRID:
    """
    @ Created by Sunhaoran 20230331
    """
    # 读图像文件
    @staticmethod
    def read_img(filename):
        dataset = gdal.Open(filename)  # 打开文件:Tiff或者ENVI文件

        im_width = dataset.RasterXSize  # 栅格矩阵的列数
        im_height = dataset.RasterYSize  # 栅格矩阵的行数

        im_geotrans = dataset.GetGeoTransform()  # 仿射矩阵
        im_proj = dataset.GetProjection()  # 地图投影信息
        im_data = dataset.ReadAsArray(0, 0, im_width, im_height)  # 将数据写成数组，对应栅格矩阵

        del dataset
        return im_proj, im_geotrans, im_data

    @staticmethod
    def get_data(filename, band):
        """
        band: 读取第几个通道的数据
        """
        dataset = gdal.Open(filename)
        band = dataset.GetRasterBand(band)
        data = band.ReadAsArray()
        return data

    # 获取经纬度信息
    @staticmethod
    def get_lon_lat_Transform( filename):
        """获取经纬度信息"""
        dataset = gdal.Open(filename)  # 打开文件:Tiff或者ENVI文件

        im_width = dataset.RasterXSize  # 栅格矩阵的列数
        im_height = dataset.RasterYSize  # 栅格矩阵的行数

        gtf = dataset.GetGeoTransform()  # 仿射矩阵
        x_range = range(0, im_width)
        y_range = range(0, im_height)
        x, y = np.meshgrid(x_range, y_range)
        lon = gtf[0] + x * gtf[1] + y * gtf[2]
        lat = gtf[3] + x * gtf[4] + y * gtf[5]
        return lon, lat, gtf

    # 获取影像四至范围
    @staticmethod
    def get_raster_extent(filename):
        """
        获取影像四至范围
        :param filename:
        :return:
        """
        dataset = gdal.Open(filename)  # 打开文件:Tiff或者ENVI文件

        im_width = dataset.RasterXSize  # 栅格矩阵的列数
        im_height = dataset.RasterYSize  # 栅格矩阵的行数

        gtf = dataset.GetGeoTransform()  # 仿射矩阵
        left_lon = gtf[0]
        up_lat = gtf[3]
        right_lon = gtf[0] + (im_width + 0.05) * gtf[1] + (im_height + 0.05) * gtf[2]
        low_lat = gtf[3] + (im_width + 0.05) * gtf[4] + (im_height + 0.05) * gtf[5]
        return left_lon, right_lon, low_lat, up_lat

    @staticmethod
    def get_raster_CGCS2000extent(filename):
        """
            获取影像2000系的四至范围
            :param filename:
            :return:
        """
        # change proj to epsg:4490
        if GRID.raster_proj(filename, "EPSG") != "4490":
            # 定义文件名
            name = os.path.basename(tempfile.NamedTemporaryFile().name)
            proj2000_img = os.path.join(os.path.dirname(filename), name + '_CGCS2000.tif')
            if os.path.exists(proj2000_img): os.remove(proj2000_img)
            # change proj
            GRID.rst_trans_proj(filename, proj2000_img, dst_ref="4490")
            filename = proj2000_img

        # get extent
        dataset = gdal.Open(filename)  # 打开文件:Tiff或者ENVI文件

        im_width = dataset.RasterXSize  # 栅格矩阵的列数
        im_height = dataset.RasterYSize  # 栅格矩阵的行数

        gtf = dataset.GetGeoTransform()  # 仿射矩阵
        left_lon = gtf[0]
        up_lat = gtf[3]
        right_lon = gtf[0] + (im_width + 0.05) * gtf[1] + (im_height + 0.05) * gtf[2]
        low_lat = gtf[3] + (im_width + 0.05) * gtf[4] + (im_height + 0.05) * gtf[5]

        del dataset
        if GRID.raster_proj(filename, "EPSG") != "4490": os.remove(proj2000_img)
        return left_lon, right_lon, low_lat, up_lat

    # # 像素坐标和地理坐标仿射变换
    @staticmethod
    def CoordTransf(Xpixel, Ypixel, GeoTransform):
        """
        像素坐标转为经纬度
        :param Xpixel: 像元位置x坐标
        :param Ypixel: 像元位置y坐标
        :param GeoTransform: 仿射矩阵
        :return: 经纬度
        """
        XGeo = GeoTransform[0] + GeoTransform[1] * Xpixel + Ypixel * GeoTransform[2]
        YGeo = GeoTransform[3] + GeoTransform[4] * Xpixel + Ypixel * GeoTransform[5]
        return XGeo, YGeo

    # 写文件，以写成tif为例
    @staticmethod
    def write_img(filename, im_proj, im_geotrans, im_data):
        """
        #gdal数据类型包括
        #gdal.GDT_Byte,
        #gdal .GDT_UInt16, gdal.GDT_Int16, gdal.GDT_UInt32, gdal.GDT_Int32,
        #gdal.GDT_Float32, gdal.GDT_Float64
        :param filename:
        :param im_proj:
        :param im_geotrans:
        :param im_data:
        :return:
        """
        # 判断栅格数据的数据类型
        # print('type---', im_data.dtype.name)
        if 'uint8' in im_data.dtype.name:
            datatype = gdal.GDT_Byte  # 范围：无符号整数（0 to 255）
        elif 'int16' in im_data.dtype.name:
            datatype = gdal.GDT_Int16  # 范围：整数（-32768 to 32767
        elif 'uint16' in im_data.dtype.name:
            datatype = gdal.GDT_UInt16  # 范围：无符号整数（0 to 65535）
        elif 'int32' in im_data.dtype.name:
            datatype = gdal.GDT_Int32  # 范围：整数（-2147483648 to 2147483647）
        elif 'uint32' in im_data.dtype.name:
            datatype = gdal.GDT_UInt32  # 范围：无符号整数（0 to 4294967295）
        elif 'float32' in im_data.dtype.name:
            datatype = gdal.GDT_Float32  # 范围：单精度浮点数，包括：1 个符号位，8 个指数位，23 个尾数位
        else:
            datatype = gdal.GDT_Float64  # 范围：双精度浮点数，包括：1 个符号位，11 个指数位，52 个尾数位

        # 判读数组维数
        # print(im_data.shape)
        if len(im_data.shape) == 3:
          im_bands, im_height, im_width = im_data.shape
        else:
          im_bands, (im_height, im_width) = 1, im_data.shape

        # 创建文件
        driver = gdal.GetDriverByName("GTiff")  # 数据类型必须有，因为要计算需要多大内存空间
        dataset = driver.Create(filename, im_width, im_height, im_bands, datatype)

        dataset.SetGeoTransform(im_geotrans)  # 写入仿射变换参数
        dataset.SetProjection(im_proj)  # 写入投影

        if im_bands == 1:
            dataset.GetRasterBand(1).WriteArray(im_data)  # 写入数组数据
        else:
            for i in range(im_bands):
                dataset.GetRasterBand(i + 1).WriteArray(im_data[i])

        del dataset

    #得到栅格的epsg
    @staticmethod
    def raster_proj(infile,out_flag=None):
        """

        :param infile:
        :param out_flag:
        :return:
        """
        d = gdal.Open(infile)
        proj = osr.SpatialReference(wkt=d.GetProjection())
        epsg = proj.GetAttrValue('AUTHORITY',1)
        if out_flag is None:
            print("请输入输出投影参数，'EPSG'或者'proj',类型是字符串")
        elif out_flag == r"EPSG":
            output =  epsg
        elif out_flag == r"proj":
            output =  proj

        return output

    # 栅格裁剪
    @staticmethod
    def raster_clip(infile,outfile,shpfile,dstNodata=None):
        """

        :param infile:
        :param outfile:
        :param shpfile:
        :param dstNodata:
        :return:
        """
        input_raster = gdal.Open(infile)
        if os.path.exists(outfile) : os.remove(outfile)
        if dstNodata is None:
            ds = gdal.Warp(outfile, input_raster, format='GTiff', cutlineDSName=shpfile,cropToCutline=True)
        else:
            ds = gdal.Warp(outfile, input_raster, format='GTiff', cutlineDSName=shpfile,# optionally you can filter your cutline (shapefile) based on attribute values
                         # cutlineWhere="FIELD = 'whatever'",# optionally you can filter your cutline (shapefile) based on attribute values
                         cropToCutline=True,dstNodata=np.nan)  # select the no data value you like
        ds = None
        return outfile

    # 投影转换
    @staticmethod
    def rst_trans_proj(infile,outfile,dst_ref=None):
        """
        给出矢量或者栅格文件或者epsg代码，并以此为标准进行输入栅格的投影转换
        infile:输入文件
        outfile:输出文件
        dstEPSG:目标epsg代码
        dst_ref:目标坐标系，矢量或者栅格文件或者epsg代码,字符串类型
        """
        if dst_ref.isdigit():
            dstEPSG = dst_ref
            # print(dstEPSG)
            # srs = ogr.osr.SpatialReference()
            # dstSRS = srs.ImportFromEPSG(4326)
            dstSRS =  "EPSG:" + dstEPSG
            # print(dstSRS)
        elif os.path.basename(dst_ref).split(".")[len(os.path.basename(dst_ref).split("."))-1] in [r"tif",r"dat",r"img"]:
            ds2 = gdal.Open(dst_ref)
            # 这里也可以直接用ds2.GetProjection()输出坐标系
            dstSRS = osr.SpatialReference(wkt=ds2.GetProjection())
            dstEPSG = dstSRS.GetAttrValue('AUTHORITY',1)
        elif os.path.basename(dst_ref).split(".")[len(os.path.basename(dst_ref).split("."))-1] in [r"shp"]:
            driver = ogr.GetDriverByName('ESRI Shapefile')
            dataset = driver.Open(dst_ref)
            # from Layer
            layer = dataset.GetLayer()
            spatialRef = layer.GetSpatialRef()
            # # from Geometry
            # feature = layer.GetNextFeature()
            # geom = feature.GetGeometryRef()
            # spatialRef = geom.GetSpatialReference()
            dstSRS = spatialRef
            dstEPSG = dstSRS.GetAttrValue('AUTHORITY', 1)
        else:
            print("出现错误，请检查输入参数！！！")
        ######
        src_epsg = GRID.raster_proj(infile,out_flag = 'EPSG')
        srcSRS = 'EPSG:' + src_epsg

        if int(dstEPSG) == int(src_epsg) :
            print("投影一致无须转换")
            outfile = infile
        else:
            options = gdal.WarpOptions(format='GTiff', srcSRS = srcSRS, dstSRS=dstSRS)
            gdal.Warp(outfile,infile,options = options)

        flag = True
        return flag,outfile


    # 对影像进行百分比拉伸
    @staticmethod
    def rst_stretch_percentclip(infile,outfile,percent=None, nodata=None,gamma = None):
        """

        :param infile:
        :param outfile:
        :param percent:
        :param nodata:
        :param gamma:
        :return:
        """
        im_proj, im_geotrans, im_data = GRID.read_img(infile)
        R_band = GRID.get_data(infile, 3)
        G_band = GRID.get_data(infile, 2)
        B_band = GRID.get_data(infile, 1)

        N_R = stretch_percentclip(R_band, percent=percent, nodata=nodata, gamma=gamma)
        N_G = stretch_percentclip(G_band, percent=percent, nodata=nodata, gamma=gamma)
        N_B = stretch_percentclip(B_band, percent=percent, nodata=nodata, gamma=gamma)
        new_data = np.array([N_R,N_G,N_B])
        GRID.write_img(outfile, im_proj, im_geotrans, new_data)

        return None