from osgeo import gdal, ogr, osr
import numpy as np
import cv2
import tempfile
import numpy.ma as ma
import os
import shapefile
import tqdm
import shutil
import geopandas as gpd
import glob
from image_base import GRID
import subprocess


# =====矢量处理基本功能函数（类）=====
# 矢量处理类

GRID = GRID()
class shpProcess:
    """
    @ Created by Sunhaoran
    """
    # 输出矢量的四至为矢量文件
    @staticmethod
    def shp_squareExtent(inShapefile,outShapefile):
        """
        根据四至，输出正方形范围
        :param inShapefile:
        :param outShapefile:
        :return:
        """
        # 结果存在就删掉
        basename = os.path.basename(outShapefile).split(".shp")[0]
        dirname = os.path.dirname(outShapefile)
        if os.path.exists(os.path.join(dirname, basename + ".shp")): os.remove(os.path.join(dirname, basename + ".shp"))
        if os.path.exists(os.path.join(dirname, basename + ".prj")): os.remove(os.path.join(dirname, basename + ".prj"))
        if os.path.exists(os.path.join(dirname, basename + ".dbf")): os.remove(os.path.join(dirname, basename + ".dbf"))
        if os.path.exists(os.path.join(dirname, basename + ".shx")): os.remove(os.path.join(dirname, basename + ".shx"))

        # Get a Layer's Extent
        inDriver = ogr.GetDriverByName("ESRI Shapefile")
        inDataSource = inDriver.Open(inShapefile, 0)
        inLayer = inDataSource.GetLayer()
        extent = inLayer.GetExtent()

        # print((extent[1] - extent[0]),(extent[3] - extent[2]))
        if (extent[1] - extent[0]) > (extent[3] - extent[2]):
            gap = (extent[1] - extent[0]) - (extent[3] - extent[2])
            extent = [extent[0], extent[1], extent[2] - gap/2, extent[3] + gap/2]
        else:
            gap = (extent[3] - extent[2]) - (extent[1] - extent[0])
            extent = [extent[0] - gap / 2, extent[1] + gap / 2, extent[2], extent[3]]

        # Create a Polygon from the extent tuple
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(extent[0], extent[2])
        ring.AddPoint(extent[1], extent[2])
        ring.AddPoint(extent[1], extent[3])
        ring.AddPoint(extent[0], extent[3])
        ring.AddPoint(extent[0], extent[2])
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)


        # Save extent to a new Shapefile
        outDriver = ogr.GetDriverByName("ESRI Shapefile")
        # Remove output shapefile if it already exists
        if os.path.exists(outShapefile):
          outDriver.DeleteDataSource(outShapefile)
        # Create the output shapefile
        outDataSource = outDriver.CreateDataSource(outShapefile)
        outLayer = outDataSource.CreateLayer("states_extent", geom_type=ogr.wkbPolygon)

        # Add an ID field
        idField = ogr.FieldDefn("id", ogr.OFTInteger)
        outLayer.CreateField(idField)

        # Create the feature and set values
        featureDefn = outLayer.GetLayerDefn()
        feature = ogr.Feature(featureDefn)
        feature.SetGeometry(poly)
        feature.SetField("id", 1)
        outLayer.CreateFeature(feature)
        # 定义投影
        # proj = osr.SpatialReference()
        # proj.ImportFromEPSG(4326)  # 4326-GCS_WGS_1984; 4490- GCS_China_Geodetic_Coordinate_System_2000
        ref_shp_ds = ogr.Open(inShapefile)
        ref_shp_layer = ref_shp_ds.GetLayer()
        ref_proj = ref_shp_layer.GetSpatialRef()
        wkt = ref_proj.ExportToWkt()
        # 写入投影
        f = open(outShapefile.replace(".shp", ".prj"), 'w')
        f.write(wkt)
        feature = None

        # Save and close DataSource
        inDataSource = None
        outDataSource = None

        flag = True
        return flag,extent

    # 根据四至输出矢量
    @staticmethod
    def extent_2shp(extent, outShapefile, epsg_val=None):
        # Create a Polygon from the extent tuple
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(extent[0], extent[2])
        ring.AddPoint(extent[1], extent[2])
        ring.AddPoint(extent[1], extent[3])
        ring.AddPoint(extent[0], extent[3])
        ring.AddPoint(extent[0], extent[2])
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)

        # Save extent to a new Shapefile
        outDriver = ogr.GetDriverByName("ESRI Shapefile")
        # Remove output shapefile if it already exists
        if os.path.exists(outShapefile):
            outDriver.DeleteDataSource(outShapefile)
        # Create the output shapefile
        outDataSource = outDriver.CreateDataSource(outShapefile)
        outLayer = outDataSource.CreateLayer("states_extent", geom_type=ogr.wkbPolygon)

        # Add an ID field
        idField = ogr.FieldDefn("id", ogr.OFTInteger)
        outLayer.CreateField(idField)

        # Create the feature and set values
        featureDefn = outLayer.GetLayerDefn()
        feature = ogr.Feature(featureDefn)
        feature.SetGeometry(poly)
        feature.SetField("id", 1)
        outLayer.CreateFeature(feature)
        # 定义投影
        # proj = osr.SpatialReference()
        # proj.ImportFromEPSG(4326)  # 4326-GCS_WGS_1984; 4490- GCS_China_Geodetic_Coordinate_System_2000
        # ref_shp_ds = ogr.Open(inShapefile)
        # ref_shp_layer = ref_shp_ds.GetLayer()
        # ref_proj = ref_shp_layer.GetSpatialRef()
        # wkt = ref_proj.ExportToWkt()

        # 4326-GCS_WGS_1984; 4490- GCS_China_Geodetic_Coordinate_System_2000
        if epsg_val is None: epsg_val = 4326
        proj = osr.SpatialReference()
        proj.ImportFromEPSG(epsg_val)
        wkt = proj.ExportToWkt()
        # 写入投影
        f = open(outShapefile.replace(".shp", ".prj"), 'w')
        f.write(wkt)

        return None

    # 输出矢量的四至为矢量文件，正常的四至范围
    @staticmethod
    def shp_extent(inShapefile, outShapefile):
        # 结果存在就删掉
        basename = os.path.basename(outShapefile).split(".shp")[0]
        dirname = os.path.dirname(outShapefile)
        if os.path.exists(os.path.join(dirname, basename + ".shp")): os.remove(os.path.join(dirname, basename + ".shp"))
        if os.path.exists(os.path.join(dirname, basename + ".prj")): os.remove(os.path.join(dirname, basename + ".prj"))
        if os.path.exists(os.path.join(dirname, basename + ".dbf")): os.remove(os.path.join(dirname, basename + ".dbf"))
        if os.path.exists(os.path.join(dirname, basename + ".shx")): os.remove(os.path.join(dirname, basename + ".shx"))

        # Get a Layer's Extent
        inDriver = ogr.GetDriverByName("ESRI Shapefile")
        inDataSource = inDriver.Open(inShapefile, 0)
        inLayer = inDataSource.GetLayer()
        extent = inLayer.GetExtent()

        # Create a Polygon from the extent tuple
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(extent[0], extent[2])
        ring.AddPoint(extent[1], extent[2])
        ring.AddPoint(extent[1], extent[3])
        ring.AddPoint(extent[0], extent[3])
        ring.AddPoint(extent[0], extent[2])
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)


        # Save extent to a new Shapefile
        outDriver = ogr.GetDriverByName("ESRI Shapefile")
        # Remove output shapefile if it already exists
        if os.path.exists(outShapefile):
            outDriver.DeleteDataSource(outShapefile)
        # Create the output shapefile
        outDataSource = outDriver.CreateDataSource(outShapefile)
        outLayer = outDataSource.CreateLayer("states_extent", geom_type=ogr.wkbPolygon)

        # Add an ID field
        idField = ogr.FieldDefn("id", ogr.OFTInteger)
        outLayer.CreateField(idField)

        # Create the feature and set values
        featureDefn = outLayer.GetLayerDefn()
        feature = ogr.Feature(featureDefn)
        feature.SetGeometry(poly)
        feature.SetField("id", 1)
        outLayer.CreateFeature(feature)
        # 定义投影
        # proj = osr.SpatialReference()
        # proj.ImportFromEPSG(4326)  # 4326-GCS_WGS_1984; 4490- GCS_China_Geodetic_Coordinate_System_2000
        ref_shp_ds = ogr.Open(inShapefile)
        ref_shp_layer = ref_shp_ds.GetLayer()
        ref_proj = ref_shp_layer.GetSpatialRef()
        wkt = ref_proj.ExportToWkt()
        # 写入投影
        f = open(outShapefile.replace(".shp", ".prj"), 'w')
        f.write(wkt)
        feature = None

        # Save and close DataSource
        inDataSource = None
        outDataSource = None

        flag = True
        return flag, extent

    # del shpfile
    @staticmethod
    def del_shpfile(file):
        dirname = os.path.dirname(file)
        basename = os.path.basename(file).split(".shp")[0]
        if os.path.exists(os.path.join(dirname, basename + ".shp")): os.remove(os.path.join(dirname, basename + ".shp"))
        if os.path.exists(os.path.join(dirname, basename + ".prj")): os.remove(os.path.join(dirname, basename + ".prj"))
        if os.path.exists(os.path.join(dirname, basename + ".dbf")): os.remove(os.path.join(dirname, basename + ".dbf"))
        if os.path.exists(os.path.join(dirname, basename + ".shx")): os.remove(os.path.join(dirname, basename + ".shx"))
        if os.path.exists(os.path.join(dirname, basename + ".cpg")): os.remove(os.path.join(dirname, basename + ".cpg"))
        return None

    @staticmethod
    def splitShape(inshp, split_field, outdir, encoding='gbk'):
        """
        按照某个字段把矢量分成一个一个
        :param inshp:
        :param split_field:
        :param outdir:
        :param encoding:
        :return:
        """
        sf = shapefile.Reader(inshp, encoding=encoding)
        fields = sf.fields
        records = sf.records()
        records_array = np.array(records)
        shapes = sf.shapes()
        shapes_array = np.array(shapes)
        field_values = []
        for rec in records:
            field_values.append(rec[split_field])
        field_values_array = np.array(field_values)
        unique_values = np.unique(field_values_array)
        if not os.path.isdir(outdir):
            os.makedirs(outdir, exist_ok=True)

        filelist = []
        for i in tqdm(range(len(unique_values))):
            temp_value = unique_values[i]
            locations = (field_values_array == temp_value)
            temp_records = (records_array[locations]).tolist()
            temp_shapes = (shapes_array[locations]).tolist()

            outshpfile = os.path.join(outdir, str(temp_value) + '.shp')
            writer = shapefile.Writer(outshpfile, encoding=encoding)
            writer.fields = fields
            for rec, shape in zip(temp_records, temp_shapes):
                writer.record(*rec)
                writer.shape(shape)
            writer.close()

            # Copy projection file
            in_prj_file = os.path.splitext(inshp)[0] + '.prj'
            o_prj_file = os.path.splitext(outshpfile)[0] + '.prj'
            shutil.copy(in_prj_file, o_prj_file)

            filelist.append(outshpfile)

        sf.close()
        return filelist
    # 矢量建立缓冲区，并输出缓冲区的矢量
    @staticmethod
    def createBuffer(inputfn,outputBufferfn,bufferDist):
        """
        bufferDist:单位和坐标系的单位一致
        """
        basename = os.path.basename(outputBufferfn).split(".shp")[0]
        dirname = os.path.dirname(outputBufferfn)
        if os.path.exists(os.path.join(dirname, basename + ".shp")): os.remove(os.path.join(dirname, basename + ".shp"))
        if os.path.exists(os.path.join(dirname, basename + ".prj")): os.remove(os.path.join(dirname, basename + ".prj"))
        if os.path.exists(os.path.join(dirname, basename + ".dbf")): os.remove(os.path.join(dirname, basename + ".dbf"))
        if os.path.exists(os.path.join(dirname, basename + ".shx")): os.remove(os.path.join(dirname, basename + ".shx"))

        inputds = ogr.Open(inputfn)
        inputlyr = inputds.GetLayer()

        shpdriver = ogr.GetDriverByName('ESRI Shapefile')
        if os.path.exists(outputBufferfn):
            shpdriver.DeleteDataSource(outputBufferfn)
        outputBufferds = shpdriver.CreateDataSource(outputBufferfn)
        bufferlyr = outputBufferds.CreateLayer(outputBufferfn, geom_type=ogr.wkbPolygon)
        featureDefn = bufferlyr.GetLayerDefn()

        for feature in inputlyr:
            ingeom = feature.GetGeometryRef()
            geomBuffer = ingeom.Buffer(bufferDist)

            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(geomBuffer)
            bufferlyr.CreateFeature(outFeature)
            outFeature = None

        ref_shp_ds = ogr.Open(inputfn)
        ref_shp_layer = ref_shp_ds.GetLayer()
        ref_proj = ref_shp_layer.GetSpatialRef()
        wkt = ref_proj.ExportToWkt()
        # 写入投影
        f = open(outputBufferfn.replace(".shp", ".prj"), 'w')
        f.write(wkt)
        flag = True
        return outputBufferfn

    # 获得矢量的投影坐标系
    def shp_proj(self,infile,out_flag=None):
        ds = ogr.Open(infile)
        driver = ds.GetDriver()
        driver_type = driver.GetName()

        # 打开一个图层
        layer_nums = ds.GetLayerCount()
        for i in range(layer_nums):
            layer = ds.GetLayerByIndex(i)
            layerDefn = layer.GetLayerDefn()
            layer_name = layerDefn.GetName()

        # 获取epsg
        crs = layer.GetSpatialRef()
        epsg = crs.GetAttrValue('AUTHORITY', 1)

        if out_flag is None:
            print("请输入输出投影参数，'epsg'或者'proj',类型是字符串")
        elif out_flag == r"epsg":
            output = epsg
        elif out_flag == r"proj":
            output = crs

        return output

    # 循环逐个输出矢量
    @staticmethod
    def loop_outsubshp(inshp,outpath):
        """
        循环逐个输出矢量的子集，文件名加序号，但是目前该算法未实现无法保留原有字段
        :param inshp:
        :param outpath:
        :return: 文件列表
        """
        inDriver = ogr.GetDriverByName("ESRI Shapefile")
        inDataSource = inDriver.Open(inshp, 0)
        inputlyr = inDataSource.GetLayer()

        filelist = []
        for fi in range(len(inputlyr)):
            # define out file name
            basename = os.path.basename(inshp).split(".shp")[0]
            outname = basename + "_" + str(fi) + ".shp"
            outfile = os.path.join(outpath, outname)

            # define output
            shpdriver = ogr.GetDriverByName('ESRI Shapefile')
            if os.path.exists(outfile):
                shpdriver.DeleteDataSource(outfile)
            outputBufferds = shpdriver.CreateDataSource(outfile)
            # outlyr = outputBufferds.CreateLayer(outfile, geom_type=ogr.wkbPolygon)
            # 为图层指定编码
            outlyr = outputBufferds.CreateLayer("CommercialHousing", srs=inputlyr.GetSpatialRef(),
                                                geom_type=ogr.wkbPolygon, options=["ENCODING=gbk2312"])
            # Get the output Layer's Feature Definition
            featureDefn = outlyr.GetLayerDefn()

            # Add input Layer Fields to the output Layer if it is the one we want
            inLayerDefn = inputlyr.GetLayerDefn()
            for i in range(0, inLayerDefn.GetFieldCount()):
                fieldDefn = inLayerDefn.GetFieldDefn(i)
                fieldName = fieldDefn.GetName()
                outlyr.CreateField(fieldDefn)

            # get output file
            feature = inputlyr.GetFeature(fi)
            ingeom = feature.GetGeometryRef()
            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(ingeom)
            outlyr.CreateFeature(outFeature)

            # Add field values from input Layer
            for i in range(0, featureDefn.GetFieldCount()):
                fieldDefn = featureDefn.GetFieldDefn(i)
                fieldName = fieldDefn.GetName()
                outFeature.SetField(featureDefn.GetFieldDefn(i).GetNameRef(), feature.GetField(i))

            outFeature = None
            #
            ref_shp_ds = ogr.Open(inshp)
            ref_shp_layer = ref_shp_ds.GetLayer()
            ref_proj = ref_shp_layer.GetSpatialRef()
            wkt = ref_proj.ExportToWkt()
            # 写入投影
            f = open(outfile.replace(".shp", ".prj"), 'w')
            f.write(wkt)

            filelist.append(outfile)

        del inDataSource, outputBufferds

        return filelist

    # 矢量求交集
    @staticmethod
    def shps_get_intersection(gdf1_file, gdf2_file, outpath=None):
        """
        gdf1_file:
        gdf2_file:
        outpath:输出路径,未定义默认输出到gdf1_file同文件下
        返回两者交集
        """
        # 得到输出文件名
        basename2 = os.path.basename(gdf2_file).split(".shp")[0]
        basename1 = os.path.basename(gdf1_file).split(".shp")[0]
        if outpath is None:
            dirname = os.path.dirname(gdf1_file)
            gdf1_subinter_file = os.path.join(dirname,basename1 + r"_subinter.shp")
            gdf_inter_file = os.path.join(dirname,basename1 + r"_" + basename2 + r"_inter.shp")
            gdf2_subinter_file = os.path.join(dirname,basename2 + r"_subinter.shp")
            gdf1_gdf2_unionfile = os.path.join(dirname,basename1 + r"_" + basename2 + r"_union.shp")
        else:
            gdf1_subinter_file = os.path.join(outpath, basename1 + r"_subinter.shp")
            gdf_inter_file = os.path.join(outpath, basename1 + r"_" + basename2 + r"_inter.shp")
            gdf2_subinter_file = os.path.join(outpath, basename2 + r"_subinter.shp")
            gdf1_gdf2_unionfile = os.path.join(outpath,basename1 + r"_" + basename2 + r"_union.shp")

        # del files
        shpProcess.del_shpfile(gdf1_subinter_file)
        shpProcess.del_shpfile(gdf_inter_file)
        shpProcess.del_shpfile(gdf2_subinter_file)
        shpProcess.del_shpfile(gdf1_gdf2_unionfile)

        # 求出交集部分
        gdf1 = gpd.read_file(gdf1_file)
        gdf2 = gpd.read_file(gdf2_file).to_crs(gdf1.crs)
        gdf_inter = gpd.overlay(gdf1, gdf2, how='intersection')
        gdf_inter.to_file(gdf_inter_file, encoding='gb2312')

        return gdf_inter_file

    # 矢量叠加、相减运算
    @staticmethod
    def shps_spatial_analysis(gdf1_file, gdf2_file, outpath=None):
        """
        gdf1_file:
        gdf2_file:
        outpath:输出路径,未定义默认输出到gdf1_file同文件下
        返回第一个剩余的的部分、两者交集、第二个剩余的部分
        """
        # 得到输出文件名
        basename2 = os.path.basename(gdf2_file).split(".shp")[0]
        basename1 = os.path.basename(gdf1_file).split(".shp")[0]
        if outpath is None:
            dirname = os.path.dirname(gdf1_file)
            gdf1_subinter_file = os.path.join(dirname,basename1 + r"_subinter.shp")
            gdf_inter_file = os.path.join(dirname,basename1 + r"_" + basename2 + r"_inter.shp")
            gdf2_subinter_file = os.path.join(dirname,basename2 + r"_subinter.shp")
            gdf1_gdf2_unionfile = os.path.join(dirname,basename1 + r"_" + basename2 + r"_union.shp")
        else:
            gdf1_subinter_file = os.path.join(outpath, basename1 + r"_subinter.shp")
            gdf_inter_file = os.path.join(outpath, basename1 + r"_" + basename2 + r"_inter.shp")
            gdf2_subinter_file = os.path.join(outpath, basename2 + r"_subinter.shp")
            gdf1_gdf2_unionfile = os.path.join(outpath,basename1 + r"_" + basename2 + r"_union.shp")

        # del files
        shpProcess.del_shpfile(gdf1_subinter_file)
        shpProcess.del_shpfile(gdf_inter_file)
        shpProcess.del_shpfile(gdf2_subinter_file)
        shpProcess.del_shpfile(gdf1_gdf2_unionfile)

        # 求出交集部分
        gdf1 = gpd.read_file(gdf1_file)
        gdf2 = gpd.read_file(gdf2_file).to_crs(gdf1.crs)
        gdf_inter = gpd.overlay(gdf1, gdf2, how='intersection')
        gdf_inter.to_file(gdf_inter_file, encoding='gb2312')

        # 求出gdf1扣出交集的剩余部分
        gdf_diff1 = gpd.overlay(gdf1, gdf2, how='difference', make_valid=True)
        gdf_diff1.to_file(gdf1_subinter_file, encoding='gb2312')

        # 求出gdf2扣出交集的剩余部分
        gdf_diff2 = gpd.overlay(gdf2, gdf1, how='difference', make_valid=True)
        gdf_diff2.to_file(gdf2_subinter_file, encoding='gb2312')

        # 求出并集
        gdf1_gdf2_union = gpd.overlay(gdf1, gdf2, how='union', make_valid=True)
        gdf1_gdf2_union.to_file(gdf1_gdf2_unionfile, encoding='gb2312')

        # shpfiles
        shpfiles = [gdf1_subinter_file,gdf_inter_file,gdf2_subinter_file,gdf1_gdf2_unionfile]
        return shpfiles

    # 重投影函数，要求输出的四个主要文件都不存在
    @staticmethod
    def shp_trans_proj(input_shp,out_shp,ref_shp=None, manual_srs=None):
        # if manual_srs == None:
        #   ref_shp_ds = ogr.Open(ref_shp)
        #   ref_shp_layer = ref_shp_ds.GetLayer()
        #   ref_proj = ref_shp_layer.GetSpatialRef()
        # else:
        #   ref_proj = manual_srs
        if  ((manual_srs == True) or(ref_shp == True)):
            if manual_srs == None:
                ref_shp_ds = ogr.Open(ref_shp)
                ref_shp_layer = ref_shp_ds.GetLayer()
                ref_proj = ref_shp_layer.GetSpatialRef()
            else:
                ref_proj = manual_srs
        else:
            ref_proj = osr.SpatialReference()
            ref_proj.ImportFromEPSG(4326)

        basename = os.path.basename(out_shp).split(".tif")[0]
        dirname = os.path.dirname(out_shp)
        if os.path.exists(os.path.join(dirname, basename + ".shp")): os.remove(os.path.join(dirname, basename + ".shp"))
        if os.path.exists(os.path.join(dirname, basename + ".prj")): os.remove(os.path.join(dirname, basename + ".prj"))
        if os.path.exists(os.path.join(dirname, basename + ".dbf")): os.remove(os.path.join(dirname, basename + ".dbf"))
        if os.path.exists(os.path.join(dirname, basename + ".shx")): os.remove(os.path.join(dirname, basename + ".shx"))

        cmd = ['ogr2ogr', '-t_srs', str(ref_proj), out_shp, input_shp]
        if not os.path.isfile(out_shp):
            DEVNULL = open(os.devnull, 'wb')
            subprocess.call(cmd, stdout=DEVNULL, stderr=subprocess.STDOUT)

    # 根据矢量转化为同投影同分辨率影像
    @staticmethod
    def shp_convert_raster(shpfile, outpath, cellsize, field=None, burn_val=None, backval=None):
        """
        根据矢量转化为同投影同分辨率影像；
        非常要注意：当指定字段时，val =None,backval=None无效
        :param shpfile: 字符串类型，需要被转化的矢量
        :param cellsize: 数值型，像元大小;如果矢量是经纬度投影，那么值可以定义为0.00016578866度；如果矢量是UTM、兰伯特、高斯等投影，则定义为30米
        :param outpath: 字符串类型，输出路径
        :param field:字符串类型，需要转化的字段
        :param burn_val:int，输入255，1，或其他值
        :param backval:int，背景值，0，255，或其他值
        :return:字符串类型，文件
        author： created by xuxuan,2022-06-25
        """
        # define outfile
        shp_basename = os.path.basename(shpfile).split(".shp")[0]
        outname = shp_basename + "_" + "ShpToRaster" + ".tif"
        outfile = os.path.join(outpath, outname)
        # 如文件存在，则删除，否则影响输出
        if os.path.exists(outfile):
            basename = os.path.basename(outfile).split(".")[0]
            find_file = os.path.join(outpath, basename + "*")
            files = glob.glob(find_file)
            for file in files:
                os.remove(file)

        # read shp
        shp = ogr.Open(shpfile)
        lyr = shp.GetLayer()
        srs = lyr.GetSpatialRef()
        # Extent
        x_min, x_max, y_min, y_max = lyr.GetExtent()
        x_ncells = int((x_max - x_min) / cellsize)
        y_ncells = int((y_max - y_min) / cellsize)

        # 判断输出栅格类型，没有就输出0，1；有例如ndvi就根据某个值同化
        gdf = gpd.read_file(shpfile)
        # 数据类型对应关系参考：https://blog.csdn.net/X_Cosmic/article/details/107182061
        if field is None:
            output_datatype = gdal.GDT_Byte
        elif gdf[field].dtypes == 'int8':
            output_datatype = gdal.GDT_Int16
        elif gdf[field].dtypes == 'int16':
            output_datatype = gdal.GDT_Int16
        elif gdf[field].dtypes == 'int32':
            output_datatype = gdal.GDT_Int32
        elif gdf[field].dtypes == 'int64':
            output_datatype = gdal.GDT_Int32
        elif gdf[field].dtypes == 'uint8':
            output_datatype = gdal.GDT_Byte  # 0-255
        elif gdf[field].dtypes == 'uint16':
            output_datatype = gdal.GDT_UInt16  # 0-65535
        elif gdf[field].dtypes == 'uint32':
            output_datatype = gdal.GDT_UInt32  # 0-4294967295
        elif gdf[field].dtypes == 'uint64':
            output_datatype = gdal.GDT_UInt32
        else:
            output_datatype = gdal.GDT_Float32

        # create output raster
        # Output
        out_driver = gdal.GetDriverByName('GTiff')
        if os.path.exists(outfile):
            out_driver.Delete(outfile)
        dst_ds = out_driver.Create(outfile, x_ncells, y_ncells, 1, output_datatype)
        dst_ds.SetGeoTransform((x_min, cellsize, 0, y_max, 0, -cellsize))
        dst_ds.SetProjection(srs.ExportToWkt())

        # set burnvalue，gdal默认是0
        if backval is not None:
            band = dst_ds.GetRasterBand(1)
            band.SetNoDataValue(backval)

        if field is None:
            if burn_val is None: burn_val = 1
            gdal.RasterizeLayer(dst_ds, [1], lyr, burn_values=[burn_val])
        else:
            OPTIONS = ['ATTRIBUTE=' + field]
            gdal.RasterizeLayer(dst_ds, [1], lyr, None, options=OPTIONS)
        dst_ds, shp, lyr = None, None, None
        return outfile

    # 求多边形内最大面积矩形的中心点
    @staticmethod
    def shp_get_MaxRectCenter(shpfile,coff,coor_class):
        """
        根据输入的矢量，求多边形内最大矩形的中心
        注意：该算法主要用于求最适宜放置标注的位置
        :param shp:str, input shp file absolute path
        :param coff:int, 系数，县级：30;市级：300；省级：3000
        :param coor_class:str, "geo" or "proj"
        :return:[float,float]数组: [long,lat], return center point's coordinate。
        """
        # 读取矢量并转化为同投影或坐标系栅格
        cashedir = mkdtemp()
        # resIMGfile = shpProcess.shp_convert_raster(shpfile, cashedir, cellsize = 0.00016578866, field=None, burn_val=None, backval=None)
        # cellsize=0.00016578866,18m
        if coor_class == "geo":
            resIMGfile = shpProcess.shp_convert_raster(shpfile, cashedir, cellsize=0.00016578866 * 2 * coff, field=None, burn_val=None,backval=None)
        elif coor_class == "proj":
            resIMGfile = shpProcess.shp_convert_raster(shpfile, cashedir, cellsize=30 * coff, field=None,burn_val=None, backval=None)
        else:
            print("在函数shp_get_MaxRectCenter中，coor_class 参数错误")
        # 利用栅格求取中心点
        lon, lat = GRID.img_get_BestLabelLocation(resIMGfile)
        shutil.rmtree(cashedir)
        return lon,lat

    @staticmethod
    def shp_outsubshp_byattri(shp, outpath, field, field_val):
        """
        把满足条件的属性，保存为新的矢量
        :param shp: 字符串型，输入文件
        :param outpath: 输出路径，需要保存的目录
        :param field: 字符串，字段名
        :param field_val: 字符串或者数值型，根据字段确定
        :return:
        """
        # 读取矢量到地理数据框geopandasfram
        gdf = gpd.read_file(shp)

        # 过滤条件
        gdf_filtered = gdf[(gdf[field] == field_val)]

        # 定义文件名并输出文件
        basename = os.path.basename(shp).split(".shp")[0]
        outname = basename + "_" + str(field_val) + ".shp"
        outfile = os.path.join(outpath,outname)
        gdf_filtered.to_file(outfile, "ESRI Shapefile", encoding='gb2312')
        return outfile
    @staticmethod
    def shp_outsubshpslist_byattri(shpfile, outpath, field):
        """
        根据属性，遍历保存为新的矢量
        :param shpfile: 字符串型，输入文件
        :param outpath: 输出路径，需要保存的目录
        :param field: 字符串，字段名
        :return:
        """
        gdf = gpd.read_file(shpfile)
        code_val_list = []
        # 求不同属性组合列表
        for ai in range(len(gdf[field])):
            code_str = gdf[field][ai]
            if code_str  not in code_val_list:
                code_val_list.append(code_str)
        # 循环输出矢量
        shps_list = []
        for ai in range(len(code_val_list)):
            outfile = shpProcess.shp_outsubshp_byattri(shpfile, outpath, "code", code_val_list[ai])
            shps_list.append(outfile)
        return shps_list