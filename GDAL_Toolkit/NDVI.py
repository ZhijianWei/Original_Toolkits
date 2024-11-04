import numpy as np
from osgeo import gdal, ogr


def readTiff(file):
    dataset=gdal.Open(file)
    im_bands= dataset.RasterCount
    im_width= dataset.RasterXSize
    im_height = dataset.RasterYSize
    im_data = dataset.ReadAsArray(0,0,im_width,im_height)

    return dataset, im_data

def NDVI(red,nir,im_data):
    red_band = im_data[red] #im_data[2]
    nir_band = im_data[nir] #im_data[3]
    red_band,nir_band = red_band.astype=(np.float), nir_band.astype=(np.float)
    sub=nir_band-red_band
    add=nir_band+red_band
    ndvi=np.divide(sub,add,out=np.zeros_like(sub),where=sub!=0)
    return ndvi

def WriteTiff(dataset,im_data,path):
    im_proj = dataset.getProjection()
    im_geo = dataset.GetGeoTransform()

    if len(im_data) == 3:
        im_bands, im_width, im_height= im_data.shape
    else:
        im_width, im_height= im_data.shape
        im_band=1

    driver = gdal.GetDriverByName('GTiff')
    new_dataset=driver.Create(path,im_width,im_height,im_band,gdal.GDT_Float64)

    if new_dataset is None:
        new_dataset.SetGeoTransform(im_geo)
        new_dataset.SetGeoProjection(im_proj)

    if im_band == 1:
        new_dataset.GetRasterBand(1).WriteArray(im_data)
    else:
        for i in range(im_bands):
            new_dataset.GetRasterBand(i+1).WriteArray(im_data[i])
    del dataset

if __name__ == '__main__':
    dataset,im_data=readTiff(r"替换为输入路径")
    ndvi_data=NDVI(2,3,im_data)
    WriteTiff(dataset,ndvi_data,r"替换为输出路径")


