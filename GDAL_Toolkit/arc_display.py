# -*- coding: utf-8 -*-
# @author: ChuTianjia

import arcpy

arcpy.env.workspace = arcpy.GetParameterAsText(0)
mxd_file = arcpy.GetParameterAsText(1)
lyr_file = arcpy.GetParameterAsText(2)
mask_path = arcpy.GetParameterAsText(3)
new_lyr_path = arcpy.GetParameterAsText(4)
png_path = arcpy.GetParameterAsText(5)
dpi = arcpy.GetParameterAsText(6)

my_mxd = arcpy.mapping.MapDocument(mxd_file)
data_frame = arcpy.mapping.ListDataFrames(my_mxd)[0]
my_lyr = arcpy.mapping.Layer(lyr_file)
layer_list = arcpy.mapping.ListLayers(my_mxd)

my_mxd.activeView = "PAGE_LAYOUT"

tif_file_list = arcpy.ListRasters("BJ_hour_*", "TIF")

for raster in tif_file_list:
    # Import the mask layer into ArcMap
    raster_file = mask_path + "\\" + raster
    arcpy.MakeRasterLayer_management(raster_file, raster.strip(".tif"))

    # Modify the style of the mask layer according to the reference layer
    arcpy.ApplySymbologyFromLayer_management(raster.strip(".tif"), lyr_file)
    new_lyr_file = new_lyr_path + "\\" + raster.strip(".tif") + ".lyr"

    # Save and import the mask layer after modifying the style
    arcpy.SaveToLayerFile_management(raster.strip(".tif"), new_lyr_file)
    arcpy.AddMessage("{0} has been saved.".format(raster.strip(".tif") + ".lyr"))

    new_lyr = arcpy.mapping.Layer(new_lyr_file)
    arcpy.mapping.AddLayer(data_frame, new_lyr, "TOP")

    # Modify the image name
    for element in arcpy.mapping.ListLayoutElements(my_mxd, "TEXT_ELEMENT"):
        if element.name == "title":
            element.text = "Interpolation Map of PM2.5 Concentration\n at {0}:00 on May 18, 2019, Beijing".format(
                raster[8:10])

    new_lyr.visible = True

    # Modify the legend (see the program usage document for details)
    max_pixel = arcpy.GetRasterProperties_management(new_lyr, "MAXIMUM").getOutput(0)[0:5]
    min_pixel = arcpy.GetRasterProperties_management(new_lyr, "MINIMUM").getOutput(0)[0:5]
    for element in arcpy.mapping.ListLayoutElements(my_mxd, "TEXT_ELEMENT"):
        if element.name == "MAX":
            element.text = "{:0>5.2f}".format(float(max_pixel))
        if element.name == "MIN":
            element.text = "{:0>5.2f}".format(float(min_pixel))

    # Export to picture format
    png_file = png_path + "\\" + raster.strip(".tif") + ".png"
    arcpy.mapping.ExportToPNG(my_mxd, png_file, resolution=dpi)
    arcpy.AddMessage("{0} has been saved.".format(raster.strip(".tif") + ".png"))

    new_lyr.visible = False
    arcpy.mapping.RemoveLayer(data_frame, new_lyr[0])
