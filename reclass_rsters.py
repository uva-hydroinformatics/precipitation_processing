# Name: reclassify_example02.py
# Description: Reclassifies the values in a raster.
# Requirements: Spatial Analyst Extension

# Import system modules
import arcpy
from arcpy import env
from arcpy.sa import *
from storm_stats_functions import data_dir, check_dir
import os
import numpy as np
from datetime import datetime


def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx


def closest_to_hour(f_names):
    no_extension = [f.replace('.tif', '') for f in f_names]
    f_hours = np.array([f.split('_')[-1] for f in no_extension])
    f_hours = f_hours.astype(int)
    check_hours = np.arange(0, 240000, step=10000)
    return [f_names[find_nearest(f_hours, hr)] for hr in check_hours]


# Set environment settings
rain_date = "20151002"
data_dir = "{}/nexrad/{}".format(data_dir, rain_date)
check_dir("{}/sum".format(data_dir))
check_dir("{}/reclass".format(data_dir))
env.workspace = data_dir
# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")
env.overwriteOutput = True

filenames = [f for f in os.listdir(data_dir) if f.endswith('.tif')]
hours = closest_to_hour(filenames)
i = 0
for fname in hours:
    rast_arr = arcpy.RasterToNumPyArray(fname)
    # check if date time is between boundaries
    if np.isfinite(rast_arr).sum() > 0:
        print "calculating for {}".format(fname)
        # Set local variables
        inRaster = Raster(fname)
        outReclassRR = Con(inRaster < 0, 0, inRaster)
        if i == 0:
            sumRaster = outReclassRR
        else:
            sumRaster += outReclassRR
        i += 1

            # Save the output
        outReclassRR.save("{}/reclass/rcls_{}".format(data_dir, fname))
    else:
        print "{} file seems to be all null".format(fname)
        continue
sumRaster *= 25.4
arcpy.Clip_management(sumRaster, rectangle="-76.232 36.735 -75.941 36.937",
                      out_raster="{}/sum/sum_{}.tif".format(data_dir, rain_date))
