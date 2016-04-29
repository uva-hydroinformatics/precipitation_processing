import arcpy
from arcpy import env

# get data from csv
table_path = 'C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/data/2014.09.08.xls'

# set up arcpy environment
env.overwriteOutput = True
env.workspace = "C:/Users/jeff_dsktp/Documents/ArcGIS"

# convert to table then to layer then to shapefile... phew!
date = "2014_09_08"
tbl = "{}.dbf".format(date)
arcpy.ExcelToTable_conversion(table_path, tbl)

spref = r"Coordinate Systems/Projected Coordinate Systems/State Plane/Nad 1983/NAD 1983 HARN StatePlane Virginia South FIPS 4502 (Meters).prj"
lyr = "{}_Layer".format(date)
arcpy.MakeXYEventLayer_management(tbl, 'X','Y', lyr, spref)

shp_path = "C:/Users/jeff_dsktp/Google Drive/Hampton Roads GIS Data/VA_Beach_Data/rain_gauges/"
arcpy.CopyFeatures_management(lyr, '{}.shp'.format(date.replace('_',"")))
arcpy.Delete_management(tbl)

# model params    ... todo: get these from R output
model_type = "Spherical"
lag_size = "1000.000000"
range = "10529.000000"
sill = "1300"
nugget = "0"

# do the kriging
out_dir = "C:/Users/jeff_dsktp/Google Drive/Hampton Roads GIS Data/VA_Beach_Data/kriging/"
out_est_file = "09082014_5"
cell_size = "61.3231323600002"
sample_type = "Variable"
sample_num = "27"
out_var_file = "sv09082014_5"
arcpy.gp.Kriging_sa("{}.shp".format(date.replace('_',"")), "f{}".format(date[:-1]), "{}{}".format(out_dir, out_est_file), "{} {} {} {} {}".format(model_type, lag_size,range, sill, nugget), cell_size, "{} {}".format(sample_type, sample_num), "{}{}".format(out_dir, out_var_file))


# Todo: make the model parameters dynamic: have R script output csv and read them in here
# Todo: put in loop for each of the time steps
# Todo: loop through each station removing one by one and see the effect on the rainfall estimation