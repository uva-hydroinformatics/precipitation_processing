import arcpy
from arcpy import env
env.workspace = "C:/Users/jeff_dsktp/Documents/ArcGIS"
table_path = 'C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/data/2014.09.08.xls'
date = "2014_09_08"
tbl = "{}.dbf".format(date)
arcpy.ExcelToTable_conversion(table_path, tbl)
spref = r"Coordinate Systems\Projected Coordinate Systems\Utm\Nad 1983\NAD 1983 UTM Zone 11N.prj"
lyr = "{}_Layer".format(date)
arcpy.MakeXYEventLayer_management(tbl, 'X','Y', lyr)
shp_path = 'C:/Users/jeff_dsktp/Google Drive/Hampton Roads GIS Data/VA_Beach_Data/rain_gauges/'
arcpy.CopyFeatures_management(lyr, '{}.shp'.format(date))
arcpy.Delete_management(tbl)

# model params
model_type = "Spherical"
lag_size = "1000.000000"
range = "10529.000000"
sill = "1300"
nugget = "0"

out_dir = "C:/Users/jeff_dsktp/Google Drive/Hampton Roads GIS Data/VA_Beach_Data/kriging/"
out_est_file = "09082014_3"
cell_size = "61.3231323600002"
sample_type = "Variable"
sample_num = "27"
out_var_file = "sv09082014_3"
arcpy.gp.Kriging_sa(date, date, "{}{}".format(out_dir, out_est_file), "{} {} {} {} {}".format(model_type, lag_size,range, sill, nugget), cell_size, "{} {}".format(sample_type, sample_num), "{}{}".format(out_dir, out_var_file))


# Todo: add spatial reference to data