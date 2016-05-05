import arcpy
from arcpy import env
from arcpy.sa import *
import pandas as pd
from storm_stats_functions import check_dir

# get data from csv
data_dir = 'C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/data/'
table_name = '2014.09.08'
ext = '.xls'
table_path = "{}{}{}".format(data_dir, table_name, ext)

# set up arcpy environment
k_dir = check_dir('C:/Users/jeff_dsktp/Google Drive/Hampton Roads GIS Data/VA_Beach_Data/kriging')
env.overwriteOutput = True

env.workspace = k_dir

# convert to table then to layer then to shapefile... phew!
date = "2014_09_08"
tbl = "{}.dbf".format(date)
arcpy.ExcelToTable_conversion(table_path, tbl)

spref = r"Coordinate Systems/Projected Coordinate Systems/State Plane/Nad 1983/NAD 1983 HARN StatePlane Virginia South FIPS 4502 (Meters).prj"
lyr = "{}_Layer".format(date)
arcpy.MakeXYEventLayer_management(tbl, 'X', 'Y', lyr, spref)

shp_path = check_dir("C:/Users/jeff_dsktp/Google Drive/Hampton Roads GIS Data/VA_Beach_Data/rain_gauges/{}".format(date))
arcpy.CopyFeatures_management(lyr, '{}.shp'.format(date.replace('_',"")))

# model params    ... todo: get these from R output
param_file = "{}{}_model_params.csv".format(data_dir, 'daily_tots')
param_df = pd.read_csv(param_file)
param_df.set_index('date', inplace=True)
d = date.replace("_",".")
model_type = param_df['type'][d]
lag_size = "1000.000000"
range = param_df['range'][d]
sill = param_df['sill'][d]
nugget = "0"
kmodel = KrigingModelOrdinary(model_type, lag_size, range, sill, nugget)
sample_num = "27"
krad = arcpy.sa.RadiusVariable(27, 20000)

# do the kriging
out_est_file = date.replace("_","")
cell_size = "61.3231323600002"
sample_type = "Variable"
out_var_file = "sv09082014"
arcpy.gp.Kriging_sa("{}.shp".format(date.replace('_',"")),
                   "f{}".format(date[:-1]),
                   out_est_file,
                   "{} {} {} {} {}".format(model_type, lag_size,range, sill, nugget),
                   cell_size,
                   "{} {}".format(sample_type, sample_num),
                   out_var_file)


# Todo: make the model parameters dynamic: have R script output csv and read them in here
# Todo: put in loop for each of the time steps
# Todo: loop through each station removing one by one and see the effect on the rainfall estimation



# todo: read in watershed shapefile then select each watershed one by one and calculate the zonal stats
shed_ply = "C:/Users/jeff_dsktp/Google Drive/Hampton Roads GIS Data/VA_Beach_Data/Problem Spots/problem_watersheds.shp"
arcpy.MakeTableView_management(shed_ply, "table_view")
num_rows = int(arcpy.GetCount_management("table_view").getOutput(0))
est_rast = "C:/Users/jeff_dsktp/Google Drive/Hampton Roads GIS Data/VA_Beach_Data/kriging/09082014_5"
out_tab = "temp_tab"
for i in range(1):
    i = 1
    arcpy.SelectLayerByAttribute_management("problem_watersheds", "NEW_SELECTION", " FID = {}".format(i+1))
    arcpy.sa.ZonalStatisticsAsTable("problem_watersheds", "FID", est_rast, out_tab, "DATA", "MEAN")
    cur = arcpy.da.SearchCursor(out_tab)
    for row in cur:
        print row[-1]


# Todo: make the model parameters dynamic: have R script output csv and read them in here
# Todo: put in loop for each of the time steps
# Todo: loop through each station removing one by one and see the effect on the rainfall estimation