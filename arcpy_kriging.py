import arcpy
from arcpy import env
from arcpy.sa import *
import pandas as pd
from storm_stats_functions import check_dir
import numpy as np


class ModelParams():
    def __init__(self, tpe, d_dir, d):
        param_file = "{}{}_model_params.csv".format(d_dir, tpe)
        param_df = pd.read_csv(param_file)
        param_df.set_index('date', inplace=True)
        d = d.replace("-", ".")

        self.model_type = param_df['type'][d]
        self.lag_size = "1000.000000"
        self.range = param_df['range'][d]
        self.sill = param_df['sill'][d]
        self.nugget = "0"                                                                 # todo get this from R script
        self.sample_num = "27"  # todo make this dynamic?
        self.krad = arcpy.sa.RadiusVariable(27, 20000)
        self.sample_type = "Variable"


#specify type
tpe = 'daily_tots'

# get data from table
data_dir = 'C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/data/'
table_name = 'dtots'
ext = '.xls'
table_path = "{}{}{}".format(data_dir, table_name, ext)
df = pd.read_excel(table_path)
dates = df.columns[4:]


# set up arcpy environment
arcpy.CheckOutExtension("spatial")
env.extent = arcpy.Extent(3705690, 1051630, 3724920, 1068584)
k_dir = check_dir('C:/Users/jeff_dsktp/Google Drive/Hampton Roads GIS Data/VA_Beach_Data/kriging')
env.overwriteOutput = True
env.workspace = k_dir

# get watershed polygon and iterate through them
shed_ply = "C:/Users/jeff_dsktp/Google Drive/Hampton Roads GIS Data/VA_Beach_Data/Problem Spots/problem_watersheds.shp"
arcpy.MakeTableView_management(shed_ply, "table_view")
num_rows = int(arcpy.GetCount_management("table_view").getOutput(0))
out_tab = "temp_tab"
means = []
for i in np.arange(1):
    # select individual watershed
    sel = "selection.shp"
    arcpy.Select_analysis(shed_ply, sel, '"FID" = {}'.format(i))
    cur = arcpy.SearchCursor(sel)
    for row in cur:
        wshed_descr = row.getValue("Descript")
        wshed_area = row.getValue("Area_sq_km")
        wshed_x = row.getValue("x")
        wshed_y = row.getValue("y")
    if wshed_area < 0.001:
        continue

    # take out nearest points and resave as same name
    j = 0

    # convert to table then to layer then to shapefile... phew!
    for date in [dates[0]]:
        date_df = df[['site_name', 'x', 'y', 'src', date]]  # get df for individual timestep
        date_df = date_df[date_df[date].notnull()]
        date_table = "{}{}.xls".format(data_dir, date)
        date_df.to_excel(date_table, index=False)

        tbl = "tab.dbf".format(date.replace("-", ""))

        arcpy.ExcelToTable_conversion(date_table, tbl, Sheet="Sheet1")

        spref = r"Coordinate Systems/Projected Coordinate Systems/State Plane/Nad 1983/" \
                r"NAD 1983 HARN StatePlane Virginia South FIPS 4502 (Meters).prj"
        lyr = "Layer".format(date)
        arcpy.MakeXYEventLayer_management(tbl, 'x', 'y', lyr, spref)

        arcpy.CopyFeatures_management(lyr, '{}.shp'.format(date.replace('-', "")))

        #  get model params
        mp = ModelParams(tpe, data_dir, date)

        # do kriging
        out_est_file = date.replace("-", "")
        cell_size = "61.3231323600002"
        out_var_file = "sv{}".format(date.replace("-", ""))
        arcpy.gp.Kriging_sa("{}.shp".format(date.replace('-', "")),
                            "f{}".format(date[:-1]),
                            out_est_file,
                            "{} {} {} {} {}".format(mp.model_type, mp.lag_size, mp.range, mp.sill, mp.nugget),
                            cell_size,
                            "{} {}".format(mp.sample_type, mp.sample_num),
                            out_var_file)

        # get mean of the kriged surface of selected watershed
        arcpy.sa.ZonalStatisticsAsTable(sel, "FID", out_est_file, out_tab, "DATA", "MEAN")
        cur = arcpy.da.SearchCursor(out_tab, "MEAN")
        for row in cur:
            means.append(row[0])


# Todo: put in loop for each of the time steps
# Todo: loop through each station removing one by one and see the effect on the rainfall estimation


