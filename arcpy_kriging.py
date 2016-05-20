import arcpy
from arcpy import env
from arcpy.sa import *
import pandas as pd
from storm_stats_functions import check_dir
import numpy as np
import matplotlib.pyplot as plt

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


def take_out_gages(df, x, y, k):
    df['dist'] = ((df['x']-x)**2 + (df['y'] - y)**2)**0.5
    df = df.sort_values('dist')
    df = df.reset_index()
    dists = [df.ix[i, 'dist'] for i in range(k)]
    df = df.ix[k:]
    return df, dists


def per_inc(df):
    res0 = df[df.num_removed == 0]
    res1 = df[df.num_removed > 0]
    a = res0['var']
    a.reset_index(inplace=True, drop=True)
    b = res1['var']
    b.reset_index(inplace=True, drop=True)
    inc = ((b-a)/a).mean()
    return inc


def plot_results(df):
    res0 = df[df.num_removed == 0]
    res1 = df[df.num_removed > 0]
    y0 = res0['var']
    y1 = res1['var']
    n = len(y0)

    fig,ax = plt.subplots()
    x0 = np.arange(0.875, 0.875+n)
    x1 = np.arange(1.125, 0.875+n)

    b0 = ax.bar(x0, y0, 0.25)
    b1 = ax.bar(x1, y1, 0.25, color='orange')

    ax.set_xticks(x1)
    ax.set_xticklabels(res1.time_stamp, rotation='vertical')
    ax.set_ylabel(r'Average semi-variance ($mm^2$)')
    ax.set_xlabel("Date/Time")
    ax.set_xlim(0.5, 20.5)
    ax.legend((b0, b1), ("With nearest stations", "Without nearest stations"), loc=0)
    fig.tight_layout()
    plt.show()

# specify type
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
means = []
j = 0
res_df = pd.DataFrame()
for i in np.arange(1):
    i = 5
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

    # take out nearest points
    ks = [0, 2]

    for k in ks:
        if k > 0:
            res = take_out_gages(df, wshed_x, wshed_y, k)
            red_df = res[0]
            dists_removed = res[1]
        else:
            dists_removed = 0
            red_df = df

        for date in dates:
            j += 1
            # convert to excel then to table then to layer then to shapefile... phew!
            check_dir(k_dir)
            date_df = red_df[['site_name', 'x', 'y', 'src', date]]  # get df for individual timestep
            date_df = date_df[date_df[date].notnull()]
            date_df.columns = ['site_name', 'x', 'y', 'src', 'z']
            date_table = "{}temp.xls".format(k_dir)
            date_df.to_excel(date_table, index=False)

            tbl = "tab{}.dbf".format(j)

            arcpy.ExcelToTable_conversion(date_table, tbl, Sheet="Sheet1")

            spref = r"Coordinate Systems/Projected Coordinate Systems/State Plane/Nad 1983/" \
                    r"NAD 1983 HARN StatePlane Virginia South FIPS 4502 (Meters).prj"
            lyr = "Layer{}".format(j)
            arcpy.MakeXYEventLayer_management(tbl, 'x', 'y', lyr, spref)

            rain_shp = 'temp{}.shp'.format(j)
            arcpy.CopyFeatures_management(lyr, rain_shp)

            #  get model params
            mp = ModelParams(tpe, data_dir, date)

            # do kriging
            out_est_file = "rainEst{}".format(j)
            cell_size = "61.3231323600002"
            out_var_file = "var{}".format(j)
            arcpy.gp.Kriging_sa(rain_shp,
                                "z",
                                out_est_file,
                                "{} {} {} {} {}".format(mp.model_type, mp.lag_size, mp.range, mp.sill, mp.nugget),
                                cell_size,
                                "{} {}".format(mp.sample_type, mp.sample_num),
                                out_var_file)

            # get mean of the kriged rain est surface of selected watershed
            out_tab = "out_tab{}".format(j)
            arcpy.sa.ZonalStatisticsAsTable(sel, "FID", out_est_file, out_tab, "DATA", "MEAN")
            cur = arcpy.da.SearchCursor(out_tab, "MEAN")
            for row in cur:
                rain_est = row[0]

            # get mean of the kriged var surface of selected watershed
            out_tab = "out_tab{}".format(j)
            arcpy.sa.ZonalStatisticsAsTable(sel, "FID", out_var_file, out_tab, "DATA", "MEAN")
            cur = arcpy.da.SearchCursor(out_tab, "MEAN")
            for row in cur:
                var_est = row[0]

            res_df = res_df.append({'watershed_descr': wshed_descr,
                                    'time_stamp': date,
                                    'num_removed': k,
                                    'dists': dists_removed,
                                    'est': rain_est,
                                    'var': var_est},
                                   ignore_index=True)

    res_df.to_csv("../../Manuscript/Data/{}_res.csv".format(wshed_descr))
# Todo: put in loop for each of the time steps
# Todo: loop through each station removing one by one and see the effect on the rainfall estimation


