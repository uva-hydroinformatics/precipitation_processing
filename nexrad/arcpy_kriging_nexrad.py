import arcpy
from arcpy import env
import pandas as pd
from precipitation_processing.storm_stats_functions import check_dir, get_data_frame_from_table, data_dir
import shutil
import os
import psutil
import stat
import datetime


class ModelParams:
    def __init__(self, tpe, d):
        if tpe == 'fifteen_min':
            tpe = 'fif'
        table = "{}_model_params".format(tpe)
        param_df = get_data_frame_from_table(table)
        param_df.set_index('date', inplace=True)

        self.model_type = param_df['type'][d]
        self.lag_size = "1000.000000"
        self.range = param_df['range'][d]
        self.sill = param_df['sill'][d]
        self.nugget = "0"  # todo get this from R script
        self.sample_num = "27"  # todo make this dynamic?
        self.krad = arcpy.sa.RadiusVariable(27, 20000)
        self.sample_type = "Variable"


def take_out_gages(df, x, y, k):
    df['dist'] = ((df['x'] - x) ** 2 + (df['y'] - y) ** 2) ** 0.5
    df = df.sort_values('dist')
    df = df.reset_index()
    dists = [df.ix[i, 'dist'] for i in range(k)]
    stations_removed = [df.ix[i, 'site_name'] for i in range(k)]
    df = df.ix[k:]
    return df, dists, stations_removed


def change_permissions_recursive(path):
    for root, dirs, files in os.walk(path, topdown=False):
        for dir in [os.path.join(root,d) for d in dirs]:
            os.chmod(dir, stat.S_IWRITE)
    for file in [os.path.join(root, f) for f in files]:
            os.chmod(file, stat.S_IWRITE)


def onerror(func, path):
    """
    Error handler for ``shutil.rmtree``.

    If the error is due to an access error (read only file)
    it attempts to add write permission and then retries.

    If the error is for another reason it re-raises the error.

    Usage : ``shutil.rmtree(path, onerror=onerror)``
    """
    import stat
    if not os.access(path, os.W_OK):
        # Is the error an access error ?
        os.chmod(path, stat.S_IWRITE)
        func(path)
    else:
        raise


def clearWSLocks(inputWS):
    '''Attempts to clear ArcGIS/Arcpy locks on a workspace.

    Two methods:
     1: if ANOTHER process (i.e. ArcCatalog) has the workspace open, that process is terminated
     2: if THIS process has the workspace open, it attempts to clear locks using arcpy.Exists, arcpy.Compact and arcpy.Exists in sequence

    Notes:
     1: does not work well with Python Multiprocessing
     2: this will kill ArcMap or ArcCatalog if they are accessing the worspace, so SAVE YOUR WORK

    Required imports: os, psutil
    '''

    # get process ID for this process (treated differently)
    thisPID = os.getpid()

    # normalise path
    _inputWS = os.path.normpath(inputWS)

    # get list of currently running Arc/Python processes
    p_List = []
    ps = psutil.process_iter()
    for p in ps:
        if ('Arc' in p.name()) or ('python' in p.name()):
            p_List.append(p.pid)

            # iterate through processes
    for pid in p_List:
        p = psutil.Process(pid)

        # if any have the workspace open
        if any(_inputWS in pth for pth in [fl.path for fl in p.open_files()]):
            print '      !!! Workspace open: %s' % _inputWS

            # terminate if it is another process
            if pid != thisPID:
                print '      !!! Terminating process: %s' % p.name
                p.terminate()
            else:
                print '      !!! This process has workspace open...'

                # if this process has workspace open, keep trying while it is open...
    while any(_inputWS in pth for pth in [fl.path for fl in psutil.Process(thisPID).open_files()]):
        print '    !!! Trying Exists, Compact, Exists to clear locks: %s' % all(
            [arcpy.Exists(_inputWS), arcpy.Compact_management(_inputWS), arcpy.Exists(_inputWS)])

    return True


def raster_mean_under_plygn(plygn, raster):
    out_tab = "out_tab"
    arcpy.sa.ZonalStatisticsAsTable(plygn, "FID", raster, out_tab, "DATA", "MEAN")
    cur = arcpy.da.SearchCursor(out_tab, "MEAN")
    for row in cur:
        val = row[0]
    arcpy.Delete_management(out_tab)
    try:
        return val
    except UnboundLocalError:
        return None


def make_date_shapefile(date, df):
    # convert to excel then to table then to layer then to shapefile... phew!
    date_df = df[['site_name', 'x', 'y', 'src', date]]  # get df for individual timestep
    date_df = date_df[date_df[date].notnull()]
    date_df.columns = ['site_name', 'x', 'y', 'src', 'z']
    global k_dir
    date_table = "{}/temp.xls".format(k_dir)
    date_df.to_excel(date_table, index=False)

    tbl = "tab.dbf"

    arcpy.ExcelToTable_conversion(date_table, tbl, Sheet="Sheet1")

    spref = r"Coordinate Systems/Projected Coordinate Systems/State Plane/Nad 1983/" \
            r"NAD 1983 HARN StatePlane Virginia South FIPS 4502 (Meters).prj"
    lyr = "Layer"
    arcpy.MakeXYEventLayer_management(tbl, 'x', 'y', lyr, spref)

    shp = 'temp.shp'
    arcpy.CopyFeatures_management(lyr, shp)
    arcpy.Delete_management(date_table)
    arcpy.Delete_management(tbl)
    return shp


def get_nexrad_file(ts):
    ts_utc = ts + datetime.timedelta(hours=5)
    # check for daylight savings:
    if (datetime.datetime(2016, 03, 13) > ts_utc > datetime.datetime(2015, 11, 01))\
            or (datetime.datetime(2015, 03, 8) > ts_utc > datetime.datetime(2014, 11, 02)) \
            or (datetime.datetime(2014, 03, 8) > ts_utc > datetime.datetime(2013, 11, 03)):
        ts_utc += datetime.timedelta(hours=1)

    d = os.path.join(data_dir, "nexrad/{}".format(ts.strftime('%Y%m%d')))
    files = os.listdir(d)
    nexrad_time_strings = [t.split("_")[-2] + t.split("_")[-1] for t in files]
    nexrad_time_strings_cln = [t.split(".")[0] for t in nexrad_time_strings]
    nexrad_times = [datetime.datetime.strptime(t, "%Y%m%d%H%M%S") for t in nexrad_time_strings_cln]
    time_diff = datetime.timedelta(hours=1)
    for i in range(len(nexrad_times)):
        diff = ts_utc - nexrad_times[i]
        if abs(time_diff.total_seconds()) > abs(diff.total_seconds()):
            time_diff = diff
            nex_file = files[i]
    try:
        return os.path.join(d, nex_file)
    except UnboundLocalError:
        return None


# specify type; should be 'fifteen_min', 'hr', or 'daily'
tpe = 'hr'

# get data from table
df = get_data_frame_from_table('hr')
a = df.ix[:, 4:]
non_zero_dates = a.columns[a.sum() > 0]

# set up arcpy environment
arcpy.CheckOutExtension("spatial")
env.extent = arcpy.Extent(3705690, 1051630, 3724920, 1068584)
k_dir = 'C:/Users/jeff_dsktp/Google Drive/Hampton Roads GIS Data/VA_Beach_Data/kriging1'
check_dir(k_dir)
change_permissions_recursive(k_dir)
os.chmod(k_dir, stat.S_IWRITE)
shutil.rmtree(k_dir, onerror=onerror, ignore_errors=False)
check_dir(k_dir)
change_permissions_recursive(k_dir)

env.overwriteOutput = True
env.workspace = k_dir

# get watershed polygon and iterate through them
shed_ply = "C:/Users/jeff_dsktp/Documents/Research/Sadler_1st_Paper/Manuscript/Data/GIS/problem_watersheds.shp"
arcpy.MakeTableView_management(shed_ply, "table_view")
num_rows = int(arcpy.GetCount_management("table_view").getOutput(0))
means = []
j = 0
res_df = pd.DataFrame(columns=['watershed_descr',
                               'time_stamp',
                               'est',
                               'var',
                               'nex_est',
                               'nex_file'])

# do it these watersheds (according to arcid)
hd = False  # makes it so there is a header for each file
for date in non_zero_dates:
    j += 1
    timestamp = datetime.datetime.strptime(date, "%Y-%m-%d %H:%M:%S")

    check_dir(k_dir)
    rain_shp = make_date_shapefile(date, df)

    #  get model params
    mp = ModelParams(tpe, date)

    # do kriging
    out_est_file = "rainEst"
    cell_size = "61.3231323600002"
    out_var_file = "var"
    arcpy.gp.Kriging_sa(rain_shp,
                        "z",
                        out_est_file,
                        "{} {} {} {} {}".format(mp.model_type, mp.lag_size, mp.range, mp.sill, mp.nugget),
                        cell_size,
                        "{} {}".format(mp.sample_type, mp.sample_num),
                        out_var_file)
    nexrad_raster = "nex_prjd.tif"
    nexrad_file_name = get_nexrad_file(timestamp)
    if nexrad_file_name:
        arcpy.ProjectRaster_management(nexrad_file_name, nexrad_raster, out_est_file, "BILINEAR")
    for i in [0, 1, 2, 3, 4, 5, 6]:
        # select individual watershed
        sel = "selection{}.shp".format(i)
        arcpy.Select_analysis(shed_ply, sel, '"FID" = {}'.format(i))
        cur = arcpy.SearchCursor(sel)
        for row in cur:
            wshed_descr = row.getValue("Descript")
            wshed_area = row.getValue("Area_sq_km")
            wshed_x = row.getValue("x")
            wshed_y = row.getValue("y")
        if wshed_area < 0.001:
            continue

        # get mean of the kriged rain est surface of selected watershed
        rain_est = raster_mean_under_plygn(plygn=sel, raster=out_est_file)

        # get mean of the kriged var surface of selected watershed
        var_est = raster_mean_under_plygn(plygn=sel, raster=out_var_file)

        # get mean of nexrad data
        if nexrad_file_name:
            nexrad_est_in = raster_mean_under_plygn(plygn=sel, raster=nexrad_raster)
            nexrad_file_name = nexrad_file_name.split("\\")[-1]
        else:
            nexrad_est_mm = None
        nexrad_est_mm = nexrad_est_in * 25.4 if nexrad_est_in else None

        a = res_df.append({'watershed_descr': wshed_descr,
                           'time_stamp': date,
                           'est': rain_est,
                           'nex_est': nexrad_est_mm,
                           'nex_file': nexrad_file_name,
                           'var': var_est},
                          ignore_index=True)

        a.to_csv("../../Data/kriging results/{}/nexrad/{}_{}_nex.csv".format(tpe, tpe, wshed_descr),
                 mode='a',
                 header=hd,
                 index=False)
        # clearWSLocks(k_dir)
        print "name: {}    date:{}    num_removed:{}".format(wshed_descr, date, 0)
        arcpy.Delete_management(sel)
    print "time step: {}".format(j)
    hd = False

    arcpy.Delete_management(nexrad_raster)
    arcpy.Delete_management(out_var_file)
    arcpy.Delete_management(out_est_file)
    arcpy.Delete_management(rain_shp)
