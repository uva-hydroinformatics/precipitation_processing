# Purpose: Aggregate data from access database for given time interval (originally day)
# Authors: J. Sadler, University of Virginia
# Email: jms3fb@virginia.edu

import numpy as np
import pandas
import pyodbc
import datetime
import matplotlib.pyplot as plt


def get_data_frame_from_table(table_name):
    print 'getting data for {}'.format(table_name)
    # set up db connection
    MDB = "C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/data/rainfall_data_master.accdb"; DRV = '{Microsoft Access Driver (*.mdb, *.accdb)}'; PWD = 'pw'

    # connect to db
    con = pyodbc.connect('DRIVER={};DBQ={};PWD={}'.format(DRV,MDB,PWD))
    cur = con.cursor()

    # run a query and get the results
    SQL = 'SELECT * FROM {};'.format(table_name) # your query goes here
    rows = cur.execute(SQL).fetchall()
    a = np.array(rows)
    df = pandas.DataFrame(a, columns=[i[0] for i in cur.description])
    cur.close()
    con.close()
    return df


def make_incremental(df, date_range):
    newdf = pandas.DataFrame()
    non_cumulative_df = pandas.DataFrame()
    for date in date_range:
        date_string = datetime.datetime.strptime(str(date), '%Y%m%d').strftime('%Y-%m-%d')
        df_date = df[date_string]
        df_date = df_date.groupby(['x', 'y'])
        for location in df_date:
            xy_df = location[1]
            cum_precip_arr = np.array(xy_df['precip_mm'])
            if cum_precip_arr[0] > 0:
                incr_precip = [cum_precip_arr[0]] #TODO: decide whether or not to keep this or make it zero to begin with
            else:
                incr_precip = [0]
            for i in range(len(cum_precip_arr)-1):
                if cum_precip_arr[i+1] >= cum_precip_arr[i]:
                    incr_precip.append(cum_precip_arr[i+1] - cum_precip_arr[i])
                else:
                    non_cumulative_df = non_cumulative_df.append(location[1].iloc[[i]])
                    incr_precip.append(cum_precip_arr[i+1])
            xy_df.loc[:, 'precip_mm'] = incr_precip
            newdf = newdf.append(xy_df)
    newdf = newdf.reset_index(inplace=False)
    non_cumulative_df.to_csv("{}.csv".format("non_cumulative"), mode='a')
    return newdf


def graph_daily_totals(df_sums, date, ylab, type):
    width = 1.
    fig,ax = plt.subplots()
    x = np.arange(len(df_sums.precip_mm))
    rects = ax.bar(x, df_sums.precip_mm, width)
    ax.set_ylabel(ylab)
    ax.set_xticks(x+width/2)
    ax.set_xticklabels(df_sums.site_name +" (" + df_sums.src + ")", rotation='vertical')
    ax.set_title(date)
    plt.tight_layout()
    plt.savefig('C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/python/storm_stats/{}{}{}.png'.format(type,"_",date ))

def get_daily_total(df, date):
    # return a dataframe with the sum of the rainfall at a given point for a given time span
    df_date = df[date]
    df_date = df_date.groupby(['x', 'y', 'site_name', 'src'])
    df_sums = df_date.sum()
    df_sums = df_sums.reset_index(inplace=False)
    return df_sums

def get_daily_max(df,date):
    graph_daily_totals(df_sums, date_string, "cumulative precip in mm")
    df_date = df_date.resample('15T', how={'precip_mm': 'sum'})
    df_date = df_date.reset_index(inplace=False)


hrsd_stations_in_study_area = ["MMPS-171",
                               "MMPS-185",
                               "MMPS-163",
                               "MMPS-255",
                               "MMPS-146",
                               "MMPS-004",
                               "MMPS-256",
                               "MMPS-140",
                               "MMPS-160",
                               "MMPS-144",
                               "MMPS-036",
                               "MMPS-093-2"]


date_range = [
              20130702,
              20131009,
              20140111,
              20140213,
              20140415,
              20140425,
              20140710,
              20140818,
              20140908,
              20140909,
              20140913,
              20141126,
              20141224,
              20150414,
              20150602,
              20150624,
              20150807,
              20150820,
              20150930,
              20151002
              ]

# prepare the data by pulling from the database and making the datetime the index
df = get_data_frame_from_table('vabeach_reformat_mm')
df['datetime'] = pandas.to_datetime(df['datetime'])
vab_df = df.set_index('datetime')
vab_df.insert(len(vab_df.columns), 'src', 'vab')

df = get_data_frame_from_table('hrsd_obs_spatial')
df['datetime'] = pandas.to_datetime(df['datetime'])
hrsd_df = df.set_index('datetime')
hrsd_df.insert(len(hrsd_df.columns), 'src', 'hrsd')

# df = get_data_frame_from_table('wu_observation_spatial')
# df['datetime'] = pandas.to_datetime(df['datetime'])
# df = df.set_index('datetime')
# inc_df = make_incremental(df, date_range)
# wu_df = inc_df.set_index('datetime')
# wu_df.insert(len(wu_df.columns), 'src', 'wu')
#
df_list = [hrsd_df, vab_df]

#combine the dfs in the list together
combined_df = pandas.DataFrame()
for df in df_list:
    combined_df = combined_df.append(df)

combined_df = combined_df[(combined_df.src == "vab") | (combined_df.site_name.str.rstrip().isin(hrsd_stations_in_study_area))]
daily_tots_df = pandas.DataFrame()
daily_means = []
for i in range(len(date_range)):
    date = date_range[i]
    date = datetime.datetime.strptime(str(date), '%Y%m%d').strftime('%Y-%m-%d')
    daily_tot = get_daily_total(combined_df, date)
    graph_daily_totals(daily_tot, date, "cumulative precip (mm)","daily total",)
    daily_means.append((date, daily_tot.precip_mm.mean()))

daily_means = zip(*daily_means)
width = 1.
fig,ax = plt.subplots()
x = np.arange(len(daily_means[0]))
rects = ax.bar(x, daily_means[1], width)
ax.set_ylabel("cumulative ave precip (mm)")
ax.set_xticks(x+width/2)
ax.set_xticklabels(daily_means[0], rotation='vertical')
ax.set_title("overall summary")
plt.tight_layout()
plt.savefig('C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/python/storm_stats/{}.png'.format("summary"))
