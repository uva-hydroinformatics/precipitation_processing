# Purpose: Aggregate data from access database for given time interval (originally day)
# Authors: J. Sadler, University of Virginia
# Email: jms3fb@virginia.edu

import numpy as np
import pandas as pd
import pyodbc
import datetime
import matplotlib.pyplot as plt
import math


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
    df = pd.DataFrame(a, columns=[i[0] for i in cur.description])
    cur.close()
    con.close()
    return df


def make_incremental(df, date_range):
    newdf = pd.DataFrame()
    non_cumulative_df = pd.DataFrame()
    for date in date_range:
        print "making incremental for {}".format(date)
        date_string = datetime.datetime.strptime(str(date), '%Y%m%d').strftime('%Y-%m-%d')
        df_date = df[date_string]
        df_date = df_date.groupby(['x', 'y'])
        for location in df_date:
            xy_df = location[1]
            cum_precip_arr = np.array(xy_df['precip_mm'])
            incr_precip = [0]
            for i in range(len(cum_precip_arr)-1):
                if cum_precip_arr[i+1] >= cum_precip_arr[i] and cum_precip_arr[i]>0 :
                    incr_precip.append(cum_precip_arr[i+1] - cum_precip_arr[i])
                else:
                    non_cumulative_df = non_cumulative_df.append(location[1].iloc[[i]])
                    incr_precip.append(0)
            xy_df.loc[:, 'precip_mm'] = incr_precip
            newdf = newdf.append(xy_df)
    newdf = newdf.reset_index(inplace=False)
    non_cumulative_df.to_csv("{}.csv".format("non_cumulative"), mode='a')
    return newdf

def autolabel(ax, rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        if math.isnan(height): continue
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%d' % int(height),
                ha='center', va='bottom')

def graph_bars(ax, x, y, title, ylab, xlabs):
    print ("graphing bars for {}".format(title))
    width = 1.
    rects = ax.bar(x, y, width)
    ax.set_ylabel(ylab)
    ax.set_xticks(x+width/2)
    ax.set_xticklabels(xlabs, rotation='vertical')
    ax.set_title(title)
    autolabel(ax, rects)

def graph_scatter(ax, x, y, title, scale, c_limits, marker_scale):
    print ("graphing scatter for {}".format(title))
    sc = ax.scatter(x, y, c=scale, cmap='Blues', s=scale**marker_scale, vmin=c_limits[0], vmax=c_limits[1])
    ax.set_title(title, fontsize=8)
    ax.tick_params(labelsize=8)
    ax.locator_params(nbins=5)
    return sc

def get_daily_aggregate(df, date, time_step):
    # return a dataframe with the sum of the rainfall at a given point for a given time span
    df_date = df[date]
    df_gp = df_date.groupby(['x', 'y', 'site_name', 'src'])
    df_agg = df_gp.resample(time_step, how={'precip_mm': 'sum'})
    return df_agg

def get_daily_tots_df(summary_df, df, date_range):
    for date in date_range:
        date = datetime.datetime.strptime(str(date), '%Y%m%d').strftime('%Y-%m-%d')
        daily_tot = get_daily_aggregate(df, date, "D")
        daily_tot = daily_tot.sum(level="site_name") #todo: add code to this so we can look at different time intervals (hourly, 15 min)

        #add to summary dataframe
        daily_tot.rename(columns={'precip_mm': date}, inplace=True)
        summary_df = summary_df.join(daily_tot[date])
    return summary_df

def plot_scatter_subplots(df, units, type, dir, title, marker_scale):
    fig, a = plt.subplots(10, 2, figsize=(10,10), sharex=True, sharey = True)
    a = a.ravel()
    #scatter subplots for all days
    c_limits = (df.iloc[:, 3:].min().min(), df.iloc[:, 3:].max().max())
    i = 0
    for col in df.columns[3:]:
        #check if value is below a certain level in mm, leave it out
        threshold = 1 #mm
        d= df[daily_tots_df[col] > threshold]

        sc = graph_scatter(a[i],
                           d.x/1000,
                           d.y/1000,
                           col,
                           d[col],
                           c_limits,
                           marker_scale)
        i += 1
    cax = fig.add_axes([0.915, 0.1, 0.025, 0.8])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label(units)
    fig.text(0.05, .5, "y [km]", rotation="vertical")
    fig.text(.5, .05, "x [km]")
    plt.tick_params(labelsize=10)
    plt.subplots_adjust(wspace=0.1, hspace=.3)
    plt.savefig("{}{}_{}".format(dir, title, type))
    # plt.show()


def get_storm_durations(df, date_range, trim_percent):
    i = 0
    for date in date_range:
        date = datetime.datetime.strptime(str(date), '%Y%m%d').strftime('%Y-%m-%d')
        df_agg = get_daily_aggregate(df, date, "15T")
        time_sums = df_agg.sum(level="datetime")
        trimmed_indices = time_sums[(time_sums['precip_mm']>trim_percent) & (time_sums['precip_mm']<(1-trim_percent))].index
        start_time = trimmed_indices.min().time()
        end_time = trimmed_indices.max().time()
        duration = (trimmed_indices.max()-trimmed_indices.min()).total_seconds()/3600
        if i == 0:
            dur_df = pd.DataFrame({'start_time': [start_time], 'end_time':[end_time], 'date': [date], 'duration (hr)': [duration]})
        else:
            dur_df = dur_df.append({'start_time': start_time, 'end_time':end_time, 'date': date, 'duration (hr)': duration}, ignore_index=True)
        i=1
    dur_df = dur_df.set_index("date")
    return dur_df

def get_daily_max_intensities(summary_df, df, date_range, time_step):
    for date in date_range:
        date = datetime.datetime.strptime(str(date), '%Y%m%d').strftime('%Y-%m-%d')
        df_agg = get_daily_aggregate(df, date, "15T")
        if time_step == 'H':
            df_agg = pd.rolling_sum(df_agg, 4).max(level='site_name')
        else:
            df_agg = df_agg.max(level="site_name")

        #add to summary dataframe
        df_agg.rename(columns={'precip_mm': date}, inplace=True)
        summary_df = summary_df.join(df_agg[date])
    return summary_df


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
df['datetime'] = pd.to_datetime(df['datetime'])
vab_df = df.set_index('datetime')
vab_df.insert(len(vab_df.columns), 'src', 'vab')

df = get_data_frame_from_table('hrsd_obs_spatial')
df['datetime'] = pd.to_datetime(df['datetime'])
hrsd_df = df.set_index('datetime')
hrsd_df.insert(len(hrsd_df.columns), 'src', 'hrsd')
hrsd_df = hrsd_df[hrsd_df.site_name.str.rstrip().isin(hrsd_stations_in_study_area)]

df = get_data_frame_from_table('wu_inc')
df['datetime'] = pd.to_datetime(df['datetime'])
wu_df = df.set_index('datetime')
wu_df.insert(len(wu_df.columns), 'src', 'wu')
wu_df.drop('site_name', axis=1, inplace=True)
wu_df.rename(columns={'site_code':'site_name'}, inplace=True)

df_list = [wu_df, hrsd_df, vab_df]

#combine the dfs in the list together
combined_df = pd.DataFrame()
for df in df_list:
    combined_df = combined_df.append(df)

#create an empty dataframe with just the site_names, xs, ys, and srcs to fill in the summary data
empty_daily_tots_df = combined_df.loc[:,['site_name', 'x', 'y', 'src']].set_index('site_name').drop_duplicates()

type = 'with_wu_filt'
fig_dir = 'C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/figures/python/storm_stats/{}/filt/'.format("with_wu")

# get the summary statistics
daily_tots_df = get_daily_tots_df(empty_daily_tots_df, combined_df, date_range)
plot_scatter_subplots(daily_tots_df, "Total cummulative precip (mm)", type, fig_dir, "total_cummulative_precip", 1)

max_daily_intensities_fifteen = get_daily_max_intensities(empty_daily_tots_df, combined_df, date_range, "15T")
plot_scatter_subplots(max_daily_intensities_fifteen, "Max daily intensity (mm/15 min)", type, fig_dir, "max_daily_intensities_15_min", 1.3)

max_daily_intensities_hour = get_daily_max_intensities(empty_daily_tots_df, combined_df, date_range, "H")
plot_scatter_subplots(max_daily_intensities_hour, "Max daily intensity (mm/hour)", type, fig_dir, "max_daily_intensities_hour", 1.1)

durations = get_storm_durations(combined_df, date_range, 0.05)

#bar chart for overall sum by stations
fig, ax = plt.subplots(1, 2, figsize=(10, 10))
sum_by_station = daily_tots_df.iloc[:,3:].sum(axis=1)
barX = np.arange(len(sum_by_station))
graph_bars(ax[0], barX, sum_by_station, title="", xlabs=sum_by_station.index, ylab="precip (mm)")

#scatter for overall sum by station
graph_scatter(ax[1],
              daily_tots_df.x,
              daily_tots_df.y,
              title="",
              scale=sum_by_station,
              c_limits=(sum_by_station.min(),
                        sum_by_station.max()),
              marker_scale=.6)
fig.suptitle("Summary by station")
plt.tight_layout()
plt.savefig("{}{}_{}.png".format(fig_dir, "overall_summary_by_station", type))

#bar graph for mean rainfall for each day
fig, ax = plt.subplots()
daily_totals = daily_tots_df.iloc[:,3:].mean()
x = np.arange(len(daily_totals))
y = daily_totals
graph_bars(ax, x, y, title="overall summary by date", xlabs=daily_totals.index, ylab="cumulative total precip (mm)")
plt.tight_layout()
plt.savefig('{}{}_{}.png'.format(fig_dir, "overall_summary_by_date", type))

#compile overall storm summary table
overall_summary_df_by_date = durations.join(pd.DataFrame({'mean_total_rainfall_volume (mm)': daily_tots_df.mean()}))
overall_summary_df_by_date['average_intensity (mm/hr)'] = overall_summary_df_by_date['mean_total_rainfall_volume (mm)']/overall_summary_df_by_date['duration (hr)']
overall_summary_df_by_date = overall_summary_df_by_date.join(pd.DataFrame({'mean_max_hourly_intensity (mm/hr)': max_daily_intensities_hour.mean()}))
overall_summary_df_by_date = overall_summary_df_by_date.join(pd.DataFrame({'max_max_hourly_intensity (mm/hr)': max_daily_intensities_hour.max()}))
overall_summary_df_by_date = overall_summary_df_by_date.join(pd.DataFrame({'mean_max_15min_intensity (mm/15 min)': max_daily_intensities_fifteen.mean()}))
overall_summary_df_by_date = overall_summary_df_by_date.join(pd.DataFrame({'max_max_15min_intensity (mm/15 min)': max_daily_intensities_fifteen.max()}))

#write to csv file
overall_summary_df_by_date.to_csv("C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/data/overall_summary_by_date.csv")
