from storm_stats_functions import *

########################################################################################################################
# Prepare Data##########################################################################################################
########################################################################################################################
combined_df = combine_data_frames()
date_range = get_date_range()
outline_poly = get_outline_polygon()

flavor = 'all data'
base_dir = 'C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/'
fig_dir = '{}figures/python/{}/'.format(base_dir, flavor)
data_dir = '{}data/{}/'.format(base_dir, flavor)


########################################################################################################################
# Plot/summarize data ##################################################################################################
########################################################################################################################
# date_range = ['2014-07-10']

# get daily summary ##
# daily_tots_df = get_daily_tots_df(combined_df, date_range)
#
# get storm durations ##
# durations = get_storm_durations(combined_df, date_range, 0.025)
#
# get subdaily summary ##
# timestep = "15T"
# t = get_subdaily_df(combined_df, date_range, timestep)
# t = combine_sub_daily_dfs(combined_df, t)

# for a in t:
#     a[1].to_csv("{}{}-{}.csv".format(check_dir(data_dir), a[0], 'fifteen_min'))
# plot subdaily data ##
# plot_subdaily_scatter(t,
#                       False,
#                       timestep,
#                       units="Precip (mm)",
#                       type="subdaily",
#                       dty=fig_dir,
#                       marker_scale=5,
#                       threshold=-1,
#                       ply=outline_poly,
#                       label=True)
# #
# # for d in t:
# #     d[1].to_csv("{}{}_{}.csv".format(check_dir(data_dir), timestep, d[0]))
#
# # plot daily data ##
timestep = "4/15/2014 11:15:00"
df = read_sub_daily('fif')
df = df.ix[:, ['x', 'y', 'src', timestep]]
plot_scatter_subplots(df,
                      units="Rainfall (mm)",
                      type=flavor,
                      dty=fig_dir,
                      title="sub_day",
                      marker_scale=15,
                      threshold=0.01,
                      ply=outline_poly)
#
# get and plot intensity data at 15 min step ##
# max_daily_intensities_fifteen = get_daily_max_intensities(combined_df, date_range, "15T")
# plot_scatter_subplots(max_daily_intensities_fifteen,
#                       units="Max daily intensity (mm/15 min)",
#                       type=flavor,
#                       dty=fig_dir,
#                       title="max_daily_intensities_15_min",
#                       marker_scale=1.3,
#                       threshold=1,
#                       ply=outline_poly)

# get and plot intensity data at hour step ##
# max_daily_intensities_hour = get_daily_max_intensities(combined_df, date_range, "H")
# plot_scatter_subplots(max_daily_intensities_hour,
#                       units="Max daily intensity (mm/hour)",
#                       type=flavor,
#                       dty=fig_dir,
#                       title="max_daily_intensities_hour",
#                       marker_scale=1.1,
#                       threshold=1,
#                       ply=outline_poly)


# plot summary scatter + bars ##
# plot_sum_by_station_bars(daily_tots_df, fig_dir, flavor, outline_poly)

# bar graph for mean rainfall for each day ##
plot_sum_by_day(daily_tots_df, 'summary.png')

# compile overall storm summary table ##
# create_summary_table(daily_tots_df,
#                      max_daily_intensities_hour,
#                      max_daily_intensities_fifteen,
#                      data_dir,
#                      "overall_storm_summary")
#