import storm_stats_functions as st
import semivariogram_cj_samp as cj
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def daily_non_zero(daily_tots_df, date):
    df = daily_tots_df[daily_tots_df[date]>0]
    df = df[['x', 'y', 'src', date]]
    return df


def time_step_variogram_df(df, date):
    P = np.array(df[['x','y',date]])
    lag = 1000
    lags = np.arange(0, 20000, lag)
    sv = cj.SV(P, lags, lag)
    svd = pd.DataFrame(sv.transpose())
    return svd


def plot_single_day(df):
    outline_poly = st.get_outline_polygon()
    flavor = 'for_krige'
    base_dir = 'C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/'
    fig_dir = '{}figures/python/{}/'.format(base_dir, flavor)

    st.plot_scatter_subplots(df,
                             units="Daily avg precip (mm)",
                             type=flavor,
                             dty=fig_dir,
                             title="daily_avg_precip",
                             marker_scale=1,
                             threshold=1,
                             label=True,
                             ply=outline_poly)

d = st.get_daily_tots_df(st.get_date_range())
dqc = st.qc_wu(d)
d_precip = dqc.ix[:, 3:]
d_sums = d_precip.sum(1)
d_days = (d_precip>0).sum(1)
d_avg = d_sums/d_days
d_avg_df = dqc.ix[:, :3].join(d_avg.to_frame('d_avg'))


dates = st.get_date_range()
for date in dates:
    a = daily_non_zero(d, date)
    sv = time_step_variogram_df(a, date)
    sv_plot = sv.plot.scatter(0, 1, title=date)
    f = sv_plot.get_figure()
    plt.savefig('scatter{}.png'.format(date))
    plot_single_day(a)
