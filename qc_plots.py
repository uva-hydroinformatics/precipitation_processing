from storm_stats_functions import get_data_frame_from_table
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def bar_summary(df):
    fig, ax = plt.subplots()
    df.sort_values('percent_out', inplace=True)
    df.reset_index(inplace=True, drop=True)

    cs = ['moccasin', 'cornflowerblue', 'midnightblue']
    hs = ['','.','/']
    srcs = ['vab', 'hrsd', 'wu']
    for i in range(len(srcs)):
        d = df[df.src == srcs[i]]
        x = d.index
        y = d['percent_out']

        if srcs[i] == 'vab':
            lab = 'City of Virginia Beach'
        elif srcs[i] == 'hrsd':
            lab = 'Hampton Roads Sanitation District'
        elif srcs[i] == 'wu':
            lab = 'Weather Underground'

        ax.bar(x, y, color=cs[i], label=lab)
    outlier_df = df[df.station == 'KVAVIRGI52']
    ax.annotate('KVAVIRGI52',
                xy=(outlier_df.index[0], outlier_df.percent_out*0.95),
                xycoords='data',
                xytext=(32, 0.27),
                textcoords='data',
                arrowprops=dict(facecolor='black', shrink=0.05),
                verticalalignment='mid'
                )
    ax.set_xticks(np.arange(0.5, 52.5, 1))
    ax.set_xticklabels("", rotation='vertical')
    # ax.set_xticklabels(df.station, rotation='vertical')
    ax.set_ylabel('Percent of measurements as outliers')
    ax.set_xlabel('Station')
    ax.set_xlim(0, len(df.station.unique()))
    ax.legend(prop={'size': 12}, loc=0)
    fig.tight_layout()
    plt.show()


def boxplot(df):
    a = df.percent_out.plot.box()
    plt.show()


def scatter_dist(df):
    """
    produces scatter plots of outliers vs distance and colored by station type
    :param df:
    :param dist_type: should be string 'longest_dist', 'ave_distance', or 'shortest_dist'
    :return: void
    """
    dist_types =['longest_dist', 'ave_distance', 'shortest_dist']
    for dist_type in dist_types:
        dfv = df[df.src == 'vab']
        dfh = df[df.src == 'hrsd']
        dfw = df[df.src == 'wu']

        ms = ['s', '^', 'o']
        cs = ['sandybrown', 'cornflowerblue', 'midnightblue']
        i = 0
        scatters = []
        for d in dfv, dfh, dfw:
            x = d[dist_type]
            y = d['percent_out']
            scatters.append(plt.scatter(x, y, marker=ms[i], color=cs[i], s=50))
            i += 1

        if dist_type == 'ave_distance':
            lab = "Average distance to nearest neighbors"
        elif dist_type == 'longest_dist':
            lab = "Longest distance to nearest neighbors"
        elif dist_type == 'shortest_dist':
            lab = "Shortest distance to nearest neighbors"

        plt.ylabel('Percent of measurements as outliers')
        plt.xlabel('{} ($m$)'.format(lab))
        plt.ylim(0)
        plt.xlim(min(df[dist_type]-100))
        plt.legend(scatters, ['City of Virginia Beach', 'Hampton Roads Sanitation District', 'Weather Underground'],
                   prop={'size': 12},
                   scatterpoints=1)
        plt.tight_layout()
        plt.show()

# read in data
df = get_data_frame_from_table('qc_summary')
bar_summary(df)
# boxplot(df)
# scatter_dist(df)