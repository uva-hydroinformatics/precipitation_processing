from storm_stats_functions import get_data_frame_from_table
import pandas as pd
import matplotlib.pyplot as plt


def bar_summary(df):
    df.sort_values('percent_out', inplace=True)
    df = df.set_index('station')
    a = df.percent_out.plot.bar()
    a.set_ylabel('Percent of measurements as outliers')
    a.set_xlabel('Station')
    f = a.get_figure()
    f.tight_layout()
    plt.show()
    st = df[df.src != 'wu']
    sf = df[df.src == 'wu']


def boxplot(df):
    a = df.percent_out.plot.box()
    plt.show()


def scatter_dist(df, dist_type):
    """
    produces scatter plots of outliers vs distance and colored by station type
    :param df:
    :param dist_type: should be string 'longest_dist', 'ave_distance', or 'shortest_dist'
    :return: void
    """
    dfv = df[df.src == 'vab']
    dfh = df[df.src == 'hrsd']
    dfw = df[df.src == 'wu']

    ms = ['s', '^', 'o']
    cs = ['limegreen', 'cornflowerblue', 'purple']
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
# bar_summary(df)
boxplot(df)
# scatter_dist(df, 'longest_dist')