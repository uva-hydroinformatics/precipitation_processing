import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def fix_lab(l):
    if l == '0':
        lab = '15 Min'
    elif l == '1':
        lab = 'Hr'
    elif l == '2':
        lab = 'Daily'
    return lab


def plot_bars(df):
    vars = ['per_diff_rain', 'per_inc_var']
    for v in vars:
        fig, ax = plt.subplots()
        cols = ['lightpink', 'indianred', 'maroon']
        for i in df.t_index.unique():
            d = df[df.t_index == i]
            x = d.index
            y = abs(d[v])
            lab = fix_lab(i)
            ax.bar(x, y, color=cols[int(i)], label=lab)
        ax.set_xticks(np.arange(1.5, 21.5, 3))
        ax.set_xticklabels([str(j) for j in range(1, 8)])
        ax.set_xlabel('Watershed', fontsize=18)
        ax.legend(loc=0)

        if v == 'per_diff_rain':
            ylab = 'Average absolute % difference \nin rainfall estimation'
        if v == 'per_inc_var':
            ylab = 'Average % increase in variance'

        ax.set_ylabel(ylab, multialignment='center', fontsize=18)
        ax.set_xlim(0,21)
        ax.tick_params(labelsize=18)
        fig.tight_layout()
        plt.show()


def time_index_col(row):
    ts = row['time_step']
    if 'daily' in ts:
        val = '2'
    elif 'hr' in ts:
        val = '1'
    elif '15' in ts:
        val = '0'
    return val


def wshed_index_col(row):
    watershed_names = ['Shore Drive and Kendall Street',
             'Ocean View Ave and Mortons Road',
             'Shore Drive and Great Neck Road',
             '21st and Baltic',
             'Shore Drive and Red Tide Road',
             'S. Rosemont and Clubhouse',
             'S. Rosemont and S. Plaza Trail']

    nm = row['wshed']
    val = watershed_names.index(nm)
    return val


def summarize_by_time_scale(df):
    fig,ax = plt.subplots()
    cols = ['lightpink', 'indianred', 'maroon']
    for t in df.t_index.unique():
        df1 = df[df.t_index == t]
        a = df1['per_inc_var'].mean()
        b = df1['per_diff_rain'].mean()
        ts = df1.iloc[0]['time_step']
        print "Average per_in_var for {}: {}".format(ts, a)
        print "Average per_diff_rain for {}: {}".format(ts, b)
        lab = fix_lab(t)
        ax.bar(int(t), a, color=cols[int(t)])
        ax.bar(int(t)+3, b, color=cols[int(t)], label=lab)
    ax.set_xticks([1.5, 4.5])
    ax.set_xticklabels(['Percent increase \n in variance', 'Percent difference in \n rainfall estimation'])
    ax.legend(loc=0)
    ax.tick_params(labelsize=18)
    plt.show()

data_dir = "../../Manuscript/Data/kriging results/"
filename = "{}all.xlsx".format(data_dir)
df = pd.read_excel(filename)
df['t_index'] = df.apply(time_index_col, axis=1)
df['w_index'] = df.apply(wshed_index_col, axis=1)
df.sort_values(['w_index', 't_index'], inplace=True)
df.reset_index(inplace=True, drop=True)
plot_bars(df)
# summarize_by_time_scale(df)
