import pandas as pd
from storm_stats_functions import read_sub_daily

times = ['daily', 'fif', 'hr']
for t in times:
    df = read_sub_daily(t)
    df = df.ix[:, 4:]
    non_zero_dates = df.columns[df.sum() > 0]
    df = df[non_zero_dates]
    mtostd = (df.std()/df.mean()).mean()
    print "std to mean ratio for {}: {}".format(t, mtostd)
