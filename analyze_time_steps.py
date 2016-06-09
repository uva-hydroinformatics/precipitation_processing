# Purpose: Takes a time series and prints some summary stats about the time intervals. Originally designed to do this
# for rainfall data from weatherunderground
# Authors: J. Sadler, M. Morsy, University of Virginia
# Email: jms3fb@virginia.edu
# Original date: 1/13/2016

import pandas as pd
import numpy as np
from storm_stats_functions import get_data_frame_from_table, qc_wu, get_date_range

dr = get_date_range()
# connect to db
table_name = "all_data"
# run a query and get the results
df = get_data_frame_from_table(table_name)
df = df[df.src=='WU']
df = qc_wu(df)
df.reset_index(inplace=True)
df.datetime = pd.to_datetime(df.datetime)
intervals = list()
for n in df.site_name.unique():
    df1 = df[df.site_name == n]
    for d in dr:
        df1 = df1.set_index(df1.datetime)
        try:
            df2 = df1[d]
            for i in range(len(df2.index)-1):
                intervals.append((df2.index[i + 1] - df2.index[i]).seconds)
        except:
            continue

# get intervals in seconds
intervals_arr = np.array(intervals)
print("mean: {}".format(str(np.mean(intervals_arr))))
print("max: {}".format(str(np.max(intervals_arr))))
print("min: {}".format(str(np.min(intervals_arr))))
print("median: {}".format(str(np.median(intervals_arr))))
print("std: {}".format(str(np.std(intervals_arr))))
