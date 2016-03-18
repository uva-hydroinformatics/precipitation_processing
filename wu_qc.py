from storm_stats_functions import *

c_df = combine_data_frames()
dr = get_date_range()
dtots = get_daily_tots_df(c_df, dr)
for d in dr:
    print "---- {} -----".format(d)
    df = dtots[dtots[d].notnull()]
    for samp in df.index:
        samp_x = df.loc[samp]['x']
        samp_y = df.loc[samp]['y']
        dists = ((samp_x-df['x'])**2 + (samp_y-df['y'])**2)**0.5
        k = 5
        dists.sort()
        closest_points = dists[1:k+1]
        ave_distance = closest_points.mean()
        ave_diff = abs(dtots.loc[closest_points.index][d].mean() - dtots.loc[samp][d])
        if ave_diff > dtots[d].std()*2:
            print 'site: ', samp
            print 'ave_distance: ', ave_distance
            print 'site value: ', dtots.loc[samp][d]
            print 'mean: ', dtots[d].mean()
            print 'ave_diff: ', ave_diff
            print ' '





