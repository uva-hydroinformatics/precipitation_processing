from storm_stats_functions import *

def idw(neigh):
    neigh['wv'] = neigh.val/neigh.dist
    pv = neigh.wv.sum()/sum(1/neigh.dist)
    return(pv)

def qc(sub_day_df):
    for date in sub_day_df:
        date_df = date[1]
        for d in date_df.columns[47:]:
            df = date_df[date_df[d].notnull()]
            print "---- {} -----".format(d)
            for samp in df.index:
                df = date_df[date_df[d].notnull()]
                samp_x = df.loc[samp]['x']
                samp_y = df.loc[samp]['y']
                dists = ((samp_x-df['x'])**2 + (samp_y-df['y'])**2)**0.5
                k = 5
                dists.sort()
                closest_neigh_dist = dists[1]
                closest_points = dists[1:1+k]
                ave_distance = closest_points.mean()
                neigh = df.loc[closest_points.index][d]
                dists = pd.concat([closest_points, neigh], axis=1, keys=['dist', 'val'])
                pred_val = idw(dists)
                neigh_ave = neigh.mean()
                neigh_std = neigh.std()
                samp_prcp = df.loc[samp][d]
                neigh_diff = abs(neigh_ave - samp_prcp)
                is_outlier = False
                if neigh_diff > neigh_std*2:
                    is_outlier = True
                if samp_prcp == 0 and neigh_ave > 10:
                    is_outlier = True
                if is_outlier:
                    print 'site: ', samp
                    print 'shortest_dist', closest_neigh_dist
                    print 'ave_distance: ', ave_distance
                    print 'site value: ', samp_prcp
                    print 'neigh mean: ', neigh_ave
                    print 'neigh std: ', neigh_std
                    print 'neigh_diff: ', neigh_diff
                    print 'global mean: ', df[d].mean()
                    print ' '


c_df = combine_data_frames()
dr = get_date_range()
dtots = get_daily_tots_df(c_df, dr)
sub_day = get_subdaily_df(c_df, [dr[6]], "15T")
qc(sub_day)




