from storm_stats_functions import *


def idw(neigh):
    neigh['wv'] = neigh.val/neigh.dist
    pv = neigh.wv.sum()/sum(1/neigh.dist)
    return pv


def qc(df, outlier_check):
    container_df = pd.DataFrame()
    stations = df.index
    stations = ['11 BEACH BOROUGH']
    for station in stations:
        print "---- {} -----".format(station)
        for d in df.columns[3:]:
            if not math.isnan(df.loc[station][d]):
                filt_df = df[df[d].notnull()]
                filt_df = filt_df[filt_df['src'] != 'wu']
                samp_x = df.loc[station]['x']
                samp_y = df.loc[station]['y']

                dists = ((samp_x-filt_df['x'])**2 + (samp_y-filt_df['y'])**2)**0.5
                k = 3
                dists.sort_values(inplace=True)
                closest_neigh_dist = dists[1]
                closest_points = dists[1:1+k]
                ave_distance = closest_points.mean()

                neigh = filt_df.loc[closest_points.index][d]
                dists = pd.concat([closest_points, neigh], axis=1, keys=['dist', 'val'])
                pred_val = idw(dists)
                neigh_ave = neigh.mean()
                neigh_std = neigh.std()
                samp_prcp = df.loc[station][d]
                neigh_diff = abs(pred_val - samp_prcp)
                if outlier_check:
                    percent_diff = neigh_diff*2/(samp_prcp+pred_val)
                else:
                    percent_diff = neigh_diff/samp_prcp
                out = 'NA'
                if neigh_diff > neigh_std*3:
                    out = 'a'
                if samp_prcp == 0 and pred_val > 10:
                    out = 'b'
                container_df = container_df.append({'datetime': d,
                                                    'true_station': station,
                                                    'recorded_val': samp_prcp,
                                                    'neighbors': closest_points.index.tolist(),
                                                    'shortest_dist': closest_neigh_dist,
                                                    'ave_distance': ave_distance,
                                                    'pred_val': pred_val,
                                                    'global_mean': filt_df[d].mean(),
                                                    'percent_diff': percent_diff,
                                                    'index': percent_diff/ave_distance,
                                                    'src': df.loc[station]['src'],
                                                    'out': out
                                                    },
                                                   ignore_index=True)
    if outlier_check:
        container_df = container_df[container_df['out'] != 'NA']
    return container_df

comb_sub_day = read_sub_daily('15_min_all')
comb_sub_day = qc_wu(comb_sub_day)
c = qc(comb_sub_day, False)





