from storm_stats_functions import *


def idw(neigh):
    neigh['wv'] = neigh.val/neigh.dist
    pv = neigh.wv.sum()/sum(1/neigh.dist)
    return pv


def analyze_from_neighbors(df, outlier_check):
    indiv_df = pd.DataFrame()
    summary_df = pd.DataFrame()
    stations = df.index
    for station in stations:
        print "---- {} -----".format(station)
        count = 0
        for d in df.columns[3:]:
            if not math.isnan(df.loc[station][d]):
                count += 1
                filt_df = df[df[d].notnull()]
                if outlier_check:
                    filt_df = filt_df[filt_df['src'] != 'wu']
                samp_x = df.loc[station]['x']
                samp_y = df.loc[station]['y']

                dists = ((samp_x-filt_df['x'])**2 + (samp_y-filt_df['y'])**2)**0.5
                k = 3
                dists.sort_values(inplace=True)
                closest_points = dists[dists < 5000]
                closest_points = closest_points[1:]
                num_neighs = closest_points.shape[0]
                if num_neighs < k:
                    closest_points = dists[1:1+k]
                closest_neigh_dist = closest_points[0]
                furthest_neigh_dist = closest_points[-1]
                ave_distance = closest_points.mean()

                neigh = filt_df.loc[closest_points.index][d]
                dists = pd.concat([closest_points, neigh], axis=1, keys=['dist', 'val'])
                pred_val = idw(dists)
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
                indiv_df = indiv_df.append({'datetime': d,
                                            'true_station': station,
                                            'recorded_val': samp_prcp,
                                            'neighbors': closest_points.index.tolist(),
                                            'num_neighbors': num_neighs,
                                            'shortest_dist': closest_neigh_dist,
                                            'longest_dist': furthest_neigh_dist,
                                            'ave_distance': ave_distance,
                                            'pred_val': pred_val,
                                            'global_mean': filt_df[d].mean(),
                                            'percent_diff': percent_diff,
                                            'index': percent_diff/ave_distance,
                                            'src': df.loc[station]['src'],
                                            'out': out
                                            },
                                           ignore_index=True)
        s_df = indiv_df[indiv_df['true_station'] == station]
        if outlier_check:
            n_out = sum(s_df.out != 'NA')
            n_poss = float(s_df.shape[0])
            indivi_summ = {'n_out': n_out, 'n_poss': n_poss, 'percent_out': (n_out/n_poss),
                           'neighbors': closest_points.index.tolist(), 'num_neighbors': num_neighs,
                           'shortest_dist': closest_neigh_dist, 'longest_dist': furthest_neigh_dist,
                           'ave_distance': ave_distance}
        else:
            indivi_summ = error_summary(s_df)
        indivi_summ['station'] = station
        indivi_summ['src'] = df.loc[station]['src']
        summary_df = summary_df.append(indivi_summ, ignore_index=True)
    if outlier_check:
        indiv_df = indiv_df[indiv_df['out'] != 'NA']
    return indiv_df, summary_df


def error_summary(error_df):
    mpf = error_df.percent_diff[(error_df.percent_diff > 0) & (np.isfinite(error_df.percent_diff))].mean()
    maxpf = error_df.percent_diff[(error_df.percent_diff > 0) & (np.isfinite(error_df.percent_diff))].max()
    dif = abs(error_df.pred_val-error_df.recorded_val)
    mdif = dif.mean()
    maxdif = dif.max()
    summary = {'mean percent error': mpf,
               'max percent error': maxpf,
               'mean mm diff': mdif,
               'max mm diff': maxdif}
    return summary

fmin = read_sub_daily('fif')
c = analyze_from_neighbors(fmin, True)


