import pandas as pd
from storm_stats_functions import get_data_frame_from_table, data_dir


def overall_summary(ts):
    df = pd.DataFrame()
    for t in ts:
        df1 = pd.read_csv('{0}kriging results/{1}/{1}_summary_single.csv'.format(data_dir, t))
        df = df.append(df1, ignore_index=True)
    df.to_csv('all_summary_single_hz.csv')


def summarize(res0, res1):
    a = res0['var']
    a.reset_index(inplace=True, drop=True)
    b = res1['var']
    b.reset_index(inplace=True, drop=True)
    inc = ((b - a) / a).mean()

    a = res0['est']
    a.reset_index(inplace=True, drop=True)
    b = res1['est']
    b.reset_index(inplace=True, drop=True)
    per_diff = abs(((b - a) / ((a + b) / 2))).mean()
    diff = abs(b-a)

    num_rem = res1.num_removed.mean()
    try:
        stations_removed = res1.iloc[1]['stations_removed']
    except:
        stations_removed = ""
    try:
        dists = res1.dists.iloc[0].replace('[', "")
    except IndexError:
        print 'whatt?'
    dists = dists.replace(']', "")
    dists = dists.split(',')
    dists = [float(d) for d in dists]
    ave_dist = sum(dists) / len(dists)

    station_list = stations_removed.replace("[", "")
    station_list = station_list.replace("]", "")
    station_list = station_list.split("'")
    # station_list = station_list.split("\\")
    sites_df = get_data_frame_from_table('sites_list')
    last_station = station_list[-2]
    src = sites_df[sites_df.site_name == station_list[-2]]['src'].values[0]
    percent_rec_df = get_data_frame_from_table('percent_recorded')
    percent_recorded = percent_rec_df[percent_rec_df.site_name == last_station]['percent_recorded'].values[0]

    return {'per_inc_var': inc,
            'per_diff_rain': per_diff,
            'num_missing': num_rem,
            'max_diff': diff.max(),
            'ave_diff': diff.mean(),
            'ave_dist': ave_dist,
            'max_dist': max(dists),
            'min_dist': min(dists),
            'dists': dists,
            'src': src,
            'stations_removed': stations_removed,
            'percent_recorded': percent_recorded
            }


def individual_summaries(names, types, single):
    for typ in types:
        sum_df = pd.DataFrame(columns=["wshed", "time_step", "per_inc_var", "per_diff_rain", "num_missing", "ave_dist",
                                       "min_dist", "max_dist", "stations_removed", "percent_neg"])
        for n in names:
            f = '{0}kriging results/{1}/{1}_{2}.csv'.format(data_dir, typ, n)
            df = pd.read_table(f, sep=',')

            #change negative estimates to 0
            neg_indices = df.est[df.est < 0].index
            df.est.loc[neg_indices] = 0

            wshd_df = pd.read_csv('{}wshed-ids.csv'.format(data_dir))
            p = wshd_df.loc[wshd_df['Description'] == n]['numrem'].values[0]

            for i in df.num_removed.unique():
                if i > 0:
                    if single:
                        s = summarize(df[df.num_removed == 0], df[df.num_removed == p])
                    else:
                        try:
                            s = summarize(df[df.num_removed == i - 1], df[df.num_removed == i])
                        except TypeError:
                            print n

                    s['id'] = wshd_df[wshd_df['Description'] == n]['ID'].values[0]
                    s['wshed'] = n
                    s['time_step'] = typ
                    s['percent_neg'] = len(df.est[df.est<0])/float(len(df.est))

                    sum_df = sum_df.append(s, ignore_index=True)
                    if single:
                        break
        if single:
            sum_df.to_csv("{0}kriging results/{1}/{1}_summary_single.csv".format(data_dir, typ), index=False)
        else:
            sum_df.to_csv("{0}kriging results/{1}/{1}_summary_multi.csv".format(data_dir, typ), index=False)


names = ['Shore Drive and Kendall Street',
         'Ocean View Ave and Mortons Road',
         'Shore Drive and Great Neck Road',
         '21st and Baltic',
         'Shore Drive and Red Tide Road',
         'S. Rosemont and Clubhouse',
         'S. Rosemont and S. Plaza Trail']


types = ['hr_zeroes']

individual_summaries(names, types, single=True)
overall_summary(types)
