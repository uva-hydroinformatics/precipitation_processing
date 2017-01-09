from precipitation_processing.storm_stats_functions import read_sub_daily, get_daily_tots_df


def get_src_count(s, firstletters):
    return sum(s.index.str.startswith(firstletters))

# Figure out how many stations are null values for the whole day
df = get_daily_tots_df(exclude_zeros=False, qc=True)
df1 = df.ix[:, 3:]
sum_by_station = df1.sum(axis=1)
print "these are the stations that did not record any rainfall for all 20 days:"
print sum_by_station[sum_by_station == 0].index

null_counts = df1[df1 == 0]
null_counts.dropna(axis=1, inplace=True, how="all")
null_counts.dropna(axis=0, inplace=True, how="all")

null_sum = (null_counts == 0).sum(axis=1)
print "here are the null counts"
print null_sum
print "- "*40
print "{} out of {} stations have daily sum of zero for at least one whole day".format(
    len(null_sum),
    len(df1.index)
)

# determine sources
wu_cnt = get_src_count(null_counts, "KVAVIR")
print "{} from WU".format(wu_cnt)
hrsd_cnt = get_src_count(null_counts, "MMPS-")
print "{} from HRSD".format(hrsd_cnt)
print "{} from VAB".format(len(null_counts) - wu_cnt - hrsd_cnt)

print "- "*40
print null_counts
