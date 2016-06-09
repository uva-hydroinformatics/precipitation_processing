from storm_stats_functions import qc_wu, read_sub_daily
import pandas as pd

# read table from db. 'fif', 'hr', or 'daily'
df = read_sub_daily('fif')

# perform qc for wu stations
# df = qc_wu(df)

#save as excel file
df.to_excel('../data/fifteen_min_no_qc.xlsx', index=True)