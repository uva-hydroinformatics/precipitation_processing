
import pyodbc
import pandas
import numpy as np
# set up db connection
MDB = "C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/data/rainfall_data_master.accdb"; DRV = '{Microsoft Access Driver (*.mdb, *.accdb)}'; PWD = 'pw'

# connect to db
con = pyodbc.connect('DRIVER={};DBQ={};PWD={}'.format(DRV,MDB,PWD))
cur = con.cursor()

table_name = "wu_observation"
# run a query and get the results
SQL = 'SELECT * FROM {};'.format(table_name) # your query goes here
rows = cur.execute(SQL).fetchall()
a = np.array(rows)
df = pandas.DataFrame(a, columns=[i[0] for i in cur.description])
times = df['datetime']
intervals = list()
for i in range(len(times)-1):
    intervals.append((times[i+1] - times[i]).seconds)
intervals_arr = np.array(intervals)
print("mean: {}".format(str(np.mean(intervals_arr))))
print("max: {}".format(str(np.max(intervals_arr))))
print("min: {}".format(str(np.min(intervals_arr))))
print("median: {}".format(str(np.median(intervals_arr))))
print("std: {}".format(str(np.std(intervals_arr))))
cur.close()
con.close()
