# This code is based off of code written by C. Johnson and was accessed at the following url:
# http://connor-johnson.com/2014/03/20/simple-kriging-in-python/
# Purpose: Create semivariogram figures given point rainfall data from rain gauge networks
# Authors: J. Sadler, M. Morsy, University of Virginia
# Email: jms3fb@virginia.edu

from pylab import *
import numpy as np
import pandas
from scipy.spatial.distance import pdist, squareform
import pyodbc


def plot_stations(z):
    fig, ax = subplots()
    ax.scatter(z.x, z.y, c=z.rain, cmap='gray')
    ax.set_aspect(1)
    # xlim(-1500,22000)
    # ylim(-1500,17500)
    xlabel('Easting [m]')
    ylabel('Northing [m]')
    title('Rain (in)')
    fig.savefig('C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfallgoogle.png')


def SVh(P, h, bw):
    """
    Experimental semivariogram for a single lag
    """
    pd = squareform(pdist(P[:, :2]))
    N = pd.shape[0]
    Z = list()
    for i in range(N):
        for j in range(i + 1, N):
            if (pd[i, j] >= h - bw) and (pd[i, j] <= h + bw):
                Z.append(((P[i, 2] - P[j, 2]) ** 2.0))
    return np.sum(Z) / (2.0 * len(Z))


def SV(P, hs, bw):
    """
    Experimental variogram for a collection of lags
    """
    sv = list()
    for h in hs:
        sv.append(SVh(P, h, bw))
    sv = [[hs[i], sv[i]] for i in range(len(hs)) if sv[i] > 0]
    return np.array(sv).T

def C( P, h, bw ):
    """
    Calculate the sill
    """
    c0 = np.var( P[:,2] )
    if h == 0:
        return c0
    return c0 - SVh( P, h, bw )

def opt(fct, x, y, C0, parameterRange=None, meshSize=1000):
    if parameterRange == None:
        parameterRange = [x[1], x[-1]]
    mse = np.zeros(meshSize)
    a = np.linspace(parameterRange[0], parameterRange[1], meshSize)
    for i in range(meshSize):
        mse[i] = np.mean((y - fct(x, a[i], C0)) ** 2.0)
    return a[mse.argmin()]


def spherical(h, a, C0):
    """
    Spherical model of the semivariogram
    """
    # if h is a single digit
    if type(h) == np.float64:
        # calculate the spherical function
        if h <= a:
            return C0 * (1.5 * h / a - 0.5 * (h / a) ** 3.0)
        else:
            return C0
    # if h is an iterable
    else:
        # calcualte the spherical function for all elements
        a = np.ones(h.size) * a
        C0 = np.ones(h.size) * C0
        return map(spherical, h, a, C0)


def cvmodel(P, model, hs, bw):
    """
    Input:  (P)      ndarray, data
            (model)  modeling function
                      - spherical
                      - exponential
                      - gaussian
            (hs)     distances
            (bw)     bandwidth
    Output: (covfct) function modeling the covariance
    """
    # calculate the semivariogram
    sv = SV(P, hs, bw)
    # calculate the sill
    C0 = C(P, hs[0], bw)
    # calculate the optimal parameters
    param = opt(model, sv[0], sv[1], C0)
    # return a covariance function
    covfct = lambda h, a=param: model(h, a, C0) #changed from original "C0 - model( h, a, C0 )
    return covfct


def get_data_frame_from_table(table_name):

    # set up db connection
    MDB = "C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/data/rainfall_data_master.accdb"; DRV = '{Microsoft Access Driver (*.mdb, *.accdb)}'; PWD = 'pw'

    # connect to db
    con = pyodbc.connect('DRIVER={};DBQ={};PWD={}'.format(DRV,MDB,PWD))
    cur = con.cursor()

    # run a query and get the results
    SQL = 'SELECT * FROM {};'.format(table_name) # your query goes here
    rows = cur.execute(SQL).fetchall()
    a = np.array(rows)
    df = pandas.DataFrame(a, columns=[i[0] for i in cur.description])
    cur.close()
    con.close()
    return df

def aggregate_time_steps(df_date, hours, date, wu):
    # return a dataframe with the sum of the rainfall at a given point for a given time span
    date_string = datetime.datetime.strptime(str(date), '%Y%m%d').strftime('%Y-%m-%d')
    df_date = df[date_string]
    df_date = df_date.groupby(['x', 'y'])
    df_date = df_date.resample('D', how={'precip_mm': 'sum'})
    df_date = df_date.reset_index(inplace=False)
    if wu == True:
        newdf = pandas.DataFrame()
        for group in df_date:
            xy_df = group[1]
            cum_precip_arr = np.array(xy_df['precip_mm'])
            if cum_precip_arr[0]>0:
                incr_precip = [cum_precip_arr[0]] #TODO: decide whether or not to keep this or make it zero to begin with
            else:
                incr_precip = [0]
            for i in range(len(cum_precip_arr)-1):
                if cum_precip_arr[i+1] >= cum_precip_arr[i]:
                    incr_precip.append(cum_precip_arr[i+1] - cum_precip_arr[i])
                else:
                    incr_precip.append(cum_precip_arr[i+1])
            xy_df.loc[:, 'precip_mm'] = incr_precip
            newdf = newdf.append(xy_df)
        newdf = newdf.groupby(['x', 'y', 'site_name'])
        newdf = newdf.resample('D', how={'precip_mm': 'sum'})
        newdf = newdf.reset_index(inplace=False)
        return newdf
    else:
        return df_date


def create_semivariogram(df, name, date_range, bw, hs):
    print 'doing data for {}'.format(name)
    i = 1
    for date in date_range:
        print 'creating semivariagram for {}'.format(str(date))
        dfdates = df['datetime']
        df_for_date = df[df['datetime'] == datetime.datetime.strptime(str(date), "%Y%m%d")]
        if len(df_for_date) == 0:
            continue
        P = np.array(df_for_date[['x', 'y', 'precip_mm']], dtype=float)
        P = P[~np.isnan(P).any(axis=1)]
        sv = SV(P, hs, bw)
        #add sv to combined semivariogram
        if i == 1:
            ave_df = sv
        else:
            for w in range(len(sv[0])):
                for j in range(len(ave_df[0])):
                    if sv[0][w] == ave_df[0][j]:
                        ave_df[1][j] = (ave_df[1][j] + sv[1][w])*0.5
        #create the model
        sp = cvmodel(P, spherical, hs, bw)
        #plot the data and the model
        subplot(10, 2, i)
        plot(sv[0], sv[1]/np.var(P[:, 2]), '.--', lw=0.8, ms=4)
        nugget = sv[1][0]
        model_results = [r + nugget for r in sp(sv[0])]
        plot(sv[0], model_results/np.var(P[:, 2]), color='r', lw=0.8)
        xlim(0, hs.max())
        # ylim(-1500,17500)
        # figure(figsize=(7, 4))
        if i == 19 or i == 20:
            xlabel('Lag [m]', fontsize=4)
        if i % 2 == 1:
            ylabel('Semivariance', fontsize=4)
        title((datetime.datetime.strptime(str(date), '%Y%m%d').strftime('%m/%d/%y')), fontsize=5)
        tick_params(labelsize=4, length=2)
        grid( color='0.65')

        i += 1
    tight_layout(pad=0.1, w_pad=0.1, h_pad=0.05)
    savefig('C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/semivariogram_figure/semivariogram_model 3{}.png'
            .format(name),
            fmt='png',
            dpi=400,
            bbox_inches='tight')
    close()

    #plot the combined semivariogram
    sp = cvmodel(P, spherical, hs, bw)
    plot(ave_df[0], ave_df[1]/np.var(P[:, 2]), '.--', lw=0.8, ms=4)
    nugget = 0
    model_results = [r + nugget for r in sp(sv[0])]
    plot(sv[0], model_results/np.var(P[:, 2]), color='r', lw=0.8)
    xlim(0, hs.max())
    xlabel('Lag [m]', fontsize=4)
    ylabel('Semivariance', fontsize=4)
    title('average semi-variogram for {}'.format(name))
    tick_params(labelsize=4, length=2)
    grid( color='0.65')
    savefig('C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/semivariogram_figure/semivariogram_model_AVE 3{}.png'
            .format(name),
            fmt='png',
            dpi=400,
            bbox_inches='tight')
    close()


date_range = [20130702,
              20131009,
              20140111,
              20140213,
              20140415,
              20140425,
              20140710,
              20140818,
              20140908,
              20140909,
              20140913,
              20141126,
              20141224,
              20150414,
              20150602,
              20150624,
              20150807,
              20150820,
              20150930,
              20151002
              ]

df_list = []
table_name_list = ['vabeach_reformat_mm']
print 'getting data from database for: {}'.format(table_name_list[0])
df = get_data_frame_from_table(table_name_list[0])
df['datetime'] = pandas.to_datetime(df['datetime'])
df = df.set_index('datetime')
df_list.append(df)
# df_list = [get_data_frame_from_table(table) in table_name_list]
# combined_df = pandas.DataFrame(columns=['x','y','rain'])
# for df in df_list:
#     combined_df = combined_df.append(df, ignore_index=True)
# df_list.append({'df': combined_df, 'name': 'combined'})

comb_df = pandas.DataFrame()
for date in date_range:
    print 'aggregating data for: {}'.format(str(date))
    comb_df = comb_df.append(aggregate_time_steps(df, 24, date, wu=False))
P = np.array(df[['x', 'y', 'precip_mm']], dtype=np.float)
P = P[~np.isnan(P).any(axis=1)]
# bandwidth, plus or minus bw meters. Here we are using the 10th percentile
# pd = squareform(pdist(P[:, :2]))
# bw = np.percentile(pd[np.nonzero(pd)], 10)
bw = 3700
# lags in bw meter increments from zero to the max distance

hs = np.arange(0, 18000, bw)
create_semivariogram(comb_df, 'vab', date_range, bw, hs)
