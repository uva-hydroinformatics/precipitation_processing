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
import matplotlib.patheffects as pe


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


def opt(fct, x, y, C0, parameterRange=None, meshSize=1000): #TODO evaulate sensitivity of meshsize, see what a change in hs does
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
    print 'getting data for {}'.format(table_name)
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


def make_incremental(df, date_range):
    newdf = pandas.DataFrame()
    for date in date_range:
        date_string = datetime.datetime.strptime(str(date), '%Y%m%d').strftime('%Y-%m-%d')
        df_date = df[date_string]
        us = len(df_date['site_name'].unique())
        df_date = df_date.groupby(['x', 'y'])
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
    return newdf


def aggregate_time_steps(df, hours, date, wu):
    # return a dataframe with the sum of the rainfall at a given point for a given time span
    date_string = datetime.datetime.strptime(str(date), '%Y%m%d').strftime('%Y-%m-%d')
    df_date = df[date_string]
    df_date = df_date.groupby(['x', 'y'])
    df_date = df_date.resample('D', how={'precip_mm': 'sum'})
    df_date = df_date.reset_index(inplace=False)
    # this poriton of the code changes cumulative values to incremental
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


def create_semivariograms(df, semivar_name, date_range):
    i = 1
    for date in date_range:
        print "creating sv for {}".format(str(date))
        indivi_df = aggregate_time_steps(df, 24, date, wu=False)

        #create array for sv and kriging
        P = create_array(indivi_df)
        pd = squareform(pdist(P[:, :2]))
        bw = np.percentile(pd[np.nonzero(pd)], 10)
        hs = np.arange(0, np.max(pd), bw)
        sv_results = SV(P, hs, bw)
        model = cvmodel(P, spherical, hs, bw)
        subplot(10, 2, i)
        plot_semivariogram(sv_results[0], sv_results[1], model(sv_results[0]), str(date), hs.max(), np.var(P[:, 2]))
        if i == 19 or i == 20:
            xlabel('Lag [m]', fontsize=4)
        if i % 2 == 1:
            ylabel('Semivariance', fontsize=4)
        title((datetime.datetime.strptime(str(date), '%Y%m%d').strftime('%m/%d/%y')), fontsize=5)
        tick_params(labelsize=4, length=2)
        grid(color='0.65')

        i += 1
    tight_layout(pad=0.1, w_pad=0.1, h_pad=0.05)
    savefig('C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/semivariogram_figure/{}.png'
            .format(semivar_name),
            fmt='png',
            dpi=400,
            bbox_inches='tight')


def plot_semivariogram(x, y, model_results, title_text, xmax, var):
    plot(x, y/var, '.--', lw=0.8, ms=4)
    nugget = y[0]
    model_results = [(r + nugget)/var for r in model_results]
    plot(x, model_results, color='r', lw=0.8)
    xlim(0, xmax)
    title('average semi-variogram for {}'.format(title_text))


def create_array(df):
        P = np.array(df[['x', 'y', 'precip_mm']], dtype=float)
        P = P[~np.isnan(P).any(axis=1)]
        return P


def create_interpolations(P, x_step, y_step, num_neigbors):
    X0, X1 = P[:, 0].min(), P[:, 0].max()
    Y0, Y1 = P[:, 1].min(), P[:, 1].max()
    Z = np.zeros((y_step, x_step))
    dx, dy = (X1-X0)/x_step, (Y1-Y0)/y_step
    for i in range(y_step):
        print i,
        for j in range(x_step):
            Z[i, j] = krige(P, spherical, hs, bw, (indivi_df.x.min()+dx*j, indivi_df.y.min()+dy*i), num_neigbors)
    return Z

def plot_interpolation(x_list, y_list, z_list, Z, font_size, plot_limits):
    marg = 500
    scale = 1000
    x_list_scaled = [i / scale for i in x_list]
    y_list_scaled = [i / scale for i in y_list]
    scatter(x_list_scaled, y_list_scaled, c=z_list, linewidths=0.75, s=20, cmap=get_cmap('Oranges'))
    for i in range(len(z_list)):
        labs = annotate(z_list[i], (x_list_scaled[i],
                             y_list_scaled[i]),
                             fontsize=font_size,
                             weight='bold',
                             # xytext=(1, 1) TODO figure out how do offset labels
                 )
    Z = np.flipud(Z)
    imshow(Z, interpolation='nearest', extent=[min(x_list_scaled), max(x_list_scaled), min(y_list_scaled), max(y_list_scaled)])
    set_cmap("Blues")
    tick_params(labelsize=font_size)
    cbar = colorbar(shrink=0.9)
    cbar.ax.tick_params(labelsize=font_size)
    xlim((plot_limits['xmin'] - marg)/scale, ((plot_limits['xmax'] + marg)/scale))
    ylim((plot_limits['ymin'] - marg)/scale, ((plot_limits['ymax'] + marg)/scale))
    # tight_layout(pad=0.1, w_pad=0.1, h_pad=0.05)
    title('{}'.format(date), fontsize=4, fontdict={'verticalalignment':'bottom'})


def krige(P, model, hs, bw, u, N):
    '''
    Input  (P)     ndarray, data
           (model) modeling function
                    - spherical
                    - exponential
                    - gaussian
           (hs)    kriging distances
           (bw)    kriging bandwidth
           (u)     unsampled point
           (N)     number of neighboring
                   points to consider
    '''

    # covariance function
    covfct = cvmodel( P, model, hs, bw )
    # mean of the variable
    mu = np.mean( P[:,2] )

    # distance between u and each data point in P
    d = np.sqrt( ( P[:,0]-u[0] )**2.0 + ( P[:,1]-u[1] )**2.0 )
    # add these distances to P
    P = np.vstack(( P.T, d )).T
    # sort P by these distances
    # take the first N of them
    P = P[d.argsort()[:N]]

    # apply the covariance model to the distances
    k = covfct( P[:,3] )
    # cast as a matrix
    k = np.matrix( k ).T

    # form a matrix of distances between existing data points
    K = squareform( pdist( P[:,:2] ) )
    # apply the covariance model to these distances
    K = covfct( K.ravel() )
    # re-cast as a NumPy array -- thanks M.L.
    K = np.array( K )
    # reshape into an array
    K = K.reshape(N,N)
    # cast as a matrix
    K = np.matrix( K )

    # calculate the kriging weights
    weights = np.linalg.inv( K ) * k
    weights = np.array( weights )

    # calculate the residuals
    residuals = P[:,2] - mu

    # calculate the estimation
    estimation = np.dot( weights.T, residuals ) + mu

    return float( estimation )




date_range = [
              20130702,
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

# prepare the data by pulling from the database and making the datetime the index
df_list = []
df = get_data_frame_from_table('vabeach_reformat_mm')
df['datetime'] = pandas.to_datetime(df['datetime'])
df_vab = df.set_index('datetime')
df_list.append(df)

df = get_data_frame_from_table('wu_observation_spatial')
df['datetime'] = pandas.to_datetime(df['datetime'])
df = df.set_index('datetime')
inc_df = make_incremental(df, date_range)
df_list.append(inc_df)
df_wu = inc_df

#combine the dfs in the list together
combined_df = pandas.DataFrame()
for df in df_list:
    combined_df = combined_df.append(df)

create_semivariograms(combined_df, "testing again", date_range)

a = 1
for date in date_range[:1]:
    print 'creating interpolation for {}'.format(str(date))
    indivi_df = aggregate_time_steps(combined_df, 24, date, wu=False)
    P = create_array(indivi_df)
    pd = squareform(pdist(P[:, :2]))
    bw = np.percentile(pd[np.nonzero(pd)], 10)
    hs = np.arange(0, np.max(pd), bw)
    Z = create_interpolations(P, x_step=50, y_step=50, num_neigbors=8)
    subplot(4,3,a)
    plot_limits = {
        'xmin': indivi_df.x.min(),
        'xmax': indivi_df.x.max(),
        'ymin': indivi_df.y.min(),
        'ymax': indivi_df.y.max()
    }
    plot_interpolation(P[:, 0].tolist(), P[:, 1].tolist(), P[:, 2].tolist(), Z, font_size=4, plot_limits=plot_limits)
    tight_layout(pad=0.1, w_pad=0.01, h_pad=0.05)
    a += 1
inter_file_name = "testing again vab"
savefig('C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/kriging/{}.png'
        .format(inter_file_name),
        fmt='png',
        dpi=800)
close()



