# This code is based off of code written by C. Johnson and was accessed at the following url:
# http://connor-johnson.com/2014/03/20/simple-kriging-in-python/
# Purpose: Create semivariogram figures given point rainfall data from rain gauge networks
# Authors: J. Sadler, M. Morsy, University of Virginia
# Email: jms3fb@virginia.edu

from pylab import *
import numpy as np
import pandas
from scipy.spatial.distance import pdist, squareform



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
                Z.append((P[i, 2] - P[j, 2]) ** 2.0)
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


def create_semivariogram(df, name, date_range, bw, hs):
    print 'doing data for {}'.format(name)
    i = 1
    sv_combined = list()
    for date in date_range:
        print 'analyzing data for {}'.format(str(date))
        df_for_date = df[df['date'] == date]
        if len(df_for_date) == 0:
            continue
        P = np.array(df_for_date[['x', 'y', 'rain']])
        P = P[~np.isnan(P).any(axis=1)]
        sv = SV(P, hs, bw)
        #add sv to combined semivariogram
        if i == 1:
            sv_combined.append(sv[0])
        sv_combined.append(sv[1])
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
    savefig('semivariogram_figure/semivariogram_model {}.png'.format(name), fmt='png', dpi=200, bbox_inches='tight')
    close()

    #plot the combined semivariogram
    sv_combined_array = np.array(sv_combined)
    sv_array_average = np.mean(sv_combined_array[1:, :], axis=0)
    plot(sv_combined_array[0], sv_array_average, '.--', lw=0.8, ms=4)
    nugget = sv[1][0]
    model_results = [r + nugget for r in sp(sv[0])]
    plot(sv[0], model_results, color='r', lw=0.8)
    xlim(0, hs.max())
    xlabel('Lag [m]', fontsize=4)
    ylabel('Semivariance', fontsize=4)
    title('average semi-variogram for {}'.format(name))
    tick_params(labelsize=4, length=2)
    grid( color='0.65')
    savefig('semivariogram_figure/semivariogram_model_AVE {}.png'.format(name), fmt='png', dpi=200, bbox_inches='tight')
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
              20151002]

df_list = []
path = 'C:/Users/jeff_dsktp/Box Sync/Sadler_1stPaper/rainfall/data/'
file_name_list = ['wu_daily_sum', 'vabeach_daily_sum_reformat_in']
ext = 'csv'
for file_name in file_name_list:
    df = pandas.read_csv('{}{}.{}'.format(path, file_name, ext))
    df_list.append({'df': df, 'name': file_name})
combined_df = pandas.DataFrame(columns=['x','y','rain'])
for df in df_list:
    combined_df = combined_df.append(df['df'], ignore_index=True)
df_list.append({'df': combined_df, 'name': 'combined'})

for df in df_list:
    P = np.array(df['df'][['x', 'y', 'rain']])
    P = P[~np.isnan(P).any(axis=1)]
    # bandwidth, plus or minus bw meters. Here we are using the 10th percentile
    pd = squareform(pdist(P[:, :2]))
    bw = np.percentile(pd[np.nonzero(pd)], 10)
    # lags in bw meter increments from zero to the max distance
    hs = np.arange(0, np.amax(pd), bw)
    create_semivariogram(df['df'], df['name'], date_range, bw, hs)
