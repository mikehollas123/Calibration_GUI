from Calibration_analysis import *
import os, sys
import csv
from appJar import gui
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import webbrowser
from scipy import stats
import seaborn as sns
import numpy as np

from pylab import *
from scipy.optimize import curve_fit


from matplotlib.colors import ListedColormap, LinearSegmentedColormap
def GetmeansStd(df):
    dict_meansvsstd = {}

    df_total_std = df
    df_total_std = df_total_std.groupby(["Charge"])['Slope'].std()

    df_total_means = df
    df_total_means = df_total_means.groupby(["Charge"])['Slope'].mean()
    df_total_Ns = df
    Ns = df_total_Ns.groupby(["Charge"])['Slope'].count()

    df_total_error = df
    df_total_error = df_total_error.groupby(["Charge"])['Slope'].sem()

    """    _95conf = [stats.norm.interval(0.95, df_total_means.values[i], df_total_error.values[i] ) for i in
     range(len(df_total_means.values))]"""

    _95conf = [stats.norm.interval(0.95, df_total_means.values[i], (df_total_error.values[i] / np.sqrt(Ns.values[i]))) for i in range(len(df_total_means.values))]


    dict_meansvsstd['total'] = {
        'means':df_total_means.values,
        'means_index':df_total_means.index,
        'std': df_total_std.values,
        'std_index': df_total_std.index,
        '3std+':df_total_means.values + df_total_std.values*3,
        '3std-': df_total_means.values - df_total_std.values * 3,
        '95+':[_95conf[i][1] for i in range(len(_95conf))],
        '95-': [_95conf[i][0] for i in range(len(_95conf))]
    }


    charge_groups = getChargeGroups(df)

    for j in range(len(dict_key_map)):
        min_charge = min(charge_groups[j])
        max_charge = max(charge_groups[j])

        dict_meansvsstd[dict_key_map[j]] = {}
        group_df = df.loc[(df['Charge'] <= max_charge) & (df['Charge'] >= min_charge)]

        df_means = group_df
        df_means = df_means.groupby(["Charge"])['Slope'].mean()

        df_std = group_df
        df_std = df_std.groupby(["Charge"])['Slope'].std()

        dict_meansvsstd[dict_key_map[j]]['means'] = df_means.values
        dict_meansvsstd[dict_key_map[j]]['means_index']= df_means.index
        dict_meansvsstd[dict_key_map[j]]['std'] = df_std.values
        dict_meansvsstd[dict_key_map[j]]['std_index'] = df_std.index

    return dict_meansvsstd
def getChargeGroups(df):
    global dict_key_map
    dict_key_map = list(set(df['standard']))
    charge_groups = []
    for each in dict_key_map:

        charge_groups.append( list(set(df.loc[df['standard'] == each,'Charge'])))



    return charge_groups
def CreateDF(basepath):
    global df
    global progress
    calibration_list = []


    dirs = os.listdir(basepath)
    total_count = len(dirs)
    count = 0
    for entry in dirs:
        try:
            standard = entry.split("_")[0]
        except:
            standard = "unknown"
        if os.path.isfile(os.path.join(basepath, entry)):
            with open(os.path.join(basepath, entry), newline = '') as calibrations:
                calibration_reader = csv.DictReader(calibrations, delimiter='\t')
                for calibration in calibration_reader:
                    del calibration['Time of Death']
                    del calibration['R Squared']
                    del calibration['Frequency']
                    del calibration['m/z']
                    del calibration['Intensity']
                    del calibration['Noise']
                    del calibration['Baseline']
                    del calibration['Resolution']
                    del calibration['File']

                    dict_temp = dict(calibration)
                    dict_temp['standard'] = standard
                    calibration_list.append(dict_temp)

        count += 1

    df =  pd.DataFrame(calibration_list)

    df[['Charge', 'Slope']] = df[['Charge', 'Slope']].apply(pd.to_numeric, errors='coerce')

    return df
def linear_regressions(df,sample_numbers=300):
    slope_grand, intercept_grand, r_value_grand, p_value_grand, std_err_grand = stats.linregress(df['Charge'],
                                                                                                 df['Slope'])
    dict_stats = {}
    dict_stats['total'] = {'slope_grand': slope_grand, 'intercept_grand': intercept_grand,
                           'r_value_grand': r_value_grand, 'p_value_grand': p_value_grand,
                           'std_err_grand': std_err_grand, 'charge': df['Charge'],
                           'slope': df['Slope']}

    charge_groups = getChargeGroups(df)



    for j in range(len(dict_key_map)):
        dict_stats[dict_key_map[j]] = {}
        min_charge = min(charge_groups[j])
        max_charge = max(charge_groups[j])

        group_slope = df.loc[(df['Charge'] <= max_charge) & (df['Charge'] >= min_charge), 'Slope']
        group_charge = df.loc[(df['Charge'] <= max_charge) & (df['Charge'] >= min_charge), 'Charge']

        slope, intercept, r_value, p_value, std_err = stats.linregress(group_charge, group_slope)

        # subsample of data using seaborn - if sample number greater than dataset just use dataset
        try:
            sub_slope = group_slope.sample(sample_numbers)
            sub_charge = group_charge.sample(sample_numbers)
            dict_stats[dict_key_map[j]]['sub_slope'] = sub_slope
            dict_stats[dict_key_map[j]]['sub_charge'] = sub_charge
        except:
            dict_stats[dict_key_map[j]]['sub_slope'] = group_slope
            dict_stats[dict_key_map[j]]['sub_charge'] = group_charge
        dict_stats[dict_key_map[j]]['slope'] = slope
        dict_stats[dict_key_map[j]]['intercept'] = intercept
        dict_stats[dict_key_map[j]]['r_value']=r_value
        dict_stats[dict_key_map[j]]['p_value'] = p_value
        dict_stats[dict_key_map[j]]['std_err'] = std_err
    return dict_stats

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), stats.sem(a)
    h = se * stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


def twoDhist(df):

    from pylab import percentile


    cmap = plt.cm.Reds
    norm = colors.PowerNorm(gamma=0.5)
    fig = plt.figure()

    ax = fig.add_subplot()

    print(np.percentile((df['Charge'], df['Slope']),95,axis=0))

    ax.hist2d(df['Charge'], df['Slope'], bins=(500, 500), norm=norm, cmap=cmap)

    ax.set_xlabel("Charge")
    ax.set_ylabel("Slope")
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
    X = df['Charge']

    ax.plot(X, dict_stats['total']['intercept_grand'] + dict_stats['total']['slope_grand'] * X, 'black', linewidth=0.2)
    #ax.plot(X,[dict_stats['total']['intercept_grand'] + dict_stats['total']['slope_grand'] * X[x] + 1.96*dict_meansvsstd['total']['std'][x] for x in range(len(X))]  )
    ax.scatter(df['Charge'],np.percentile((df['Charge'], df['Slope']),95,axis=0))

    #ax.plot(dict_meansvsstd['total']['means_index'],dict_meansvsstd['total']['95+'])
    #ax.plot(dict_meansvsstd['total']['means_index'], dict_meansvsstd['total']['95-'])

    #sns.set(rc={"lines.linewidth": 0.2})
    #sns.regplot(x=dict_meansvsstd['total']['means_index'], y=dict_meansvsstd['total']['means'],ci=95, scatter=False, color='r',ax=ax)
    #sns.regplot(x=dict_meansvsstd['total']['means_index'], y=dict_meansvsstd['total']['3std-'],scatter=False, color='r',ax=ax)
    ax.scatter(dict_meansvsstd['total']['means_index'], dict_meansvsstd['total']['means'], color='g', marker='+')


    ax.legend()
    plt.show()


def test2():
    fig = plt.figure()

    ax = fig.add_subplot()

    x = df['Charge']
    y = df['Slope']

    z = np.polyfit(x,y,1)
    p = np.poly1d(z)
    fit = p(x)

    c_y = [np.min(fit), np.max(fit)]
    c_x = [np.min(x), np.max(x)]

    p_y = z[0] * x + z[1]

    # calculate the y-error (residuals)
    y_err = y - p_y

    # create series of new test x-values to predict for
    p_x = np.arange(np.min(x), np.max(x) + 1, 1)

    # now calculate confidence intervals for new test x-series
    mean_x = np.mean(x)  # mean of x
    n = len(x)  # number of samples in origional fit
    t = 2.31  # appropriate t value (where n=9, two tailed 95%)
    s_err = np.sum(np.power(y_err, 2))  # sum of the squares of the residuals

    confs = t * np.sqrt((s_err / (n - 2)) * (1.0 / n + (np.power((p_x - mean_x), 2) /
                                                        ((np.sum(np.power(x, 2))) - n * (np.power(mean_x, 2))))))
    print(confs)
    # now predict y based on test x-values
    p_y = z[0] * p_x + z[1]

    # get lower and upper confidence limits based on predicted y and confidence intervals
    lower = p_y - abs(confs)
    upper = p_y + abs(confs)

    ax.scatter(x,y)
    ax.plot(c_x,c_y)
    ax.plot(p_x,lower)
    ax.plot(p_x, upper)
    plt.show()


def test3():
    fig = plt.figure()

    ax = fig.add_subplot()

    x = df['Charge']
    y = df['Slope']

    res =stats.theilslopes(y,x,0.9)
    lsq_res = stats.linregress(x, y)

def test():
    fig = plt.figure()

    ax = fig.add_subplot()

    x = df['Charge']
    y = df['Slope']

    f = lambda x, *p: polyval(p, x)
    p, cov = curve_fit(f, x, y, [1, 1])
    xi = linspace(np.min(x), np.max(x), 100)

    ps = np.random.multivariate_normal(p, cov, 10000)
    ysample = np.asarray([f(xi, *pi) for pi in ps])
    lower = percentile(ysample, 5, axis=0)
    upper = percentile(ysample, 95, axis=0)

    y_fit = poly1d(p)(xi)
    print(upper - y_fit)

    ax.plot(x, y, 'bo')
    ax.plot(xi, y_fit, 'r-')
    ax.plot(xi, lower, 'b--')
    ax.plot(xi, upper, 'b--')

    plt.show()
def Grouped_charge_histogram(df):
    # Let's make a plot
    fig = plt.figure()
    ax = fig.add_subplot()

    df_bycharge = df.groupby(["Charge"])


    for charge,slope in df_bycharge:


        hist = ax.hist(slope['Slope'], bins=100, range=None, histtype='step', align='mid', orientation='vertical', rwidth=None, label="+{0}".format(slope['Charge'].max()),
            log=False, stacked=False)

        max_y = max(hist[0])
        for j in range(len(hist[0])):
            if hist[0][j] == max_y:
                max_x = hist[1][j]




        ax.text(max_x,max_y,"+{0}".format(slope['Charge'].max()))


    ax.set_title("Charge Grouped Observations Histogram")
    ax.set_xlabel("Slope")
    ax.set_ylabel("Frequency")
    #ax.legend()

    plt.show()

if __name__ == '__main__':

    df= CreateDF(r"C:\Data\BOB\othercalibrationfolders\UHMRchargestatecalibration")
    dict_stats = linear_regressions(df)
    dict_meansvsstd = GetmeansStd(df)
    #twoDhist(df)
    Grouped_charge_histogram(df)
    #test()

