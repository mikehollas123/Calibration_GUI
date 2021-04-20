# Richard LeDuc, 20200522
# & Michael Hollas 20200612

# Parses I2MS calibration data and summarizes the results

# Imports
import os, sys
import time
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
from matplotlib.lines import Line2D
from scipy.stats import norm


#Globals
axis_text_size = 24
tick_text_size = 20
title_text_size = 24
line_width = 1

class Mikes_CI_and_PI_calc:

    def __init__(self, X, Y):
        import numpy as np
        import scipy
        x = np.array(X)
        y = np.array(Y)
        n = len(x)
        M = x[:, np.newaxis] ** [0, 1]  # Linear regression
        reg, R, a, b = scipy.linalg.lstsq(M, Y)  # least squares fitting
        self.Intercept, self.Slope = reg

        y_prime = x * self.Slope + self.Intercept
        var = y - y_prime
        var_sq = var * var
        s_est = np.sqrt(sum(var_sq) / y.size - 2)

        x_minus_xmean = x - np.mean(x)
        sum_x2 = sum(x_minus_xmean * x_minus_xmean)

        self.Standard_error = s_est / np.sqrt(sum_x2)

        self.t_value = scipy.stats.t.ppf(0.95, n- 1)

        """Confidence Interval for the Slope
        lower limit: slope - (t.95)(sb)
    upper limit: b + (t.95)(sb)
        """
        self.CI_plus = self.Slope + self.t_value * self.Standard_error
        self.CI_minus = self.Slope - self.t_value * self.Standard_error
        x_pre = np.linspace(x.min(), x.max(), n)
        x_pre_xminusmean = x_pre - np.mean(x_pre)
        se_y = s_est * np.sqrt(1 / n + (x_pre_xminusmean * x_pre_xminusmean) / sum(x_pre_xminusmean * x_pre_xminusmean))
        pre_se_y = s_est * np.sqrt(
            1 + 1 / n + (x_pre_xminusmean * x_pre_xminusmean) / sum(x_pre_xminusmean * x_pre_xminusmean))

        self.Regression_XY = (x_pre, x_pre * self.Slope + self.Intercept)
        self.Upper_CI_XY = (x_pre, x_pre * self.Slope + self.Intercept + se_y * self.t_value)
        self.Lower_CI_XY = (x_pre, x_pre * self.Slope + self.Intercept - se_y * self.t_value)
        self.Upper_PI_XY = (x_pre, x_pre * self.Slope + self.Intercept + pre_se_y * self.t_value)
        self.Lower_PI_XY = (x_pre, x_pre * self.Slope + self.Intercept - pre_se_y * self.t_value)


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
                    del calibration['Intensity']
                    del calibration['Noise']
                    del calibration['Baseline']
                    del calibration['Resolution']
                    del calibration['File']

                    dict_temp = dict(calibration)
                    dict_temp['standard'] = standard
                    calibration_list.append(dict_temp)

        count += 1
        progress = round((count/total_count)*50)
        app.registerEvent(updatprogress)

    df =  pd.DataFrame(calibration_list)

    df[['Charge', 'Slope','m/z']] = df[['Charge', 'Slope','m/z']].apply(pd.to_numeric, errors='coerce')

    return df

def Violin_Slope_dist(df,x_lims=0):

    X = dict_stats['total']['charge']
    fig.clf()
    ax = fig.add_subplot()

    if x_lims != 0:
        ax.set_xlim(min(x_lims),max(x_lims))

    plot_data = Mikes_CI_and_PI_calc(df['Charge'], df['Slope'])

    df_bycharge = df.groupby(["Charge"])
    data = []
    pos = []
    for charge, slope in df_bycharge:
        ch = slope['Charge'].max()
        data.append(np.array(slope["Slope"]/1000000))
        pos.append(ch)
        mean = np.mean(slope["Slope"]/1000000)
        quartile1, medians, quartile3 = np.percentile(np.array(slope["Slope"]/1000000), [25, 50, 75], axis=0)
        ax.vlines(ch,quartile1,quartile3,linewidth=3)
        ax.scatter(ch, mean, marker='o', color='black',s=50)

    violin_parts = ax.violinplot(data, pos, showmeans=False, showextrema=False, showmedians=False)
    ax.set_title("Charge Grouped Observations Violin Plot",fontsize=title_text_size)
    ax.set_xlabel("Charge",fontsize=axis_text_size)
    ax.set_ylabel("Slope (x$10^6$)",fontsize=axis_text_size)
    ax.tick_params(axis='y',labelsize= tick_text_size)
    ax.tick_params(axis='x',  labelsize=tick_text_size)
    for pc in violin_parts['bodies']:
        pc.set_color('red')
        pc.set_edgecolor('red')
        pc.set_alpha(0.4)

    print(plot_data.Slope,plot_data.Intercept)
    ax.plot(plot_data.Regression_XY[0],plot_data.Regression_XY[1]/1000000,color='black',linewidth=2)
    ax.plot(plot_data.Lower_CI_XY[0],plot_data.Lower_CI_XY[1]/1000000,color='lightsteelblue')

    ax.plot(plot_data.Upper_CI_XY[0],plot_data.Upper_CI_XY[1]/1000000,color='lightsteelblue')
    ax.plot(plot_data.Lower_PI_XY[0], plot_data.Lower_PI_XY[1]/1000000,color='plum')

    ax.plot(plot_data.Upper_PI_XY[0], plot_data.Upper_PI_XY[1]/1000000,color='plum')

    customlegend = [Line2D([0], [0], color='black', lw=4), Line2D([0], [0], color='plum', lw=4),Line2D([0], [0], color='lightsteelblue', lw=4),Line2D([0], [0], marker='o', color='black', label='Scatter',
                          markerfacecolor='black', markersize=10)]
    ax.legend(customlegend, ["Regression", "Prediction Interval","Confidence Interval","Mean"])



    app.queueFunction(app.refreshPlot, "p1")

def Grouped_charge_histogram(df):
    # Let's make a plot
    fig.clf()
    ax = fig.add_subplot()

    df_bycharge = df.groupby(["Charge"])


    for charge,slope in df_bycharge:

        ch = slope['Charge'].max()
        hist = ax.hist(slope['Slope'], bins=100, range=None, histtype='step', align='mid', orientation='vertical', rwidth=None, label="+{0}".format(ch),
            log=False, stacked=False)


        max_y = max(hist[0])
        for j in range(len(hist[0])):
            if hist[0][j] == max_y:
                max_x = hist[1][j]

        slope_from_Charge = dict_stats['total']['intercept_grand'] + dict_stats['total']['slope_grand'] * ch
        ax.vlines(slope_from_Charge,0,max_y,'r')

        ax.vlines(np.mean(slope['Slope']),0,max_y )
        ax.text(max_x,max_y,"+{0}".format(slope['Charge'].max()))


    ax.set_title("Charge Grouped Observations Histogram",fontsize=title_text_size)
    ax.set_xlabel("Slope",fontsize=axis_text_size)
    ax.set_ylabel("Frequency",fontsize=axis_text_size)
    ax.tick_params(axis='y',labelsize= tick_text_size)
    ax.tick_params(axis='x',  labelsize=tick_text_size)
    customlegend = [Line2D([0],[0],color='black',lw=4),Line2D([0],[0],color='red',lw=4)]
    ax.legend(customlegend,["Mean","From Regression"])

    app.queueFunction(app.refreshPlot, "p1")



def Grouped_charge_histogram_norm(df):
    # Let's make a plot
    fig.clf()
    ax = fig.add_subplot()

    df_bycharge = df.groupby(["Charge"])


    for charge,slope in df_bycharge:


        hist = ax.hist(slope['Slope'], bins=100, range=None, density=True, histtype='step', align='mid', orientation='vertical', rwidth=None, label="+{0}".format(slope['Charge'].max()),
            log=False, stacked=True)

        max_y = max(hist[0])
        for j in range(len(hist[0])):
            if hist[0][j] == max_y:
                max_x = hist[1][j]



        ax.text(max_x,max_y,"+{0}".format(slope['Charge'].max()))


    ax.set_title("Charge Grouped Observations Histogram \n- Normalized by Sum",fontsize=title_text_size)
    ax.set_xlabel("Slope",fontsize=axis_text_size)
    ax.set_ylabel("Frequency",fontsize=axis_text_size)
    ax.tick_params(axis='y',labelsize= tick_text_size)
    ax.tick_params(axis='x',  labelsize=tick_text_size)
    #ax.legend()


    app.queueFunction(app.refreshPlot, "p1")

def nice_histogram(df):
    # Let's make a plot
    fig.clf()
    ax = fig.add_subplot(label="1")
    ax2 = fig.add_subplot(label="2", frame_on=False)

    ax2.xaxis.tick_bottom()
    ax2.yaxis.tick_right()
    ax2.xaxis.set_label_position('bottom')
    ax2.yaxis.set_label_position('right')
    ax2.set_xlabel("Charge", color="r")
    ax2.set_ylabel("Frequency", color="r")
    # Nice histograph across all charges

    stats = linear_regressions(df)
    grand_slope = stats['total']['slope_grand']
    grand_intercept = stats['total']['intercept_grand']

    X = ax.hist(df['Slope'], bins=1000, range=None, histtype='step', align='mid', orientation='vertical', rwidth=None,
                          log=False, color='r', label="Slope", stacked=False)
    start = (X[1][0]-grand_intercept)/grand_slope
    end = (X[1][-1]-grand_intercept)/grand_slope

    print(start)
    print(end)

    Y = ax2.hist(df['Charge'], bins=1000, range=(start,end), histtype='step', align='mid', orientation='vertical', rwidth=None,
            log=False, color=None, label="Charge", stacked=False)

    start = Y[1]
    end = Y[0]

    print(start)

    ax.set_title("Observation Histogram",fontsize=title_text_size)
    ax.set_xlabel("Slope",fontsize=axis_text_size)
    ax.set_ylabel("Frequency",fontsize=axis_text_size)

    ax.tick_params(axis='y', labelsize=tick_text_size)
    ax.tick_params(axis='x', labelsize=tick_text_size)
    ax.legend()
    ax2.legend()
    app.queueFunction(app.refreshPlot, "p1")
# bins='rice' # Another binning option

# 2d Histograpms

def twoDhist(df):

    fig.clf()
    cmap = plt.cm.Reds
    norm = colors.PowerNorm(gamma=0.5)

    ax = fig.add_subplot()
    ax.hist2d(df['Charge'], df['Slope'], bins=(1000, 1000), norm=norm, cmap=cmap)
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
    X = df['Charge']
    ax.set_title(dirname,fontsize=title_text_size)
    ax.plot(X, dict_stats['total']['intercept_grand'] + dict_stats['total']['slope_grand'] * X, 'black', linewidth=0.2)


    #sns.set(rc={"lines.linewidth": 0.2})
    sns.regplot(x=dict_meansvsstd['total']['means_index'], y=dict_meansvsstd['total']['3std+'],scatter=False, color='r',ax=ax)
    sns.regplot(x=dict_meansvsstd['total']['means_index'], y=dict_meansvsstd['total']['3std-'],scatter=False, color='r',ax=ax)

    ax.set_xlabel("Charge", size = axis_text_size)
    ax.set_ylabel("Slope", size=axis_text_size)
    ax.tick_params(axis='y', labelsize=tick_text_size)
    ax.tick_params(axis='x', labelsize=tick_text_size)

    app.queueFunction(app.refreshPlot, "p1")
    #app.queueFunction(sns.set,rc={"lines.linewidth": line_width})
#---------------------------------------
# Linear regression

def getChargeGroups(df):
    global dict_key_map
    dict_key_map = list(set(df['standard']))
    charge_groups = []
    for each in dict_key_map:

        charge_groups.append( list(set(df.loc[df['standard'] == each,'Charge'])))



    return charge_groups

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

def Chargevsslope(dict_stats):
    X = dict_stats['total']['charge']
    Y = dict_stats['total']['slope']
    fig.clf()
    ax = fig.add_subplot()

    ax.plot(X, dict_stats['total']['intercept_grand'] + dict_stats['total']['slope_grand'] * X, 'r', label='Total data')
    for group in dict_stats:
        if group != 'total':
            #ax.scatter(dict_stats[group]['sub_charge'], dict_stats[group]['sub_slope'], marker='+')
            ax.plot(X, dict_stats[group]['intercept'] + dict_stats[group]['slope']* X, label='{0}'.format(group))

    # Plot subsample of data using seaborn

    #sns.regplot(x=sub_high_charge, y=sub_high_slope, color='g', marker='+')

    ax.set_title('Relationship between Charge and Slope', size=title_text_size)
    ax.set_xlabel('Charge', size=axis_text_size)
    ax.set_ylabel('Slope', size=axis_text_size)

    ax.tick_params(axis='y', labelsize=tick_text_size)
    ax.tick_params(axis='x', labelsize=tick_text_size)
    ax.legend()
    app.queueFunction(app.refreshPlot, "p1")


def residPlots(dict_stats,key):
    fig.clf()
    ax = fig.add_subplot()
    sns.residplot(x=dict_stats[key]['sub_charge'], y=dict_stats[key]['sub_slope'], ax=ax)
    ax.set_title("{0} Charge Residual".format(key),size=title_text_size)
    ax.tick_params(axis='y', labelsize=tick_text_size)
    ax.tick_params(axis='x', labelsize=tick_text_size)
    ax.set_xlabel('Charge', size=axis_text_size)
    ax.set_ylabel('Slope Residual', size=axis_text_size)
    app.queueFunction(app.refreshPlot, "p1")

def GetmeansStd(df):
    dict_meansvsstd = {}

    df_total_std = df
    df_total_std = df_total_std.groupby(["Charge"])['Slope'].std()
    df_total_means = df
    df_total_means = df_total_means.groupby(["Charge"])['Slope'].mean()
    dict_meansvsstd['total'] = {
        'means':df_total_means.values,
        'means_index':df_total_means.index,
        'std': df_total_std.values,
        'std_index': df_total_std.index,
        '3std+':df_total_means.values + df_total_std.values*3,
        '3std-': df_total_means.values - df_total_std.values * 3
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


def MeanvsStdCharges(dict_meansvsstd):

    """
    # Compare means and their std
    # plt.hist(df_means.values, bins=(75,75), label='Mean')
    X = df_means.values
    Y = df_std.values
    """
    fig.clf()
    ax = fig.add_subplot()
    for key in dict_meansvsstd:
        if key != 'total':

            ax.scatter(dict_meansvsstd[key]['means'], dict_meansvsstd[key]['std'], label='{0} Charge'.format(key),  marker='.')
    ax.set_title("Mean Slope vs Std Slope For Each Charge",size=title_text_size)
    ax.set_xlabel("Mean Slope", size = axis_text_size)
    ax.set_ylabel("Mean Standard Deviation", size=axis_text_size)
    ax.tick_params(axis='y', labelsize=tick_text_size)
    ax.tick_params(axis='x', labelsize=tick_text_size)
    ax.legend()
    app.queueFunction(app.refreshPlot, "p1")


def meanSlopevsCharge3std(dict_meansvsstd):
    fig.clf()
    ax = fig.add_subplot()
    ax.set_title('Three STD from Mean Slope by Charge', size=title_text_size)
    ax.set_xlabel('Charge', size=axis_text_size)
    ax.set_ylabel('Slope', size=axis_text_size);

    sns.regplot(x=dict_meansvsstd['total']['means_index'], y=dict_meansvsstd['total']['3std+'], color='r',ax=ax)
    sns.regplot(x=dict_meansvsstd['total']['means_index'], y=dict_meansvsstd['total']['3std-'], color='r',ax=ax)

    ax.scatter(dict_meansvsstd['total']['means_index'], dict_meansvsstd['total']['means'], color='g', marker='+')
    ax.tick_params(axis='y', labelsize=tick_text_size)
    ax.tick_params(axis='x', labelsize=tick_text_size)

    app.queueFunction(app.refreshPlot, "p1")
def QQPlot(df,charge=37):



    # We need a plot of the difference between means! Basically, the standard error of the summed observed slopes needs to be
    # less than the difference between means of two charge states. this difference is constant over the linear range, so we
    # only need to take the average difference and reduce our standard error below that.




    fig.clf()
    ax = fig.add_subplot()
    # QQ Plots for specific charge states
    data = df.loc[df['Charge']==charge, 'Slope']

    stats.probplot(data, plot=ax)

    shaperoTest = stats.shapiro(data)
    ks_statistic, p_value = stats.kstest(data, 'norm')

    ax.set_title('QQ Plot Charge: {0}'.format(charge), size=title_text_size)
    ax.set_xlabel("Theoretical Quantiles",size = axis_text_size)
    ax.set_ylabel("Ordered Values", size=axis_text_size)
    ax.tick_params(axis='y', labelsize=tick_text_size)
    ax.tick_params(axis='x', labelsize=tick_text_size)
    ax.text(0,max(data) - (max(data)- min(data)),"Total number of points: {0} \n Shapiro-Wilk test p Value {1}\n Kolmogorov Smirnov test p Value {2}".format(data.count(),shaperoTest.pvalue,p_value))
    app.queueFunction(app.refreshPlot, "p1")

    # In general across most of the range the data is normally distributed. Some charge states break down more than others.

def treeclick(title,atr):
    selected = app.getTreeSelected("tree1")
    if selected != None:

        try:
            QQplot_charge = int(selected[0].split("_")[-1])
            app.thread(QQPlot,df,QQplot_charge)

        except:
            ""
        if selected[0] == 'test':
            print("ha")
        if selected[0] == "meanSlopeVsCharge":
            app.thread(meanSlopevsCharge3std,dict_meansvsstd)

        if selected[0] == "Mean_Vs_std":
            app.thread(MeanvsStdCharges, dict_meansvsstd)
        if selected[1] == "Resid":
            app.thread(residPlots, dict_stats,selected[0])
        if selected[0] == "ChargeVsSlope":
            app.thread(Chargevsslope, dict_stats)

        if selected[0] == "Grouped_charge_histogram":
            print(selected[0])
            app.thread(Grouped_charge_histogram, df)

        if selected[0] == "Grouped_charge_histogram_norm":
            print(selected[0])
            app.thread(Grouped_charge_histogram_norm, df)

        if selected[0] == "twoDHistogram":

            app.thread(twoDhist, df)

        if selected[0] == "Violin_Slope_dist":
            app.thread(Violin_Slope_dist, df)

def openTree(xml):

    with app.panedFrame("Left",row=0,column=0):

        app.tree("tree1", xml, dbl=treeclick, editable=False)


def email():
    webbrowser.open('mailto:michael.hollas@northwestern.edu', new=1)

def updatprogress():
    app.setMeter("progress",progress)

def load_data(basepath):
    global dict_key_map
    global df
    global dict_stats
    global dict_meansvsstd
    global progress
    global dirname
    progress = 0
    app.registerEvent(updatprogress)
    app.showSubWindow("prog")


    dirname = basepath.split("/")[-1]
    CreateDF(basepath)


    """    charge = df['Charge']
    slope = df['Slope']
    with open("C:\Data\BOB\othercalibrationfolders\SlopevsCharge_{0}.csv".format(dirname),"w") as filewrite:
        for i in range(len(charge)):
            filewrite.write("{0},{1}\n".format(slope[i],charge[i]))"""

    charge_groups = getChargeGroups(df)



    dict_stats = linear_regressions(df)
    progress = 75
    app.registerEvent(updatprogress)
    dict_meansvsstd = GetmeansStd(df)

    progress = 100
    app.registerEvent(updatprogress)

    QQplots = "<QQPlots>"
    for k in range(len(charge_groups)):

        charges = sorted(charge_groups[k])
        QQplots += "<{0}>".format(dict_key_map[k])
        for l in range(len(charges)):
            QQplots += "<_{0}></_{0}>".format(charges[l])
        QQplots += "</{0}>".format(dict_key_map[k])
    QQplots += "</QQPlots>"


    resids = "<ResidPlots>"
    for each in dict_key_map:
        resids += '<{0} id="Resid"></{0}>'.format(each)
    resids += '</ResidPlots>'

    xml = """
            <{0}>
            <Histograms><Grouped_charge_histogram></Grouped_charge_histogram>
            <Grouped_charge_histogram_norm></Grouped_charge_histogram_norm><Violin_Slope_dist></Violin_Slope_dist></Histograms>
            {2}
            <ChargeVsSlope></ChargeVsSlope>
            <Mean_Vs_std></Mean_Vs_std>

            {1}
            </{0}>
            """.format(dirname, QQplots, resids)

    openTree(xml)
    progress = 100
    app.registerEvent(updatprogress)
    app.queueFunction(app.hideSubWindow, "prog")

def file_press(menu):
    if menu =="Open":
        basepath = app.directoryBox("dir1")
        #basepath = "C:/Data/BOB/othercalibrationfolders/UHMRchargestatecalibration"
        app.thread(load_data,basepath)






def about_press(menu):
    if menu == "Help":
        app.showSubWindow("help")
    if menu == "Source code":
        filename = os.path.join(os.path.dirname(sys.executable), 'SourceCode.txt')
        os.startfile(filename)

    if menu == "Version":
        app.showSubWindow("version")

def view_press(menu):
    if menu == "Plot Customisation":
        app.showSubWindow("plot_custom")
def applypress(button):

    if button == "Apply":
        global title_text_size
        global axis_text_size
        global tick_text_size
        global line_width
        try:
            title_text_size = int(app.getEntry("Title Font Size"))
            axis_text_size= int(app.getEntry("Axis Font Size"))
            tick_text_size= int(app.getEntry("Tick Font Size"))
            line_width = int(app.getEntry("Line Width"))
            app.thread(treeclick, app.getTreeSelected("tree1"), "")

        except:
            app.errorBox("Input Error","Input must be an integer")


class Click_and_drag():
    def __init__(self,fig, func):

        self.fig = fig
        self.ax = self.fig.axes[0]
        self.func = func
        self.press = False


        self.c1 = self.ax.figure.canvas.mpl_connect('button_press_event', self.onpress)
        self.c2 = self.ax.figure.canvas.mpl_connect('button_release_event', self.onrelease)
        self.selected_min_x = 0
        self.selected_max_x = 0
        self.selected_min_y = 0
        self.selected_max_y = 0

        print("connecting...")

    def onclick(self, event):

        if self.selected_min_x > 0 and self.selected_max_x > 0:
            self.ax = self.fig.axes[0]
            x_lims = sorted([self.selected_max_x,self.selected_min_x])

            X = self.ax.lines[0].get_xdata()
            Y = self.ax.lines[0].get_ydata()

            index = 0
            start_index = 0
            end_index = 0
            while index < len(X):
                if x_lims[0] < X[index] and start_index == 0:
                    start_index = index
                    print(start_index)
                elif x_lims[1] < X[index] and end_index == 0:
                    end_index = index
                    print(end_index)
                index += 1
            print(Y[start_index:end_index])
            y_min = min(Y[start_index:end_index])
            y_max = max(Y[start_index:end_index])
            self.ax.set_ylim(y_min,y_max)
            self.ax.set_xlim(x_lims[0], x_lims[1])
            self.ax.figure.canvas.draw_idle()

    def onpress(self, event):


        if event.button == 1:

            self.click_time = time.clock()
            self.press = True
            if event.xdata != None:
                self.selected_min_x = event.xdata

        else:
            self.ax.autoscale()

            self.ax.figure.canvas.draw_idle()

    def onrelease(self, event):

            if self.press == True and time.clock()-self.click_time > 1:
                self.onclick(event)
                self.press = False
                if event.xdata != None:
                    self.selected_max_x = event.xdata


def func(event):
    print(event.xdata, event.ydata)



if __name__ == '__main__':
    progress = 0
    dict_key_map = []
    df = []
    dict_meansvsstd = {}
    dict_stats = {}

    version = "0.1 alpha build"  # add version here
    # Create app - name and size
    app = gui("I2MS Calibration Diagnostics {0}".format(version), "1200x800",showIcon=False)

    file_menus = ["Open"]
    about_menu = ["Version", "Help", "Source code"]

    view_menu = ["Plot Customisation"]
    app.addMenuList("file", file_menus, file_press)
    app.addMenuList("About", about_menu, about_press)
    app.addMenuList("View", view_menu, view_press)

    # Version Window
    app.startSubWindow("version", "Version")
    app.setSize(300, 200)
    app.label("I2MS Viewer {0}".format(version))
    app.label("Created by Mike Hollas")
    app.stopSubWindow()

    # help window
    app.startSubWindow("help", "Help")
    app.setSize(300, 200)
    app.startFrame("1")
    app.message("\nFor any help please email michael.hollas@northwestern.edu\n")
    app.addButton("email", email)
    app.stopFrame()
    app.stopSubWindow()


    # Create Empty matplotlib plot
    with app.panedFrame("Right",row=0,column=1):

        fig = app.addPlotFig("p1", showNav=True)
        ax = fig.add_subplot()
        #click = Click_and_drag(fig,func)








    #progress window
    app.startSubWindow("prog","Progress")
    app.setTransparency(50)
    app.addMeter("progress")
    app.stopSubWindow()

    #Plot Customization
    app.startSubWindow("plot_custom","Plot Customisation")

    app.addLabelEntry("Title Font Size")
    app.addLabelEntry("Axis Font Size")
    app.addLabelEntry("Tick Font Size")
    app.addLabelEntry("Line Width")
    app.setEntry("Title Font Size",title_text_size)
    app.setEntry("Axis Font Size", axis_text_size)
    app.setEntry("Tick Font Size", tick_text_size)
    app.setEntry("Line Width", line_width)



    app.addButton("Apply",applypress)


    app.stopSubWindow()


    app.go()

