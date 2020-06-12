# Richard LeDuc, 20200522
# Parses I2MS calibration data and summarizes the results

# Imports
import os, sys
import csv
from appJar import gui
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import webbrowser

import numpy as np

from scipy import stats
import seaborn as sns


# Globals
#  data set:
basepath = r'C:\Data\BOB\othercalibrationfolders\UHMRchargestatecalibration'



def CreateDF(basepath):
    global df
    global progress
    calibration_list = []

    # Actual Work
    dirs = os.listdir(basepath)
    total_count = len(dirs)
    count = 0
    for entry in dirs:
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
                    calibration_list.append(dict(calibration))
        count += 1
        progress = round((count/total_count)*50)
        app.registerEvent(updatprogress)

    df =  pd.DataFrame(calibration_list)
    df[['Charge', 'Slope']] = df[['Charge', 'Slope']].apply(pd.to_numeric, errors='coerce')
    return df

def nice_histogram(df):
    # Let's make a plot
    fig.clf()
    ax = fig.add_subplot()


    # Nice histograph across all charges
    ax.hist(df['Slope'], bins=1000, range=None, density=True, histtype='step', align='mid', orientation='vertical', rwidth=None,
                          log=False, color=None, label=None, stacked=False)
    app.queueFunction(app.refreshPlot, "p1")
# bins='rice' # Another binning option

#exit()

#---------------------------------------
# Histogram by charge state
#c1 = df.loc[df['Charge']==9, 'Slope']
#c2 = df.loc[df['Charge']==10, 'Slope']
#c3 = df.loc[df['Charge']==11, 'Slope']

#c4 = df.loc[df['Charge']==67, 'Slope']
#c5 = df.loc[df['Charge']==68, 'Slope']
#c6 = df.loc[df['Charge']==69, 'Slope']
#c7 = df.loc[df['Charge']==70, 'Slope']
#c8 = df.loc[df['Charge']==71, 'Slope']
#c9 = df.loc[df['Charge']==72, 'Slope']
#c10 = df.loc[df['Charge']==73, 'Slope']
#c11 = df.loc[df['Charge']==74, 'Slope']
#c12 = df.loc[df['Charge']==75, 'Slope']
#c13 = df.loc[df['Charge']==32, 'Slope']
#c14 = df.loc[df['Charge']==33, 'Slope']
#c15 = df.loc[df['Charge']==34, 'Slope']
#c16 = df.loc[df['Charge']==35, 'Slope']
#c17 = df.loc[df['Charge']==36, 'Slope']
#c18 = df.loc[df['Charge']==37, 'Slope']


#kwargs = dict(alpha=0.5, bins=400)

##plt.hist(c1, **kwargs, color='g', label='9')
##plt.hist(c2, **kwargs, color='b', label='10')
##plt.hist(c3, **kwargs, color='r', label='11')

#plt.hist(c1, **kwargs, label='9')
#plt.hist(c2, **kwargs, label='10')
#plt.hist(c3, **kwargs, label='11')
#plt.hist(c4, **kwargs, label='67')
#plt.hist(c5, **kwargs, label='68')
#plt.hist(c6, **kwargs, label='69')
#plt.hist(c7, **kwargs, label='70')
#plt.hist(c8, **kwargs, label='71')
#plt.hist(c9, **kwargs, label='72')
#plt.hist(c10, **kwargs, label='73')
#plt.hist(c11, **kwargs, label='74')
#plt.hist(c12, **kwargs, label='75')
#plt.hist(c13, **kwargs, label='32')
#plt.hist(c14, **kwargs, label='33')
#plt.hist(c15, **kwargs, label='34')
#plt.hist(c16, **kwargs, label='35')
#plt.hist(c17, **kwargs, label='36')
#plt.hist(c18, **kwargs, label='37')


#plt.gca().set(title='Frequency Histogram of Slope by Charge', ylabel='Frequency')
###plt.xlim(50,75)
##plt.legend();

#plt.show()

#---------------------------------------
# 2d Histograpms

def twoDhist(df):

    fig.clf()
    ax = fig.add_subplot()
    ax.hist2d(df['Charge'], df['Slope'], bins=(75, 75), norm=colors.PowerNorm(gamma=1./2.), cmap=plt.cm.Greys)
    #ax.set_colorbar()
    app.queueFunction(app.refreshPlot, "p1")

#---------------------------------------
# Linear regression

def getChargeGroups(df):
    charges = sorted(list(set(df['Charge'])))
    charge_groups = []
    Current_group = [charges[0]]
    for i in range(1, len(charges)):
        if charges[i] - charges[i - 1] > 5:
            charge_groups.append(Current_group)
            Current_group = [charges[i]]
        else:
            Current_group.append(charges[i])
    charge_groups.append(Current_group)

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
    for x in range(len(dict_key_map)-1,len(charge_groups)):
        charge_groups[len(dict_key_map)-1] += charge_groups[x]

    J = len(dict_key_map)
    if len(dict_key_map)> len(charge_groups):
        J = len(charge_groups)


    for j in range(J):
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
            ax.scatter(dict_stats[group]['sub_charge'], dict_stats[group]['sub_slope'], marker='+')
            ax.plot(X, dict_stats[group]['intercept'] + dict_stats[group]['slope']* X, label='{0} charge data'.format(group))

    # Plot subsample of data using seaborn

    #sns.regplot(x=sub_high_charge, y=sub_high_slope, color='g', marker='+')

    ax.set_title('Relationship between Charge and Slope', size=24)
    ax.set_xlabel('Charge', size=18)
    ax.set_ylabel('Slope', size=18);
    ax.legend()
    app.queueFunction(app.refreshPlot, "p1")


def residPlots(dict_stats,key):
    fig.clf()
    ax = fig.add_subplot()
    sns.residplot(x=dict_stats[key]['sub_charge'], y=dict_stats[key]['sub_slope'], ax=ax)
    ax.set_title("{0} Charge Residual".format(key))
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
    for x in range(len(dict_key_map) - 1, len(charge_groups)):
        charge_groups[len(dict_key_map) - 1] += charge_groups[x]

    J = len(dict_key_map)
    if len(dict_key_map) > len(charge_groups):
        J = len(charge_groups)
    for j in range(J):
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
    ax.set_title("Mean Slope vs Std Slope For Each Charge")
    ax.legend()
    app.queueFunction(app.refreshPlot, "p1")


def meanSlopevsCharge3std(dict_meansvsstd):
    fig.clf()
    ax = fig.add_subplot()
    ax.set_title('Three STD from Mean Slope by Charge', size=18)
    ax.set_xlabel('Charge', size=14)
    ax.set_ylabel('Slope', size=14);

    sns.regplot(x=dict_meansvsstd['total']['means_index'], y=dict_meansvsstd['total']['3std+'], color='r',ax=ax)
    sns.regplot(x=dict_meansvsstd['total']['means_index'], y=dict_meansvsstd['total']['3std-'], color='r',ax=ax)

    ax.scatter(dict_meansvsstd['total']['means_index'], dict_meansvsstd['total']['means'], color='g', marker='+')

    app.queueFunction(app.refreshPlot, "p1")
def QQPlot(df,charge=37):
    # We need a plot of the difference between means! Basically, the standard error of the summed observed slopes needs to be
    # less than the difference between means of two charge states. this difference is constant over the linear range, so we
    # only need to take the average difference and reduce our standard error below that.
    fig.clf()
    ax = fig.add_subplot()
    # QQ Plots for specific charge states


    stats.probplot(df.loc[df['Charge']==charge, 'Slope'], plot=ax)
    ax.set_title('QQ Plot Charge = {0}'.format(charge))
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

        if selected[0] == "meanSlopeVsCharge":
            app.thread(meanSlopevsCharge3std,dict_meansvsstd)

        if selected[0] == "MeanVstd":
            app.thread(MeanvsStdCharges, dict_meansvsstd)
        if selected[1] == "Resid":
            app.thread(residPlots, dict_stats,selected[0])
        if selected[0] == "ChargeVsSlope":
            app.thread(Chargevsslope, dict_stats)

        if selected[0] == "nice_hist":
            app.thread(nice_histogram, df)

        if selected[0] == "twoDHistogram":

            app.thread(twoDhist, df)

def openTree(xml):

    app.startFrame("Left", row=0, column=0)

    app.tree("tree1", xml, dbl=treeclick, editable=False)
    app.stopFrame()

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
    progress = 0
    app.registerEvent(updatprogress)
    app.showSubWindow("prog")


    dirname = basepath.split("/")[-1]
    CreateDF(basepath)


    charge_groups = getChargeGroups(df)
    for x in range(len(dict_key_map) - 1, len(charge_groups)):
        charge_groups[len(dict_key_map) - 1] += charge_groups[x]

    if len(dict_key_map) > len(charge_groups):
        dict_key_map = dict_key_map[:len(charge_groups)]


    dict_stats = linear_regressions(df)
    progress = 75
    app.registerEvent(updatprogress)
    dict_meansvsstd = GetmeansStd(df)

    progress = 100
    app.registerEvent(updatprogress)

    charges = sorted(list(set(df['Charge'])))
    QQplots = "<QQPlots>"
    for each in charges:
        QQplots += "<_{0}></_{0}>".format(each)
    QQplots += "</QQPlots>"

    resids = "<ResidPlots>"
    for each in dict_key_map:
        resids += '<{0} id="Resid"></{0}>'.format(each)
    resids += '</ResidPlots>'

    xml = """
            <{0}>
            <Histograms><nice_hist></nice_hist>
            <twoDHistogram></twoDHistogram></Histograms>
            {2}
            <ChargeVsSlope></ChargeVsSlope>
            <MeanVstd></MeanVstd>
            <meanSlopeVsCharge></meanSlopeVsCharge>

            {1}
            </{0}>
            """.format(dirname, QQplots, resids)

    openTree(xml)
    progress = 100
    app.registerEvent(updatprogress)
    app.queueFunction(app.hideSubWindow, "prog")

def file_press(menu):
    if menu =="Open":

        basepath = app.directoryBox("opendir")
        app.thread(load_data,basepath)






def about_press(menu):
    if menu == "Help":
        app.showSubWindow("help")
    if menu == "Source code":
        filename = os.path.join(os.path.dirname(sys.executable), 'SourceCode.txt')
        os.startfile(filename)

    if menu == "Version":
        app.showSubWindow("version")

if __name__ == '__main__':
    progress = 0
    dict_key_map = ['low', 'mid', 'high', 'very_high']
    df = []
    dict_meansvsstd = {}
    dict_stats = {}

    version = "0.1 alpha build"  # add version here
    # Create app - name and size
    app = gui("I2MS Calibration Diagnostics {0}".format(version), "1200x800", showIcon=False)

    file_menus = ["Open"]
    about_menu = ["Version", "Help", "Source code"]
    app.addMenuList("file", file_menus, file_press)
    app.addMenuList("About", about_menu, about_press)

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
    app.startFrame("Right", row=0, column=1)
    fig = app.addPlotFig("p1", showNav=True)
    app.stopFrame()


    #progress window
    app.startSubWindow("prog","Progress")
    app.setTransparency(50)
    app.addMeter("progress")
    app.stopSubWindow()




    app.go()