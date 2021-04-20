import os, sys
import csv

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import webbrowser
from scipy import stats
import seaborn as sns
import numpy as np
from matplotlib.lines import Line2D
from i2ms_KDE import *

def CreateDF(basepath):

    calibration_list = []
    dirs = os.listdir(basepath)

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

                    del calibration['Noise']
                    del calibration['Baseline']
                    del calibration['Resolution']
                    del calibration['File']

                    dict_temp = dict(calibration)
                    dict_temp['standard'] = standard
                    calibration_list.append(dict_temp)



    df =  pd.DataFrame(calibration_list)

    df[['Charge', 'm/z','Slope']] = df[['Charge', 'm/z','Slope']].apply(pd.to_numeric, errors='coerce')

    return df

if __name__ == '__main__':

    import time
    basepath = r"C:\Data\BOB\othercalibrationfolders\test"

    df = CreateDF(basepath)
    raw_data = []
    charge = 67
    peak_iterations = 100
    peak = 697.561


    res = (140000,200)


    start_range = peak-(0.5/charge)
    end_range = peak+(0.5/charge)
    sampling_rate = 0.01/charge
    mid_mz = peak



    current_res = res[0] * np.sqrt(res[1]/mid_mz)


    theo_width = (mid_mz/current_res)/2
    prediction_bandwith5 = (((mid_mz/current_res)*0.1963)-0.0002)

    print("{0},{1},{2}".format(mid_mz,current_res,prediction_bandwith5))

    data = []
    pos = []
    #bandwidth_norm_kde = prediction_bandwith5

    for i in range(-peak_iterations,peak_iterations):
        try:
            mzrange = (start_range+(1.0030*i)/charge,end_range+(1.0030*i)/charge)
            bins = 25

            #spectra = [row["m/z"] for index, row in df.iterrows() if mzrange[0] < row["m/z"] < mzrange[1] ]
            slopes = [row["Slope"] for index, row in df.iterrows() if mzrange[0] < row["m/z"] < mzrange[1] ]

            #ordered_spectra = sorted(spectra)
            #raw_data += ordered_spectra
            from scipy.stats import norm
            from scipy.stats import cauchy
            from scipy.stats import gamma
            import matplotlib.mlab as mlab

            #x_d = np.arange(mzrange[0], mzrange[1], sampling_rate)

            #(mu, sigma) = norm.fit(ordered_spectra)
            data.append(slopes)
            pos.append(np.mean(mzrange))


            #n, bins, patches =plt.hist(ordered_spectra,bins=bins, histtype='step',label="histogram bins = {0}".format(bins))
            #hist_max = max(n)
            plt.xlabel("m/z")
            #min_ = ordered_spectra[0]
            #max_ = ordered_spectra[-1]

            """y3 = []
            x3 = []

            start_index = 0
            for i in range(len(x_d)):

                y_add, start_index = (SinC_multi(x_d[i], ordered_spectra, SinC_bandwidth, start_index))
                if y_add > 0:
                    y3.append(y_add)
                    x3.append(x_d[i])


            max2 = max(y3)
            y3 = [y3[x]*hist_max/max2 for x in range(len(y3))]
            plt.plot(x3, y3,"blue")

            y4 = []
            x4 = []

            start_index = 0
            for i in range(len(x_d)):

                y_add, start_index = (Lorentzian_multi3(x_d[i], ordered_spectra, bandwidth_L_kde,start_index))
                if y_add > 0:
                    y4.append(y_add)
                    x4.append(x_d[i])

            max2 = max(y4)
            y4 = [y4[x] * hist_max / max2 for x in range(len(y4))]
            plt.plot(x4, y4, "black")"""


            """            y5 = []
            x5 = []

            start_index = 0
            for i in range(len(x_d)):

                y_add, start_index = (norm_KDE(x_d[i], ordered_spectra, res,start_index))
                if y_add > 0:
                    y5.append(y_add)
                    x5.append(x_d[i])

            max2 = max(y5)"""
            #y5 = [y5[x] * hist_max / max2 for x in range(len(y5))]
            #plt.plot(x5, y5, "black")
            #plt.vlines(np.mean(mzrange),(hist_max/2)-hist_max*0.05,(hist_max/2)+hist_max*0.05)
            #plt.vlines([np.mean(mzrange)+theo_width,np.mean(mzrange)-theo_width], (hist_max/2)-hist_max*0.025,(hist_max/2)+hist_max*0.025,color="r",alpha=0.25)
            #plt.hlines((hist_max/2),np.mean(mzrange)-theo_width,np.mean(mzrange)+theo_width,color="r",alpha=0.25)



        except:
            ""





    #plt.scatter(raw_data,[-0.2]*len(raw_data),marker="+",color="blue",label="raw data")
    #customlegend = [Line2D([0], [0], color='black', lw=4),Line2D([0], [0], color='blue', lw=4),Line2D([0], [0], color='gray', lw=4),Line2D([0], [0], color='blue', marker="+", lw=0,label='Scatter')]
    #plt.legend(customlegend, ["Lorentzian KDE bandwidth={0}".format(bandwidth_L_kde),"Sinc KDE bandwidth={0}".format(SinC_bandwidth), "norm KDE bandwidth={0}".format(bandwidth_norm_kde),"Raw Data"])
    #customlegend = [Line2D([0], [0], color='black', lw=4),Line2D([0], [0], color='blue', marker="+", lw=0,label='Scatter')]
    #plt.legend(customlegend, ["norm KDE bandwidth={0}".format(bandwidth_norm_kde),"Raw Data"])

    plt.violinplot(data)

    plt.show()
