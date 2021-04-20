import os, sys
import csv

import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import webbrowser
from scipy import stats
import seaborn as sns
import numpy as np
from matplotlib.lines import Line2D


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
                    del calibration['Slope']
                    del calibration['Noise']
                    del calibration['Baseline']
                    del calibration['Resolution']
                    del calibration['File']

                    dict_temp = dict(calibration)
                    dict_temp['standard'] = standard
                    calibration_list.append(dict_temp)



    df =  pd.DataFrame(calibration_list)

    df[['Charge', 'm/z','Intensity']] = df[['Charge', 'm/z','Intensity']].apply(pd.to_numeric, errors='coerce')

    return df

if __name__ == '__main__':
    basepath = "C:/Data/BOB/othercalibrationfolders/UHMRchargestatecalibration"

    df = CreateDF(basepath)
    df_bycharge = df.groupby(["Charge"])
    n=30
    color=iter(cm.rainbow(np.linspace(0,1,n)))

    for charge, slope in df_bycharge:
        spectra = [(row["m/z"],1) for index, row in slope.iterrows() ]
        ordered_spectra = sorted(spectra)
        ch = slope['Charge'].max()
        min = ordered_spectra[0][0]
        max= ordered_spectra[-1][0]

        bins = np.arange(min,max,0.2)

        da_tol = 0.00/ch

        index = 0
        bin_index = 0
        gaps = []
        while index < len(ordered_spectra)-1:
            gap = ordered_spectra[index+1][0] - ordered_spectra[index][0]

            if gap <= da_tol:
                #combine peaks
                newpeak1 = (ordered_spectra[index][0],ordered_spectra[index][1] + ordered_spectra[index+1][1])
                newpeak2 = (ordered_spectra[index+1][0],ordered_spectra[index][1] + ordered_spectra[index+1][1])
                ordered_spectra.pop(index+1)
                ordered_spectra.pop(index)
                try:
                    if ordered_spectra[index][1] > ordered_spectra[index+1][1]:

                        ordered_spectra.insert(index,newpeak1)
                    else:
                        ordered_spectra.insert(index,newpeak2)
                except:
                    index += 1
                    continue

            else:
                index += 1

        """    hist = plt.hist(gaps, bins=1000, range=None, align='mid', orientation='vertical',
                       rwidth=None,
                       log=False, stacked=False)"""
        c=next(color)
        plt.vlines([ordered_spectra[x][0] for x in range(len(ordered_spectra))],0,[ordered_spectra[x][1] for x in range(len(ordered_spectra))],colors=c, label=ch)
    plt.legend()
    plt.show()