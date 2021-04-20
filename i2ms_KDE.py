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


def norm_KDE(x,mids,resolution,start_index,use_predictor=True):

    y = 0
    mid_idex = start_index
    start_index_new = start_index


    while mid_idex < len(mids):

        if mids[mid_idex] > x + 0.2:
            break
        if mids[mid_idex] < x -0.2:
            mid_idex += 1
            start_index_new = mid_idex-1
        else:
            if use_predictor == True:
                current_res = resolution[0] * np.sqrt(resolution[1]/ mids[mid_idex])
                prediction_bandwith5 = (((mids[mid_idex] / current_res) * 0.1963) - 0.0002)
            else:
                prediction_bandwith5 = resolution
            a = (1 / (prediction_bandwith5 * np.sqrt(2 * np.pi)))
            b = ((x - mids[mid_idex])/prediction_bandwith5*np.sqrt(2))
            y += a*np.exp(-0.5*b*b)
            mid_idex += 1

    return (y,start_index_new )

def Lorentzian_multi3(x,mids,width,start_index):
    y = 0
    mid_idex = start_index
    start_index_new = start_index
    width2 = width * width

    while mid_idex < len(mids):

        if mids[mid_idex] > x + 2:
            break
        if mids[mid_idex] < x -2:
            mid_idex += 1
            start_index_new = mid_idex-1
        else:
            y += (1 / np.pi) * ((width2) / (((x - mids[mid_idex]) * (x - mids[mid_idex])) + (width2)))
            mid_idex += 1

    return (y,start_index_new )

def SinC_multi(x,mids,width,start_index):
    y = 0
    mid_idex = start_index
    start_index_new = start_index
    width2 = width * width

    while mid_idex < len(mids):

        if mids[mid_idex] > x + 2:
            break
        if mids[mid_idex] < x -2:
            mid_idex += 1
            start_index_new = mid_idex-1
        else:

            y += sinc(x,mids[mid_idex],width)
            mid_idex += 1

    return (y,start_index_new )

def sinc(x,mid,width):
    b = (mid - x)/width*5.5
    if b == 0:
        return 1.0
    return abs((np.sin(b)/b))


def Lorentzian_pdf(x,mid,width,weight):

    return weight * (1 / np.pi) * ((width * width) / (((x - mid) * (x - mid)) + (width * width)))
if __name__ == '__main__':

    import time
    basepath = r"C:\Data\BOB\othercalibrationfolders\test"

    df = CreateDF(basepath)


    spectra = [row["m/z"] for index, row in df.iterrows() ]
    ordered_spectra = sorted(spectra)
    print("data loaded")
    start_t = time.time()
    min = ordered_spectra[0]
    max= ordered_spectra[-1]

    x_d = np.arange(min-10, max+10, 0.002)

    y = []
    x = []
    count = 0
    start_index = 0
    for i in range(len(x_d)):

        y_add, start_index = (norm_KDE(x_d[i], ordered_spectra, (140000,200),start_index))
        if y_add > 0:
            y.append(y_add)
            x.append(x_d[i])
        count += 1
        if count % int(len(x_d)/1000) == 0:
            print(float(count/len(x_d))*100)
    print(100)
    print(time.time() - start_t)
    plt.plot(x, y)
    plt.show()



    """   from scipy.stats import norm


    x = [ordered_spectra[x][0] for x in range(len(ordered_spectra))]
    x_d = np.arange(min, x[10000], 0.02)
    density = norm(x[2000], scale=0.025).pdf(x_d)

    count=0"""
    """    for xi in x[1:100000]:


    density += norm(xi, scale=0.025).pdf(x_d)
    count +=1
    if count % 10000:
        print(float(count / 100000) * 100)"""


    """    mz = [ordered_spectra[x][0] for x in range(len(ordered_spectra))]
    weights = [ordered_spectra[x][1] for x in range(len(ordered_spectra))]
    X = np.concatenate((mz,weights))[:, np.newaxis]
    from sklearn.neighbors import KernelDensity

    x_d = np.arange(min, max, 0.02)[:, np.newaxis]
    kde = KernelDensity(bandwidth=0.005, kernel='gaussian')
    kde.fit(X)
    print("fit_done")



    logprob = kde.score_samples(x_d)
    

    plt.vlines(x_d,0,np.exp(logprob))
    plt.show()"""