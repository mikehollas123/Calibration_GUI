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

from mpl_toolkits.mplot3d import Axes3D


def CreateDF(basepath):
    calibration_list = []
    dirs = os.listdir(basepath)

    for entry in dirs:
        try:
            standard = entry.split("_")[0]
        except:
            standard = "unknown"
        if os.path.isfile(os.path.join(basepath, entry)):
            with open(os.path.join(basepath, entry), newline='') as calibrations:
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

    df = pd.DataFrame(calibration_list)

    df[['Charge', 'm/z', 'Slope']] = df[['Charge', 'm/z', 'Slope']].apply(pd.to_numeric, errors='coerce')

    return df


def _2D_KDE(data, mzbandwidth, mzsamplerate, slopebandwidth, slopesamplerate):
    a = (1 / (mzbandwidth * np.sqrt(2 * np.pi)))
    b = ((0) / mzbandwidth * np.sqrt(2))
    density = a * np.exp(-0.5 * b * b)

    mz_data = data[0]
    slope_data = data[1]

    mzcutoff = density / 100

    a = (1 / (slopebandwidth * np.sqrt(2 * np.pi)))
    b = ((0) / slopebandwidth * np.sqrt(2))
    density = a * np.exp(-0.5 * b * b)

    slopecutoff = density / 100

    # print(slopecutoff*mzcutoff)
    min_mzdata = min(mz_data)
    max_mzdata = max(mz_data)
    mz_d = np.arange(min_mzdata, max_mzdata, mzsamplerate)

    min_slopedata = min(slope_data)
    max_slopedata = max(slope_data)
    slope_d = np.arange(min_slopedata, max_slopedata, slopesamplerate)

    matrix = np.zeros((len(mz_d), len(slope_d)))

    # mzspace first!

    for i in range(len(mz_data)):

        mz = mz_data[i]
        slope = slope_data[i]

        mz_index = (int)((mz - min_mzdata) / mzsamplerate)
        slope_index = (int)((slope - min_slopedata) / slopesamplerate)

        current_mz_index = mz_index
        current_slope_index = slope_index
        # print(f"{mz_data[i]} {slope_data[i]}  {mz_index} {slope_index}")
        while current_slope_index >= 0:

            slope_to_Check = slope - (slope_index - current_slope_index) * slopesamplerate

            a = (1 / (slopebandwidth * np.sqrt(2 * np.pi)))
            b = ((slope_to_Check - slope) / slopebandwidth * np.sqrt(2))
            slopedensity = a * np.exp(-0.5 * b * b)
            current_mz_index = mz_index

            # look left in mzspace
            while current_mz_index >= 0:

                MZtoCheck = mz - (mz_index - current_mz_index) * mzsamplerate

                a = (1 / (mzbandwidth * np.sqrt(2 * np.pi)))
                b = ((MZtoCheck - mz) / mzbandwidth * np.sqrt(2))
                mzdensity = a * np.exp(-0.5 * b * b)

                density = slopedensity * mzdensity

                # print(f"{MZtoCheck} vs {data[i]} {data[i] - MZtoCheck} {index - current_index} {density}")
                if mzdensity < mzcutoff:
                    break
                # print(f"adding to matrix {current_slope_index} {current_mz_index} {density}")
                matrix[current_mz_index][current_slope_index] += density
                current_mz_index -= 1

            # look right in mz space
            current_mz_index = mz_index + 1
            while current_mz_index < len(mz_d):

                MZtoCheck = mz - (current_mz_index - mz_index) * mzsamplerate

                a = (1 / (mzbandwidth * np.sqrt(2 * np.pi)))
                b = ((MZtoCheck - mz) / mzbandwidth * np.sqrt(2))
                mzdensity = a * np.exp(-0.5 * b * b)

                density = slopedensity * mzdensity

                # print(f"{MZtoCheck} vs {data[i]} {data[i] - MZtoCheck} {index - current_index} {density}")
                if mzdensity < mzcutoff:
                    break
                # print(f"adding to matrix {current_slope_index} {current_mz_index} {density}")
                matrix[current_mz_index][current_slope_index] += density
                current_mz_index += 1

            current_slope_index -= 1

        current_slope_index = slope_index + 1
        while current_slope_index < len(slope_d):

            slope_to_Check = slope - (slope_index - current_slope_index) * slopesamplerate

            a = (1 / (slopebandwidth * np.sqrt(2 * np.pi)))
            b = ((slope_to_Check - slope) / slopebandwidth * np.sqrt(2))
            slopedensity = a * np.exp(-0.5 * b * b)

            current_mz_index = mz_index

            # look left in mzspace
            while current_mz_index >= 0:

                MZtoCheck = mz - (mz_index - current_mz_index) * mzsamplerate

                a = (1 / (mzbandwidth * np.sqrt(2 * np.pi)))
                b = ((MZtoCheck - mz) / mzbandwidth * np.sqrt(2))
                mzdensity = a * np.exp(-0.5 * b * b)

                density = slopedensity * mzdensity

                # print(f"{MZtoCheck} vs {data[i]} {data[i] - MZtoCheck} {index - current_index} {density}")
                if mzdensity < mzcutoff:
                    break
                # print(f"adding to matrix {current_slope_index} {current_mz_index} {density}")
                matrix[current_mz_index][current_slope_index] += density
                current_mz_index -= 1

            # look right in mz space
            current_mz_index = mz_index + 1
            while current_mz_index < len(mz_d):

                MZtoCheck = mz - (current_mz_index - mz_index) * mzsamplerate

                a = (1 / (mzbandwidth * np.sqrt(2 * np.pi)))
                b = ((MZtoCheck - mz) / mzbandwidth * np.sqrt(2))
                mzdensity = a * np.exp(-0.5 * b * b)

                density = slopedensity * mzdensity

                # print(f"{MZtoCheck} vs {data[i]} {data[i] - MZtoCheck} {index - current_index} {density}")
                if mzdensity < mzcutoff:
                    break
                # print(f"adding to matrix {current_slope_index} {current_mz_index} {density}")
                matrix[current_mz_index][current_slope_index] += density
                current_mz_index += 1

            current_slope_index += 1
        print(float(i / len(mz_data) * 100))

    X, Y = np.meshgrid(slope_d, mz_d)

    print(X.shape)
    print(Y.shape)

    print(len(mz_d))
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    print(matrix.shape)
    surf = ax.plot_surface(X, Y, matrix, rstride=1, cstride=1, cmap='hot')
    plt.show()


if __name__ == '__main__':
    import time

    basepath = r"C:\Data\BOB\othercalibrationfolders\test2"

    df = CreateDF(basepath)

    spectra = [row['Slope'] for index, row in df.iterrows()]
    ordered_spectra = sorted(spectra)

    mzspectra = [row['m/z'] for index, row in df.iterrows()]
    mzordered_spectra = sorted(mzspectra)

    _2D_KDE((mzordered_spectra, ordered_spectra), 0.02, 0.1, 5000, 2000)

    """    bins = 50
    print("data loaded")
    start_t = time.time()
    X,Y =better_norm_KDE(ordered_spectra,1000,1000)
    plt.plot(X,Y)

    n, bins, patches = plt.hist(ordered_spectra, bins=bins, histtype='step', density=True,
                                label="histogram bins = {0}".format(bins))
    plt.show()"""

    """ min = ordered_spectra[0]
    max = ordered_spectra[-1]




    print(n)
    hist_max = np.max(n)
    print(bins)
    plt.scatter(ordered_spectra, [-0.000002] * len(ordered_spectra), marker="+", color="blue", label="raw data")

    x_d = np.arange(min, max, 1000)

    y = []
    x = []
    count = 0
    start_index = 0
    for i in range(len(x_d)):

        y_add, start_index = (norm_KDE(x_d[i], ordered_spectra, 20000, start_index))
        if y_add > 0:
            y.append(y_add)
            x.append(x_d[i])
        count += 1
        if count % int(len(x_d) / 1000) == 0:
            print(float(count / len(x_d)) * 100)
    print(100)
    print(time.time() - start_t)

    plt.plot(x, (y/np.max(y))*hist_max)
    plt.show()"""

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