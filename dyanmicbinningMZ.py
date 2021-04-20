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
import pickle
from matplotlib.lines import Line2D

from mpl_toolkits.mplot3d import Axes3D


def SimpleDF(basepath,minmz,maxmz,overwrite=False,minslope=1000,maxslope=4000000):
    dirs = os.listdir(basepath)

    mzs = []
    slopes = []

    if "data.pkl" in dirs and not overwrite:
        print("loading pickle")
        with open(basepath + "/data.pkl","rb") as writePKL:
            mzspectra,slpopespectra = pickle.load(writePKL)

    else:
        for entry in dirs:

            if entry.split('.')[-1] == "txt" or entry.split('.')[-1] == "csv" :
                if os.path.isfile(os.path.join(basepath, entry)):
                    with open(os.path.join(basepath, entry), newline='') as calibrations:
                        line = calibrations.readline()
                        line = calibrations.readline()
                        while  line:
                            try:
                                mz = float(line.split('\t')[3])
                                slope = float(line.split('\t')[5])

                                if  minslope< slope < maxslope:

                                    mzs.append(mz)
                                    slopes.append(slope)
                            except:
                                ""
                            line = calibrations.readline()

                print(f"{entry }file Read")


        spectra = [(mzs[i],slopes[i]) for i in range(len(mzs)) if minmz <= mzs[i] <= maxmz]
        ordered_spectra = sorted(spectra)

        mzspectra = [ordered_spectra[x][0] for x in range(len(ordered_spectra))]

        slpopespectra = [ordered_spectra[x][1] for x in range(len(ordered_spectra))]
        with open(basepath + "/data.pkl","wb") as writePKL:
            pickle.dump((mzspectra,slpopespectra),writePKL)
        print("creating pickle")
    return mzspectra,slpopespectra



def better_norm_KDE(data,bandwidth,sample_rate):


    a = (1 / (bandwidth * np.sqrt(2 * np.pi)))
    b = ((0) / bandwidth * np.sqrt(2))
    density = a * np.exp(-0.5 * b * b)

    cutoff = density/100


    min_data=min(data)
    max_data=max(data)
    X_d = np.arange(min_data, max_data, sample_rate)
    sample_count = len(X_d)
    Y = [0]*len(X_d)

    for i in range(len(data)):

        index = (int)((data[i] - min_data) / sample_rate)

        current_index = index

        #look left
        while current_index >= 0:

            MZtoCheck = data[i] - (index - current_index)*sample_rate


            a = (1 / (bandwidth * np.sqrt(2 * np.pi)))
            b = ((MZtoCheck - data[i]) / bandwidth * np.sqrt(2))
            density  = a * np.exp(-0.5 * b * b)
            #print(f"{MZtoCheck} vs {data[i]} {data[i] - MZtoCheck} {index - current_index} {density}")
            if density < cutoff:
                break

            Y[current_index] += density
            current_index -= 1

        current_index = index + 1
        while current_index < len(X_d):

            MZtoCheck = X_d[( current_index)]

            a = (1 / (bandwidth * np.sqrt(2 * np.pi)))
            b = ((MZtoCheck - data[i]) / bandwidth * np.sqrt(2))
            density = a * np.exp(-0.5 * b * b)
            if density < cutoff:
                break
            Y[current_index] += density
            current_index += 1


        print(float(i)/len(data)*100)

    return (X_d,Y)


if __name__ == '__main__':

    import time

    #basepath = r"C:\Data\BOB\othercalibrationfolders\test"
    IgGs = r"C:\Data\I2MS\IgG\20200721_jpm907_2208-rbd_002_data_csv"
    #df = CreateDF(basepath)

    #mzspectra ,slpopespectra = SimpleDF(IgGs,200,4000,True)
    print("data loaded")

    #mz_data = np.array(mzspectra,dtype="float32")
    #slpope_data = np.array(slpopespectra,dtype="float32")

    resolution = (140000, 200)
    NormalKDEBandwithRelationship = 0.1963

    min_mzdata=min(mz_data)
    #max_mzdata=max(mz_data)
    bins = []
    current_bin_size = min_mzdata / (resolution[0] * np.sqrt(resolution[1] / min_mzdata)) * NormalKDEBandwithRelationship
    start_of_bin = min_mzdata
    end_of_bin = min_mzdata + current_bin_size


    while end_of_bin < max_mzdata:
        start_of_bin = end_of_bin
        current_bin_size = start_of_bin / (resolution[0] * np.sqrt(resolution[1] / start_of_bin)) * NormalKDEBandwithRelationship
        end_of_bin = start_of_bin + current_bin_size
        bins.append((start_of_bin, end_of_bin))

    j=0

    dict_bins = {}
    current_bin_data = []
    for i in range(len(mz_data)):

        if bins[j][0] <= mz_data[i] < bins[j][1]:
            current_bin_data.append()