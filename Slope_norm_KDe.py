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

slope_mean = 49573.12491063046
intercept = 36787.940439823215 #intiaial


intercept = 31538.62741


def bay(data):
    mzspectra, slpopespectra = data

    data = [np.array((data[0][i], data[1][i])) for i in range(len(data[0]))]
    data = np.array(data)

    from sklearn.naive_bayes import GaussianNB
    model = GaussianNB()
    model.fit(data, y);


def Clustering(data):
    slpopespectra,mzspectra  = data


    data = [np.array((data[1][i],data[0][i])) for i in range(len(data[0]))]
    data =  np.array(data)

    clusters = []
    from sklearn.cluster import KMeans
    for i in range(1, 10):
        km = KMeans(n_clusters=i).fit(data)
        clusters.append(km.inertia_)


    fig, ax = plt.subplots(figsize=(12, 8))
    sns.lineplot(x=list(range(1, 10)), y=clusters, ax=ax)
    ax.set_title('Searching for Elbow')
    ax.set_xlabel('Clusters')
    ax.set_ylabel('Inertia')

    plt.show()
    print(data)

    kmeans = KMeans(n_clusters=2)
    kmeans.fit(data)
    y_kmeans = kmeans.predict(data)

    plt.scatter(data[:, 0], data[:, 1], c=y_kmeans, s=1, cmap='viridis')

    centers = kmeans.cluster_centers_
    plt.scatter(centers[:, 0], centers[:, 1], c='black', s=200, alpha=0.5);

    plt.show()


    """    clusters = []
    from sklearn.cluster import KMeans
    for i in range(1, 10):
        km = KMeans(n_clusters=i).fit(data)
        clusters.append(km.inertia_)

    print(clusters)


    fig, ax = plt.subplots(figsize=(12, 8))
    sns.lineplot(x=list(range(1, 10)), y=clusters, ax=ax)
    ax.set_title('Searching for Elbow')
    ax.set_xlabel('Clusters')
    ax.set_ylabel('Inertia')

    plt.show()


    km3 = KMeans(n_clusters=2).fit(data)
    plt.figure(figsize=(12, 8))
    sns.scatterplot(data[0], data[1], hue=km3.labels_,
                    palette=sns.color_palette('hls', 3))
    plt.title('KMeans with 2 Clusters')
    plt.show()"""


def SimpleDF(basepath,minmz,maxmz,overwrite=False,minslope=1000,maxslope=40000000):
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

def _2D_KDE_fixed(data,mzbandwidth,mzsamplerate,slopebandwidth,slopesamplerate):


    mz_data = data[0]
    slope_data = data[1]

    resolution = (140000,200)
    NormalKDEBandwithRelationship = 0.1963

    #print(slopecutoff*mzcutoff)

    min_mzdata=min(mz_data)
    max_mzdata=max(mz_data)


    print(min_mzdata,max_mzdata)
    mz_d = np.arange(min_mzdata, max_mzdata, mzsamplerate)

    min_slopedata = min(slope_data)
    max_slopedata = max(slope_data)
    slope_d = np.arange(min_slopedata, max_slopedata, slopesamplerate)

    matrix = np.zeros((len(mz_d),len(slope_d)))
    print(matrix.shape)



    #mzspace first!
    print("starting")
    for i in range(len(mz_data)):
        mz = mz_data[i]
        slope = slope_data[i]


        a = (1 / (mzbandwidth * np.sqrt(2 * np.pi)))
        b = ((0) / mzbandwidth * np.sqrt(2))
        density = a * np.exp(-0.5 * b * b)


        mzcutoff = density / 100


        a = (1 / (slopebandwidth * np.sqrt(2 * np.pi)))
        b = ((0) / slopebandwidth * np.sqrt(2))
        density = a * np.exp(-0.5 * b * b)

        slopecutoff = density / 100


        mz_index = (int)((mz - min_mzdata) / mzsamplerate)



        slope_index = (int)((slope - min_slopedata) / slopesamplerate)

        current_mz_index = mz_index
        current_slope_index = slope_index
        #print(f"{mz_data[i]} {slope_data[i]}  {mz_index} {slope_index}")
        while current_slope_index >= 0:

            slope_to_Check = slope - (slope_index - current_slope_index) * slopesamplerate

            a = (1 / (slopebandwidth * np.sqrt(2 * np.pi)))
            b = ((slope_to_Check - slope) / slopebandwidth * np.sqrt(2))
            slopedensity = a * np.exp(-0.5 * b * b)
            current_mz_index = mz_index
            if slopedensity < slopecutoff:
                break
            #look left in mzspace
            while current_mz_index >= 0:

                MZtoCheck = mz - (mz_index - current_mz_index)*mzsamplerate


                a = (1 / (mzbandwidth * np.sqrt(2 * np.pi)))
                b = ((MZtoCheck - mz) / mzbandwidth * np.sqrt(2))
                mzdensity  = a * np.exp(-0.5 * b * b)



                density = slopedensity*mzdensity

                #print(f"{MZtoCheck} vs {data[i]} {data[i] - MZtoCheck} {index - current_index} {density}")
                if mzdensity < mzcutoff:
                    break
                #print(f"adding to matrix {current_slope_index} {current_mz_index} {density}")
                matrix[current_mz_index][current_slope_index] += density
                current_mz_index -= 1


            #look right in mz space
            current_mz_index = mz_index + 1
            while current_mz_index < len(mz_d):

                MZtoCheck = mz - (current_mz_index - mz_index)*mzsamplerate

                a = (1 / (mzbandwidth * np.sqrt(2 * np.pi)))
                b = ((MZtoCheck - mz) / mzbandwidth * np.sqrt(2))
                mzdensity = a * np.exp(-0.5 * b * b)



                density = slopedensity * mzdensity

                # print(f"{MZtoCheck} vs {data[i]} {data[i] - MZtoCheck} {index - current_index} {density}")
                if mzdensity < mzcutoff:
                    break
                #print(f"adding to matrix {current_slope_index} {current_mz_index} {density}")
                matrix[current_mz_index][current_slope_index] += density
                current_mz_index += 1

            current_slope_index -= 1

        current_slope_index = slope_index + 1
        while current_slope_index < len(slope_d):

            slope_to_Check = slope - (slope_index - current_slope_index) * slopesamplerate

            a = (1 / (slopebandwidth * np.sqrt(2 * np.pi)))
            b = ((slope_to_Check - slope) / slopebandwidth * np.sqrt(2))
            slopedensity = a * np.exp(-0.5 * b * b)
            if slopedensity < slopecutoff:
                break
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
                #print(f"adding to matrix {current_slope_index} {current_mz_index} {density}")
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
                #print(f"adding to matrix {current_slope_index} {current_mz_index} {density}")
                matrix[current_mz_index][current_slope_index] += density
                current_mz_index += 1

            current_slope_index += 1
        if i % 1000 == 0:
            print(float(i/len(mz_data)*100))






    with open(r"C:\Data\BOB\Data for Mike\UbCAEnoMyoMIX_pingoutput_txt\matrix.pkl", "wb") as writePKL:
        pickle.dump((matrix,slope_d,mz_d), writePKL)

    return matrix,slope_d,mz_d



def _2D_KDE(data,mzsamplerate,slopebandwidth,slopesamplerate):


    mz_data = data[0]
    slope_data = data[1]

    resolution = (140000,200)
    NormalKDEBandwithRelationship = 0.1963

    #print(slopecutoff*mzcutoff)

    min_mzdata=min(mz_data)
    max_mzdata=max(mz_data)

    mz_d = np.arange(min_mzdata, max_mzdata, mzsamplerate)

    min_slopedata = min(slope_data)
    max_slopedata = max(slope_data)
    slope_d = np.arange(min_slopedata, max_slopedata, slopesamplerate)

    matrix = np.zeros((len(mz_d),len(slope_d)))




    #mzspace first!
    print("starting")
    for i in range(len(mz_data)):
        mz = mz_data[i]
        slope = slope_data[i]
        mzbandwidth = mz/(resolution[0]* np.sqrt(resolution[1]/mz))*NormalKDEBandwithRelationship

        a = (1 / (mzbandwidth * np.sqrt(2 * np.pi)))
        b = ((0) / mzbandwidth * np.sqrt(2))
        density = a * np.exp(-0.5 * b * b)


        mzcutoff = density / 1000


        a = (1 / (slopebandwidth * np.sqrt(2 * np.pi)))
        b = ((0) / slopebandwidth * np.sqrt(2))
        density = a * np.exp(-0.5 * b * b)

        slopecutoff = density / 1000


        mz_index = (int)((mz - min_mzdata) / mzsamplerate)



        slope_index = (int)((slope - min_slopedata) / slopesamplerate)

        current_mz_index = mz_index
        current_slope_index = slope_index
        #print(f"{mz_data[i]} {slope_data[i]}  {mz_index} {slope_index}")
        while current_slope_index >= 0:

            slope_to_Check = slope - (slope_index - current_slope_index) * slopesamplerate

            a = (1 / (slopebandwidth * np.sqrt(2 * np.pi)))
            b = ((slope_to_Check - slope) / slopebandwidth * np.sqrt(2))
            slopedensity = a * np.exp(-0.5 * b * b)
            current_mz_index = mz_index
            if slopedensity < slopecutoff:
                break
            #look left in mzspace
            while current_mz_index >= 0:

                MZtoCheck = mz - (mz_index - current_mz_index)*mzsamplerate


                a = (1 / (mzbandwidth * np.sqrt(2 * np.pi)))
                b = ((MZtoCheck - mz) / mzbandwidth * np.sqrt(2))
                mzdensity  = a * np.exp(-0.5 * b * b)



                density = slopedensity*mzdensity

                #print(f"{MZtoCheck} vs {data[i]} {data[i] - MZtoCheck} {index - current_index} {density}")
                if mzdensity < mzcutoff:
                    break
                #print(f"adding to matrix {current_slope_index} {current_mz_index} {density}")
                matrix[current_mz_index][current_slope_index] += density
                current_mz_index -= 1


            #look right in mz space
            current_mz_index = mz_index + 1
            while current_mz_index < len(mz_d):

                MZtoCheck = mz - (current_mz_index - mz_index)*mzsamplerate

                a = (1 / (mzbandwidth * np.sqrt(2 * np.pi)))
                b = ((MZtoCheck - mz) / mzbandwidth * np.sqrt(2))
                mzdensity = a * np.exp(-0.5 * b * b)



                density = slopedensity * mzdensity

                # print(f"{MZtoCheck} vs {data[i]} {data[i] - MZtoCheck} {index - current_index} {density}")
                if mzdensity < mzcutoff:
                    break
                #print(f"adding to matrix {current_slope_index} {current_mz_index} {density}")
                matrix[current_mz_index][current_slope_index] += density
                current_mz_index += 1

            current_slope_index -= 1

        current_slope_index = slope_index + 1
        while current_slope_index < len(slope_d):

            slope_to_Check = slope - (slope_index - current_slope_index) * slopesamplerate

            a = (1 / (slopebandwidth * np.sqrt(2 * np.pi)))
            b = ((slope_to_Check - slope) / slopebandwidth * np.sqrt(2))
            slopedensity = a * np.exp(-0.5 * b * b)
            if slopedensity < slopecutoff:
                break
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
                #print(f"adding to matrix {current_slope_index} {current_mz_index} {density}")
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
                #print(f"adding to matrix {current_slope_index} {current_mz_index} {density}")
                matrix[current_mz_index][current_slope_index] += density
                current_mz_index += 1

            current_slope_index += 1
        if i % 1000 == 0:
            print(float(i/len(mz_data)*100))






    with open(r"C:\Data\BOB\Data for Mike\UbCAEnoMyoMIX_pingoutput_txt\matrix.pkl", "wb") as writePKL:
        pickle.dump((matrix,slope_d,mz_d), writePKL)

    return matrix,slope_d,mz_d


def LoadMatrix(matrixFile):
    with open(matrixFile,"rb") as matrixfile:
        matrix, slope_d, mz_d = pickle.load(matrixfile)
    return matrix, slope_d, mz_d

def ShowFigs(matrix, slope_d, mz_d):
    X, Y = np.meshgrid(slope_d, mz_d)

    fig, ax = plt.subplots()
    axim    =ax.contour(X,Y,matrix,cmap='YlOrRd')

    cb = fig.colorbar(axim)

    ax.set_xlabel("Slope",fontsize=24)
    ax.set_ylabel("m/z (Da.)",fontsize=24)
    ax.tick_params(axis='y',labelsize= 20)
    ax.tick_params(axis='x',  labelsize=20)
    plt.show()




def better_norm_KDE(data,bandwidth,sample_rate):


    a = (1 / (bandwidth * np.sqrt(2 * np.pi)))
    b = ((0) / bandwidth * np.sqrt(2))
    density = a * np.exp(-0.5 * b * b)

    cutoff = density/1000


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



def GetSubsample(matrix, slope_d, mz_d,minmz,maxmz,minslope,maxslope):
    minmzindex = 0
    maxmzindex = 0
    minslopeindex = 0
    maxslopeindex = 0

    for i in range(len(mz_d)):
        if mz_d[i] >= minmz and minmzindex==0:
            minmzindex = i

        if mz_d[i] >= maxmz and maxmzindex==0:
            maxmzindex = i


    for i in range(len(slope_d)):
        if slope_d[i] >= minslope and minslopeindex==0:
            minslopeindex = i

        if slope_d[i] >= maxslope and maxslopeindex==0:
            maxslopeindex = i

    if maxslopeindex == 0:
        maxslopeindex = len(slope_d)-1


    if maxmzindex == 0:
        maxmzindex = len(mz_d)-1

    matrix2 = matrix[minmzindex:maxmzindex,minslopeindex:maxslopeindex]


    return matrix2, slope_d[minslopeindex:maxslopeindex],mz_d[minmzindex:maxmzindex]

def GetPeaks(matrix,cutoff=0.001):
    peaks = []
    total = (len(mz_d)-2)*(len(slope_d)-2)
    count = 0
    print("find peaks and remove noise")
    for i in range(1,len(matrix)-1):
        for j in range(1,len(matrix[i])-1):
            if matrix[i][j] < cutoff:
                matrix[i][j] = 0
            elif matrix[i][j] > matrix[i+1][j] and matrix[i][j] > matrix[i-1][j] and matrix[i][j] > matrix[i][j+1] and matrix[i][j] > matrix[i][j-1]:
                #print(f"peak found @ mz:{mz_d[i]} slope:: {slope_d[j]} intensity: {matrix[i][j]}")
                peaks.append((mz_d[i],slope_d[j],matrix[i][j]))
            count +=1
            if count % 1000 == 0:
                print((float(count)/float(total))*100)
    print("done")
    with open(r"C:\Data\BOB\Data for Mike\UbCAEnoMyoMIX_pingoutput_txt\peaks.pkl","wb") as peakspkl:
        pickle.dump(peaks,peakspkl)
    return peaks


if __name__ == '__main__':

    import time

    #basepath = r"C:\Data\BOB\othercalibrationfolders\test"
    IgGs = r"C:\Data\I2MS\IgG\20200721_jpm907_2208-rbd_002_data_csv"
    test = r"C:\Data\BOB\Data for Mike\UbCAEnoMyoMIX_pingoutput_txt"
    #df = CreateDF(basepath)

    mzspectra ,slpopespectra = SimpleDF(test,200,4000,True)


    mzspectra = np.array(mzspectra)
    slpopespectra = np.array(slpopespectra)

    #Clustering((mzspectra ,slpopespectra))
   # _2D_KDE((mzspectra ,slpopespectra),0.02,100000,5000)
    #matrix, slope_d, mz_d = LoadMatrix(r"C:\Data\BOB\Data for Mike\UbCAEnoMyoMIX_pingoutput_txt\matrix.pkl")
    #ShowFigs(matrix, slope_d, mz_d)
    print("data loaded")

    """data = matrix.flatten()
    plt.hist(data,10000)
    plt.show()"""
    #ShowFigs(matrix, slope_d, mz_d)
    """    print("remove noise")
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] < 0.003:
                matrix[i][j] = 0
    print("noise removed")"""

    """    #ShowFigs(matrix, slope_d, mz_d)
    peaks = []
    total = (len(mz_d)-2)*(len(slope_d)-2)
    count = 0
    print("find peaks and remove noise")
    for i in range(1,len(matrix)-1):
        for j in range(1,len(matrix[i])-1):
            if matrix[i][j] < 0.003:
                matrix[i][j] = 0
            elif matrix[i][j] > matrix[i+1][j] and matrix[i][j] > matrix[i-1][j] and matrix[i][j] > matrix[i][j+1] and matrix[i][j] > matrix[i][j-1]:
                #print(f"peak found @ mz:{mz_d[i]} slope:: {slope_d[j]} intensity: {matrix[i][j]}")
                peaks.append((mz_d[i],slope_d[j],matrix[i][j]))
            count +=1
            print((float(count)/float(total))*100)
    print("done")

    print(len(peaks))


    with open(r"C:\Data\I2MS\IgG\20200721_jpm907_2208-rbd_002_data_csv\peaks.pkl","wb") as peakspkl:
        pickle.dump(peaks,peakspkl)"""

    #intensity = [sum([matrix[i][j] for j in range(len(matrix[i]))])for i in range(len(matrix))]



    #fig, (ax1,ax2) = plt.subplots(1,2)
    #ax1.plot(mz_d, intensity)




    #matrix2, slope_d2, mz_d2 = GetSubsample(matrix,slope_d,mz_d,900,980,1167000,3000000)



    #ShowFigs(matrix, slope_d, mz_d)



    """    
    data = matrix.flatten()
    median = np.median(data)
    print(median)
    data = [data[i] for i in range(len(data)) if 0 < data[i] < median]

    plt.hist(data,10000)
    plt.show()"""



    #peaks = GetPeaks(matrix,0.000) #create peaks

    #load peaks
    with open(r"C:\Data\BOB\Data for Mike\UbCAEnoMyoMIX_pingoutput_txt\peaks.pkl","rb") as peakspkl:
        peaks = pickle.load(peakspkl)



    X = [peaks[i][0] for i in range(len(peaks))]
    Y = [peaks[i][1] for i in range(len(peaks))]
    Z = [peaks[i][2] for i in range(len(peaks))]
    maxint = max(Z)

    Z = [(Z[i]/maxint)*100 for i in range(len(Z))]


    Xbins = list(range(int(min(X))-1,int(max(X)),200))
    Ybins = list(range(int(min(Y)),int(max(Y)),200000))

    total= len(Xbins)*len(Ybins)
    count = 0
    all_peaks = []
    for i in range(len(Xbins)-1):
        for j in range(len(Ybins)-1):
            Xmax=Xbins[i+1]
            Xmin = Xbins[i]
            Ymax = Ybins[j + 1]
            Ymin = Ybins[j]


            current_region = [(X[x],Y[x],Z[x]) for x in range(len(Z)) if Xmin < X[x] < Xmax and Ymin < Y[x] < Ymax]
            if len(current_region)> 1:
                noise = np.median([current_region[x][2] for x in range(len(current_region))])
                #region_filted= [current_region[x] for x in range(len(current_region)) if noise*1 < current_region[x][2]]
                region_filted = [current_region[x] for x in range(len(current_region)) ]
                count += 1
                print(float((count)/total))
                all_peaks+= region_filted

    X = [all_peaks[i][0] for i in range(len(all_peaks))]
    Y = [all_peaks[i][1] for i in range(len(all_peaks))]
    Z = [all_peaks[i][2] for i in range(len(all_peaks))]
    with open(r"C:\Data\BOB\Data for Mike\test\2dkdeoutx1.csv","w") as writefile:
        writefile.write("Scan\tIonnumber\tSegnumber\tM/Z\tFrequency\tSlope\tSlope R Squared\tTime of Birth\tTime of Death\tIntensity\n")
        for i in range(len(X)):
            writefile.write(f"0\t0\t0\t{X[i]}\t0\t{Y[i]}\t1\t0\t2\t{Z[i]}\n")


    #noise = [Z[i] for i  in range(len(Z)) if  Z[i] < max(Z)/1000 ]
    print("done")

    n, bins, patches = plt.hist(Z,1000)
    plt.show()


    
    max_int = max(Z)


    plt.vlines(X,0,Z)
    plt.show()


    fig, ax = plt.subplots()
    #axim = ax.contour(X2,Y2,matrix2, cmap='YlOrRd')

    #cb = fig.colorbar(axim)

    ax.set_xlabel("Slope", fontsize=24)
    ax.set_ylabel("m/z (Da.)", fontsize=24)
    ax.tick_params(axis='y', labelsize=20)
    ax.tick_params(axis='x', labelsize=20)

    sizes = [(Z[i]/max_int)*1000 for i in range(len(Z)) ]

    plt.scatter(Y,X,sizes,color="g",marker="+")


    startsinglemz = 812.5
    startCharge = 36
    startSlope = 1980000
    deltaslope = 50000


    X2 = []
    Y2=[]
    mass = (startsinglemz*startCharge)-(1.00728*startCharge)
    for i in range(20):
        try:
            currentCharge = startCharge - i
            currentMZ = (mass+(currentCharge*1.00728))/currentCharge
            currentSlope = startSlope - i*deltaslope
            X2.append(currentSlope)
            Y2.append(currentMZ)
        except:
            ""
    plt.scatter(X2,Y2,color="r",s=200,marker="o")




    startsinglemz = 832.90
    startCharge = 50
    startSlope = 2540000
    deltaslope = 50000


    X3 = []
    Y3=[]
    mass = (startsinglemz*startCharge)-(1.00728*startCharge)
    for i in range(35):
        currentCharge = startCharge - i
        currentMZ = (mass+(currentCharge*1.00728))/currentCharge
        currentSlope = startSlope - i*deltaslope
        X3.append(currentSlope)
        Y3.append(currentMZ)
    plt.scatter(X3,Y3,color="b",s=200,marker="o")

    plt.show()

    #peaks_in_window = [peaks[i] for i in range(len(peaks)) if 920 < peaks[i][0] < 955 if 1100000 < peaks[i][1] < 1500000] #LC
    #peaks_in_window = [peaks[i] for i in range(len(peaks)) if 1100< peaks[i][0] < 1600 if 1600000 < peaks[i][1] < 2500000] # part of heavy chain 4-5 charge groups
    peaks_in_window = [peaks[i] for i in range(len(peaks)) ]
    print(len(peaks_in_window))

    X = [(peaks_in_window[i][0]) for i in range(len(peaks_in_window))]
    Y = [peaks_in_window[i][1] for i in range(len(peaks_in_window))]
    Z = [peaks_in_window[i][2] for i in range(len(peaks_in_window))]

    X_Mass = []


    plt.vlines(X,0,Z)
    plt.show()
    """    window = 3

    Y2 = [sum(peaks_in_window[i+j][1] for j in range(window+1))/window for i in range(len(peaks_in_window)-window+1)]
    last_point = (peaks_in_window[-1][1]+peaks_in_window[-2][1]+peaks_in_window[-3][1])/3
    Y2.append(last_point)
    last_point = (peaks_in_window[-2][1]+peaks_in_window[-3][1]+peaks_in_window[-4][1])/3
    Y2.append(last_point)"""



    for i in range(len(Y)):
        slope = Y[i]

        charge = round(( slope - intercept)/slope_mean)


        X_Mass.append( ((peaks_in_window[i][0]) * charge) - 1.00728 * charge)


    fig, (ax1,ax2) = plt.subplots(2,1)

    X = [((peaks_in_window[i][0])*26) -1.00728*26 for i in range(len(peaks_in_window))]

    ax1.vlines(X,0,Z)



    ax2.vlines(X_Mass,0,Z)
    plt.show()

    """    intensity = [sum([matrix[i][j] for j in range(len(matrix[i]))]) for i in range(len(matrix))]

    print(len(mz_d))
    print(len(intensity))

    ax2.plot(mz_d, intensity)
    plt.show()"""


    #_2D_KDE((mzspectra ,slpopespectra),0.5,8000,10000)
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