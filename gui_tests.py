
# Imports
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

def treeclick(title,atr):
    selected = app.getTreeSelected("tree1")

xml = """            <stuff><ChargeVsSlope></ChargeVsSlope>
            <Mean_Vs_std></Mean_Vs_std></stuff>
            """


app=gui("Grid Demo", "300x300")
app.setSticky("news")
app.setExpand("both")
app.setFont(14)

app.addLabel("l1", "row=0\ncolumn=0")

fig = app.addPlotFig("p1", showNav=True,row=0,column=2,colspan=4,rowspan=2)

app.setLabelBg("l1", "red")

app.go()