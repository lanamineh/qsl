#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

# Read filename from command line
filename = "pauliX.csv"
if(len(sys.argv) == 2):
    filename = sys.argv[1]

# Read the file into a dataframe
df = pd.read_csv(filename, comment='#', delim_whitespace=True)

# Read first line of file for some context
f = open(filename, "r")
heading = filename + " " + f.readline()
f.close()

# Plotting the graph
cols_list = list(range(1, len(df.columns)))
df.plot(0, y=cols_list, rot=1, logy=True, grid=True, style='o-', title=heading)

plt.show()
