#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

# Plotting graphs for the simulator talk

# Read filenames
filename = "qsim_cnot.csv"
# Read the file into a dataframe
df = pd.read_csv(filename, comment='#', delim_whitespace=True)


# Plotting the graph
cols_list = list(range(1, len(df.columns)))
df.plot(0, y=cols_list, rot=1, logy=True, grid=True, style='o-', \
        xlabel='Number of qubits', ylabel='Time taken per gate (seconds)',
        title='Testing the controlled Not gate', linewidth=0.7)

plt.show()
