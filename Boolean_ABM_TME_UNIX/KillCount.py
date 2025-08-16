import pandas as pd
import numpy as np
import math
import os
import matplotlib.pyplot as plt
import seaborn as sns
import sys

'''
This script finds the average kill count per day over 10 replicates of a given simulation
User gives the folder location of the 10 sets
Script outputs a csv of these averages as well as a png figure
'''

filepath = sys.argv[1] # like "NoProlif_All_base"

df = []

for i in range (1,11):
    file = filepath + "/set_" + str(i) + "/TimeSeries.csv"

    df1 = pd.read_csv(file)
    
    df.append(df1)

dfm = pd.concat(df[0:10]).groupby(level=0).mean()


p = sns.lineplot(data=dfm, x="time_step",y="killCount",label="Kill Count")
sns.lineplot(data=dfm, x="time_step", y="cd8_activate",ax=p,label="CD8 Activated")
sns.lineplot(data=dfm, x="time_step", y="cd8_suppressed",ax=p,label="CD8 Suppressed")
p.set_ylim(0, dfm["cancer"][0])
p.legend()
p.set_title(str(filepath))

p_fig = p.get_figure()
p_fig.savefig(filepath+"/" + filepath+ "_Avg.png")

dfm.to_csv(filepath+"/" + filepath+ "_Avg.csv")