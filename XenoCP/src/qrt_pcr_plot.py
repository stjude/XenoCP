#!/usr/bin/env python3

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


sns.set(style="ticks", palette="pastel")

# Load the dataset
qrt_pcr = pd.read_csv('./XenoCP/data/qRT-PCR.csv')
tips = sns.load_dataset("tips")

ax = sns.boxplot(x="Endpoint tumors analyzed (mice)", y="Expression (logCq)", data=qrt_pcr)
ax = sns.swarmplot(x="Endpoint tumors analyzed (mice)", y="Expression (logCq)", data=qrt_pcr, color=".25")
plt.title('qRT-PCR')
plt.show()
plt.savefig('./XenoCP/data/qRT-PCR.png')

