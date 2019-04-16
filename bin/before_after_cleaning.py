#!/usr/bin/env python3

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


sns.set(style="ticks", palette="pastel")

# Load the dataset
before_after_cleaning = pd.read_csv('../data/before_after_cleaning.csv')

ax = sns.catplot(x='Status', y='FPKM', hue='Type', kind='box', height=5, aspect=1.3,
                 data=before_after_cleaning, legend=False)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

#plt.title('H3F3A FPKM Comparison')
#plt.show()
plt.savefig('../data/noclean_vs_clean.png', dpi=1000, bbox_inches="tight")
