#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates a figure to reference when specifying what color
to use of the default color palatte.
"""

import numpy as np
import matplotlib.pyplot as plt

pltcolors = plt.rcParams['axes.prop_cycle'].by_key()['color']

fig, axs = plt.subplots()

labels = range(len(pltcolors))  # indices
x = range(len(pltcolors))  # the label locations
values = np.full(len(pltcolors), 1)

for i in range(len(pltcolors)):
    axs.bar(x[i], values[i], color=pltcolors[i])
    axs.annotate(pltcolors[i], (x[i], 0.1), rotation=90)
axs.set_title('default matplotlib colors')
axs.set_xticks(x)  # make sure each bar has a tick
axs.tick_params(axis="x", which="both", bottom=False, top=False)
plt.tick_params(labelleft=False, left=False)
axs.set_xticklabels(labels)
axs.set_xlabel('index value')




