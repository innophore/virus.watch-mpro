#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 08:53:34 2022

Produces a single frame output with the data from mpro_escape_movie.py
Decoupled from the main code in order to properly reset memory usage 
after each frame

@author: andreas.krassnigg
"""

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
import sys


# get number of command-line-arguments given
n = len(sys.argv)
# 0 is this file's name, but 1 is the one we want
if n>0:
    file_suffix = sys.argv[1]
else:
    file_suffix = ""


# load input dict from file
in_file = open("temp_"+file_suffix+".pkl", "rb")
frame_data  = pickle.load(in_file)
in_file.close()


datum       = frame_data['date']
datum_index = frame_data['date_index']
this_data   = frame_data['data']
positions   = frame_data['positions']
letters_new = frame_data['letters']
refsize     = frame_data['refsize']
refsize2    = frame_data['refsize2']
color_list  = frame_data['colors']
labels      = frame_data['labels']
path        = frame_data['path']
subpath     = frame_data['subpath']

if subpath != file_suffix:
    print("Warning: subpath vs. file_suffix mismatch!")


figure_size = (100, 25)
# create the figure for this timeframe
fig = plt.figure(figsize=figure_size, dpi=100)
ax = fig.add_subplot(111)
width = 0.45  # the width of the bars: can also be len(x) sequence

# loop over all letters for mutation at a given position
for let_ind, a_let in enumerate(letters_new):
# for let_ind, a_let in enumerate(letters_new[-1]):
    
    # grab only relevant part of data for bars
    bar_data = this_data[:, let_ind]
    
    # replot the bar plots for this letter
    # bars01 += plt.bar(positions, bar_data, width,
    #                       bottom=np.sum(this_data[:, :let_ind], axis=-1), color=color_list[let_ind])
    
    plt.bar(positions, bar_data, width,
                          bottom=np.sum(this_data[:, :let_ind], axis=-1),
                          label=a_let + ' / ' + 
                          '{:.1f}'.format(round((np.sum(bar_data) / refsize2) * 100, 1)) + 
                          '%', color=color_list[let_ind])
    
# replot the title to reflect the date
plt.title(subpath+' - unique mutations - %s' % str(datum)[0:4] + '-' + str(datum)[4:6] + '-' + str(datum)[6:],
  fontsize=100)
    
# plot the legend
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=50)
handles1, labels1 = ax.get_legend_handles_labels()
ax.legend(handles1[::-1], labels1[::-1], loc='center left', bbox_to_anchor=(1, 0.5), fontsize=50)

ax.set_ylim([0, 22])
plt.xlabel('Positions', size=50)
plt.ylabel('Number of mutations', size=80)
pos_y = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22]
plt.yticks(pos_y, size=80)
plt.xticks(positions, labels, size=200 / (refsize ** (1 / 2)))


# plt.figtext(.15, .65, "Average mut. p. position: %s\nSTD: %s\nMedian: %s\nAverage mut. p. unique sequence: "
#                     "%s\nSTD: %s\nMedian: %s\n" % (average_single, std_single, median_single, av_gene_single, std_gene_single, med_gene_single), size=50)
# plt.figtext(.15, .65, "Average mut. p. position: %s\nSTD: %s\nMedian: %s" % (
# average_single, std_single, median_single), size=50)
fig.set_size_inches(figure_size[0], figure_size[1], forward=True)

plt.savefig(os.path.join(path, 'PLOTS', 'movies', subpath, str(datum)+'.png'))
# plt.savefig(os.path.join(path, str(datum)+'.png'))
plt.close(fig)



