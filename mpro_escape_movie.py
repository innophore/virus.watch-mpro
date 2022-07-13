#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2022

Takes two input files with sequence information and transforms them
into a collection of data in an .npy file
plus a set of frame pngs to be combined into a movie via ffmpeg



@author: lena.parigger
@author: andreas.krassnigg
@author: tobias.schopper


"""

import pandas as pd
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import os
import numpy as np
from datetime import date, timedelta
import datetime
# import shutil
# import cv2
from tqdm import tqdm
# import math
# import multiprocessing
# from multiprocessing import Pool
# from matplotlib.pyplot import figure
# import matplotlib.pyplot as plt
# from datetime import datetime
# import matplotlib.dates as mdates
import matplotlib.animation as animation
import pickle
import time
import sys
from multiprocessing import Pool

# ------------------------------------------------------------
# in-house imports
from utils.data_processing import build_first_occurrence_array

from mpro_escape_input import the_inputs


# get number of command-line-arguments given
n = len(sys.argv)
# 0 is this file's name, but 1 is the one we want
if n > 1:
    # take protein name from the command line
    subpath = sys.argv[1]
else:
    # take protein name from input file
    subpath = the_inputs['protein']



# define main dirpath
path = the_inputs['path']



# whether or not to use insertions and deletions as letters below
use_insertions_and_deletions = the_inputs['use_ins_dels']

# which column to use for mutation data
mutations_column_name = the_inputs['mutation_column']

CORRECT_DATES = the_inputs['correct_dates']


# set what output to create
MAKE_PNGS = the_inputs['make_pngs']

# specify frames to produce (more useful for debugging)
START_FRAME = the_inputs['start_frame']
USE_FRAMES = the_inputs['use_frames']


VERBOSE = the_inputs['verbose']     # set to False to suppress outputs 

SILENT = the_inputs['silent']    # shut off all outputs completely


if not SILENT: print("Escape-movie creation started on Protein " + subpath)

# create output file names
data_temp_file = subpath + "_all_data"

temp_file_suffix = "_2"

RECREATE_DATA = True  

# if the data file is already there, don't recreate the data and write out a warning
if os.path.exists(os.path.join(path, data_temp_file + ".npy")):
    RECREATE_DATA = False
    if not SILENT: print("WARNING: data target file already exists. Skipping data recreation.")
    if not SILENT: print("If you want to recreate the data file, delete it first.")

# check, if target directory for video frame images exists, and if not, create it
if not os.path.isdir(os.path.join(path, 'PLOTS', 'movies', subpath)):
    os.mkdir(os.path.join(path, 'PLOTS', 'movies', subpath), mode=0o750)



time_1 = time.time()

# the following two sets of data are connected by a set of IDs present in each file

# load meta information about sequences, like protein, country, date, location
allseq = open(os.path.join(path, 'ALLSEQ', subpath + '_allseq.txt')).readlines()

# load data for sequences and sequences themselves, like length of protein, mutations, abundance
mutlist = pd.read_csv(os.path.join(path, 'MUTLISTS', subpath + '_complete_mutation_list.csv'))
if VERBOSE: print('Length of mutation list: ' + str(len(mutlist)))
# define wildtype sequence
refseq = mutlist['Sequence'][0]
if refseq[-1] == '*':
    refseq = refseq[:-1]

# length of reference sequence, as it is
refsize2 = len(refseq)
# mod it if last letter is a star and cut that off
refsize = len(refseq) - 1 if refseq[-1] == '*' else len(refseq)
# get total number of all sequences (not unique, but all)
allseq_number = mutlist['Abundance'].sum()

# get rid of the the wildtype entry and all other NaN entries
# mutlist = mutlist[1:]  # [1:] bc wildtype is nan; former version of this line
mutlist = mutlist[mutlist[mutations_column_name].notna()]
if VERBOSE: print('Length of mutation list without nan in (high qualitiy) mutations: ' + str(len(mutlist)))
# charnge format
mutlist = mutlist.astype({'ID': str})

mutlist.drop_duplicates('ID', inplace=True)
if VERBOSE: print('Length of mutation list with dropped IDs (should be same as above): ' + str(len(mutlist)))


# create an array with all integers from 1 to number of positions
positions = np.arange(1, len(refseq) + 1)

time_2 = time.time()

if VERBOSE: print("Time spent loading data from files:", np.round(time_2 - time_1, 2), "sec\n")

# create a list of all AA at the positions of the original AAS
labels = []
for amino_acid in refseq:
    labels.append(amino_acid)

# create a list of IDs per date
correct_identifiers = []  # write allprot date (if retrievable) and ID in csv
correct_dates = []
correct_location = []
for line in allseq:
    try:
        the_line = line.split('|')
        the_id, the_location, date1 = the_line[:3]
        if CORRECT_DATES:
            # we use corrected dates
            date1 = the_line[-2]
            cor_da = int(date1.replace('-', ''))

            if date1 == 'xxxx-xx-xx' or cor_da < 20191224:
                pass
            else:
                correct_dates.append(cor_da)
                correct_identifiers.append(the_id.split('_')[0][1:])

                if the_location != '':
                    loc_string = the_location.split('/')
                    if len(loc_string) > 1:
                        correct_location.append(loc_string[1])
                    else:
                        correct_location.append('other')
                else:
                    correct_location.append('other')

        else:
            # we don't use corrected dates

            if the_line[2].startswith('EPI'):
                date1 = the_line[3]
            else:
                pass

            if date1[:2] == '20':
                the_year, the_month, the_day = date1.split('-')
                if the_month != '00' and the_month != '0' and the_day != '00' and the_day != '0':

                    if len(the_month) == 2 and len(the_day) == 2:  # if date = 2021-01-02
                        cor_da = int(date1.replace('-', ''))
                        if cor_da >= 20191224:  # the sequences should not occur before the wild type
                            correct_dates.append(cor_da)
                            correct_identifiers.append(the_id.split('_')[0][1:])
                    elif len(the_month) == 1 and len(the_day) == 1:  # if date = 2021-1-2
                        cor_da = int(date1.replace('-', '0'))
                        if cor_da >= 20191224:  # the sequences should not occur before the wild type
                            correct_dates.append(cor_da)
                            correct_identifiers.append(the_id.split('_')[0][1:])
                    elif len(the_month) == 2 and len(the_day) == 1:  # if date = 2021-01-2
                        cor_da = int(date1.replace('-', '')[:-1] + '0' + date1[-1])
                        if cor_da >= 20191224:  # the sequences should not occur before the wild type
                            correct_dates.append(cor_da)
                            correct_identifiers.append(the_id.split('_')[0][1:])
                    elif len(the_month) == 1 and len(the_day) == 2:  # if date = 2021-1-02
                        cor_da = int(date1.replace('-', '')[:-3] + '0' + date1.replace('-', '')[-3:])
                        if cor_da >= 20191224:  # the sequences should not occur before the wild type
                            correct_dates.append(cor_da)
                            correct_identifiers.append(the_id.split('_')[0][1:])
                    if the_location != '':
                        loc_string = the_location.split('/')
                        if len(loc_string) > 1:
                            correct_location.append(loc_string[1])
                        else:
                            correct_location.append('other')
                    else:
                        correct_location.append('other')

    except IndexError:
        pass
    except ValueError:
        pass


# save memory
del allseq 

time_3 = time.time()

if VERBOSE: print("Time spent getting identifiers and dates list:", np.round(time_3 - time_2, 2), "sec\n")

if VERBOSE: print("Number of corect dates:", len(correct_dates))
if VERBOSE: print("Number of corect identifiers:", len(correct_identifiers))

csv_all = pd.DataFrame(list(zip(correct_identifiers, correct_dates, correct_location)),
                       columns=['ID', 'Date', 'Location'])
# print(csv_all)
if VERBOSE: print("Number of combined id-dates", len(csv_all))
csv_dropped_dates = csv_all.drop_duplicates('Date')

if VERBOSE: print("Number of combined id-dates after dropping duplicates", len(csv_dropped_dates))
dates_with_entries = csv_dropped_dates['Date'].to_numpy()

# create a list of dates which have entries
dates_with_entries_numbers = dates_with_entries.astype(int)
sorted_dates_with_entries_numbers = np.sort(dates_with_entries_numbers)
# if VERBOSE: print("Dates with entries:\n", sorted_dates_with_entries_numbers)

first_date = str(sorted_dates_with_entries_numbers[0])
last_date = str(sorted_dates_with_entries_numbers[-1])

# create a list of all dates from first to last entry
sdate = date(int(first_date[:4]),
             int(first_date[4:6]),
             int(first_date[6:]))  # start date
edate = date(int(last_date[:4]),
             int(last_date[4:6]),
             int(last_date[6:]))  # end date

# all_dates = pd.date_range(sdate, edate - timedelta(days=1), freq='d').to_native_types()
all_dates = pd.date_range(sdate, edate, freq='d').to_native_types()
all_dates_numbers = [int(da.replace('-', '')) for da in all_dates]

if VERBOSE: print("Checking correct end date:", all_dates_numbers[-1])

# define a few things necessary in the animation process
number_of_frames = len(all_dates_numbers)

# reset the number of frames, if all are wanted
if USE_FRAMES is None:
    USE_FRAMES = number_of_frames

if use_insertions_and_deletions:
    # list of possible AA letters
    letters_new = ['Y', 'W', 'V', 'T', 'S', 'R', 'Q', 'P', 'N', 'M', 'L', 'K', 'I', 'H', 'G', 'F', 'E', 'D', 'C', 'A',
                   '-', '+', '*'][::-1]

    # list of colors to be assigned to the exchanges
    color_list = ['black', 'dimgray', 'lightgray', 'lightcoral', 'maroon', 'red', 'sienna', 'peru', 'orange', 'gold',
                  'olivedrab', 'darkolivegreen', 'greenyellow', 'darkseagreen', 'mediumaquamarine', 'turquoise',
                  'lightseagreen', 'teal', 'blue', 'mediumpurple', 'mediumorchid', 'hotpink', 'pink']  # [::-1]

else:
    # list of possible AA letters
    letters_new = ['Y', 'W', 'V', 'T', 'S', 'R', 'Q', 'P', 'N', 'M', 'L', 'K', 'I', 'H', 'G', 'F', 'E', 'D', 'C', 'A',
                   '*'][::-1]

    # list of colors to be assigned to the exchanges
    color_list = ['black', 'lightcoral', 'maroon', 'red', 'sienna', 'peru', 'orange', 'gold',
                  'olivedrab', 'darkolivegreen', 'greenyellow', 'darkseagreen', 'mediumaquamarine', 'turquoise',
                  'lightseagreen', 'teal', 'blue', 'mediumpurple', 'mediumorchid', 'hotpink', 'pink']  # [::-1]

interesting_countries_list = ['USA', 'Italy', 'Spain', 'India', 'England', 'Mexico', 'Japan', 'other']

time_4 = time.time()

if VERBOSE: print("Time spent checking and sorting dates:", np.round(time_4 - time_3, 2), "sec\n")

# do this only if necessary
if RECREATE_DATA:

    complete_data = np.zeros((len(all_dates_numbers), refsize2, len(letters_new))).astype(int)
    complete_data, _ = build_first_occurrence_array(csv_all, mutlist, protein_length=len(refseq))

    time_5a = time.time()

    if VERBOSE: print("Time spent sifting through data:", np.round(time_5a - time_4, 2), "sec\n")

    np.save(os.path.join(path, data_temp_file), complete_data)

    time_5 = time.time()

    if VERBOSE: print("Time spent saving data:", np.round(time_5 - time_5a, 2), "sec\n")

else:
    pass
    # just load the data
    complete_data = np.load(os.path.join(path, data_temp_file + ".npy"), allow_pickle=True)

    time_5 = time.time()
    if VERBOSE: print("Time spent loading data:", np.round(time_5 - time_4, 2), "sec\n")

# get data array size
data_shape = np.shape(complete_data)
if RECREATE_DATA:
    if VERBOSE: print("Created data of size:", data_shape)
else:
    if VERBOSE: print("Loaded data of size:", data_shape)

if VERBOSE: print("testing sums of data", np.sum(complete_data, axis=(0, 1)))
if VERBOSE: print(np.sum(complete_data))


if MAKE_PNGS:

    # loop over all dates
    for datum_index, datum in tqdm(enumerate(all_dates_numbers[START_FRAME:START_FRAME + USE_FRAMES])):

        # define everything that changes for each new frame in the animation
        time_a = time.time()

        # print out progress of the loop
        if VERBOSE: print("generating frame number " + str(START_FRAME + datum_index) + "/" + str(data_shape[0]))

        # create the figure for this timeframe via external call

        # prepare and save data to file first

        # get the current date string
        datum = all_dates_numbers[datum_index]

        save_dict = {'data': complete_data[datum_index, :, :],
                     'date': datum,
                     'date_index': datum_index,
                     'positions': positions,
                     'letters': letters_new,
                     'refsize2': refsize2,
                     'colors': color_list,
                     'refsize': refsize,
                     'labels': labels,
                     'path': path,
                     'subpath': subpath,
                     'countries': interesting_countries_list,

                     }

        out_file = open("temp_" + subpath + ".pkl", "wb")
        pickle.dump(save_dict, out_file)
        out_file.close()

        # call external routine via system to avoid memory blowup
        os.system("python create_single_frame.py '" + subpath + "'")

        time_b = time.time()

        this_timing = time_b - time_a
        remaining_frames = (START_FRAME + USE_FRAMES - datum_index)
        if VERBOSE: print(
            "frame number " + str(START_FRAME + datum_index) + " took " + str(np.round(this_timing, 2)) + "sec\n")
        if VERBOSE: print("Based on this, it will take another " + str(
            np.round(this_timing * remaining_frames, 2)) + " seconds to finish rendering the remaining " + str(
            remaining_frames) + " of " + str(USE_FRAMES) + " frames in this movie.")
        if VERBOSE: print("Based on this, it would take " + str(
            np.round(this_timing * data_shape[0], 2)) + " seconds to render all " + str(
            data_shape[0]) + " possible frames.")

    time_7 = time.time()

    if VERBOSE: print("Time spent creating frame pictures:", np.round(time_7 - time_5, 2), "sec\n")

if VERBOSE: print("Timing summary for escape-movie creation for protein " + subpath + ":\n")
if VERBOSE: print("Time spent loading data from files:", np.round(time_2 - time_1, 2), "sec\n")
if VERBOSE: print("Time spent getting identifiers and dates list:", np.round(time_3 - time_2, 2), "sec\n")
if VERBOSE: print("Time spent checking and sorting dates:", np.round(time_4 - time_3, 2), "sec\n")
if RECREATE_DATA:
    if VERBOSE: print("Time spent sifting through data:", np.round(time_5a - time_4, 2), "sec\n")
    if VERBOSE: print("Time spent saving data:", np.round(time_5 - time_5a, 2), "sec\n")
else:
    if VERBOSE: print("Time spent loading data:", np.round(time_5 - time_4, 2), "sec\n")
if MAKE_PNGS:
    if VERBOSE: print("Time spent creating frame pictures:", np.round(time_7 - time_5, 2), "sec\n")

if not SILENT: print("All done!")

#
