#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to help in data preprocessing.

@author: tobias.schopper
"""
import numpy as np
import pandas as pd
from datetime import datetime
from tqdm import tqdm
from multiprocessing import Pool
from tqdm.contrib.concurrent import process_map, thread_map
def build_first_occurrence_array(csv_all: pd.DataFrame, mutlist: pd.DataFrame, protein_length: int = 306):
    csv_all.sort_values(by="Date")
    mutations = np.ones((protein_length, 21)) * 23450420  # initialize mutations so that they have not yet occurred.
    # here we will be storing only the first date of when a mutation occured.
    letters_new = ['Y', 'W', 'V', 'T', 'S', 'R', 'Q', 'P', 'N', 'M', 'L', 'K', 'I', 'H', 'G', 'F', 'E', 'D', 'C', 'A',
                   '*'][::-1]

    mapping = {l: i for i, l in enumerate(letters_new)}  # maps letters to numbers
    back_mapping = {i: l for i, l in enumerate(letters_new)}  # maps numbers back to letters
    ID_mut_mapping = {ID: muts for ID, muts in
                      mutlist[["ID", "high_quality_exchanges"]].__array__()}  # maps IDs to mutations

    unique_dates = set(csv_all["Date"])  # give all unique dates
    max_date, min_date = max(unique_dates), min(unique_dates)  # get maximum and minimum dates

    first_date = datetime.strptime(str(min_date), "%Y%m%d")  # convert dates to datetime
    last_date = datetime.strptime(str(max_date), "%Y%m%d")

    all_dates = pd.date_range(first_date, last_date, freq='d').to_native_types()  # get daterange for final array
    all_dates_numbers = [int(da.replace('-', '')) for da in all_dates]

    for i, entry in tqdm(enumerate(csv_all.iloc)):  # for element in dataframe
        ID = entry["ID"]
        date = entry["Date"]
        if ID in ID_mut_mapping.keys():  # check if ID is actually in there
            muts = ID_mut_mapping[ID]  # Get mutation from the ID-mutation map
            for mut in muts.split(", "):  # split properly in case of multiple mutations
                mut = mut.rstrip()
                pos = int(mut[1:-1])
                mut = mapping[mut[-1]]
                if mutations[pos - 1][mut] > date:  # if date in array after current date, replace by earliest
                    mutations[pos - 1][mut] = date  # "mutations" was initialized so that every date would be the last.

    muts_by_day = []  # expand to mutations per (day, position, letter)
    for date in all_dates_numbers:
        muts_by_day.append((mutations < date)[None, ...])  # muts_by_day at position [date] is now a boolean
    muts_by_day = np.vstack(muts_by_day)  # array with "true" at (date,pos,letter) if it has occured.
    return muts_by_day, back_mapping  # return the (date,pos,letter) boolean array, as well as the
    # dictionary for mapping indices back to numbers.#
