#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 10:15:49 2022

Input file for mpro_escape_movie.py

to be expanded

@author: andreaskrassnigg
"""


all_protein_names = ['E', 'M', 'N', 'NS3', 'NS6', 'NS7a', 'NS7b', 'NS8', 'NS9b',  # 0 - 8
                     'NS9c', 'NS10', 'NSP1', 'NSP2', 'NSP3', 'NSP4', 'NSP5', 'NSP6',      # 9 - 15
                     'NSP7', 'NSP8', 'NSP9', 'NSP10', 'NSP11', 'NSP12', 'NSP13',  # 16 - 22
                     'NSP14', 'NSP15', 'NSP16', 'Spike']                          # 23 - 26


# run all proteins
run_protein_names = all_protein_names     

# specify a sublist to be run at once
# run_protein_names = ['NSP3', 'NSP12', 'NSP13', 'Spike']    


all_stages = ['escape_movie', 'time_evolution_analysis']

# run all stages
# run_stages = all_stages     

run_stages = ['escape_movie']

# pick one of these and put it in the dict

# the dict should contain all meaninful switches for the code over time
the_inputs = {'all_protein_names' : all_protein_names,   # to be able to run a set of proteins specified here
              'run_protein_names' : run_protein_names,   # this is the list to be run
              'run_stages'        : run_stages,  # which stages to include in an entire-workflow run
              'verbose'           : True,   # whether or not to generate optional outputs
              'silent'            : False,  # shuts off all outputs completely
              'protein'           : 'NSP5',
              'path'              : './Final_data_analysis',  
              'drop_last_date'    : True,
              'use_ins_dels'      : False,
              'mutation_column'   : 'high_quality_exchanges',   # 'high_quality_exchanges', 'Mutations'
              'correct_dates'     : False,    # use corrected dates from allseq file
              'make_pngs'         : False,  # make a png image for each date/frame
              'start_frame'       : 0,     # specify frames to produce (more useful for debugging)
              'use_frames'        : None,  # None means use all. Specify a number like 875, 829 or something else, if needed
              'num_processes'     : 8,     # if N > 1, lets N processes run in parallel when creating the images.
              'Compression'       : True   # whether or not to compress the data to save RAM.
              }