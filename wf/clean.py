import pandas as pd
import numpy as np
import sys
import os
import string
import statistics
import gzip

def averageDuplicates(big_list):
  barcodes_match = {}
  final = {}
  ee = []
  for i in big_list:
    ee.extend(i)
    for x in ee:
      if x[0] not in barcodes_match.keys():
        barcodes_match[x[0]] = [x[1]]
      else:
        if x[1] not in barcodes_match[x[0]]:
          barcodes_match[x[0]].append(x[1])
  for i,j in barcodes_match.items():
    mean = statistics.mean(j)
    final[i] = mean
  return final

def filter_sc(singlecell_path, position_path):
  # reformat data, remove headers, apply custom column names for dataframes, add -1 to positions
  # remove off tixels
    singlecell = pd.read_csv(singlecell_path, header=[0,1])
    singlecell = pd.DataFrame(singlecell[:].values)
    
    singlecell.columns = ['barcode', 'total', 'duplicate', 'chimeric', 'unmapped',	'lowmapq', 'mitochondrial',	'nonprimary',	'passed_filters',	'is__cell_barcode',	'excluded_reason',	'TSS_fragments',	'DNase_sensitive_region_fragments',	'enhancer_region_fragments', 'promoter_region_fragments', 'on_target_fragments', 'blacklist_region_fragments', 'peak_region_fragments', 'peak_region_cutsites']
    positions = pd.read_csv(position_path, header=None)
    positions.columns = ['barcode', 'on_off', 'row', 'col', 'y', 'x']
    add_value = positions.loc[:,'barcode'].apply(lambda x: x + "-1")
    positions['barcode'] = add_value
    merged = pd.merge(positions.astype(object), singlecell.astype(object), how ='outer', on = 'barcode')
    final_value = merged[merged['on_off'] == 1]
    return final_value

def get_reductions(singlecell, axis_id, deviations):
  # Return table with barcode|barcode_index|adjust where "adjust" is the new
  # value to reduce outlier lanes to; table to be used to reduce fragments.tsv

  # calculate axis medians
  deviations = int(deviations)
  str_length = singlecell[axis_id].unique().tolist()
  all_indexes = {}
  for i in str_length: 
    pre_list = np.where(singlecell[axis_id] == i)
    indexes = pre_list[0].tolist()
    pre_sort = singlecell.iloc[indexes]['passed_filters'].tolist()
    pre_sort.sort()
    all_indexes[i] = statistics.median(pre_sort)
  
  apply_medians = singlecell.loc[:,axis_id].apply(lambda x: all_indexes[x])
  singlecell[axis_id] = apply_medians

  mean = statistics.mean(list(all_indexes.values()))
  std = statistics.stdev(list(all_indexes.values()))

  # identify lanes more than x standard deviations above mean
  upper_limit = mean + deviations * std
  

  # Filter singlecell table to only outliers
  singlecell = singlecell[singlecell[axis_id] > upper_limit]

  # Add "adjust" column containing value to reduce reads to
  singlecell = singlecell.assign(adjust = lambda x: x.passed_filters * (mean / x[axis_id]))
  final = singlecell[['barcode', 'adjust']].copy()
  return final

def get_diag_reductions(singlecell, deviations):
    # Return reduction table for diagonal if median of diagonal counts an outlier
    # compared to either rows or columns.

    deviations = int(deviations)
    
    str_length_r = singlecell['row'].unique().tolist()
    all_indexes_r = {}
    for i in str_length_r: 
      pre_list = np.where(singlecell['row'] == i)
      indexes = pre_list[0].tolist()
      pre_sort = singlecell.iloc[indexes]['passed_filters'].tolist()
      pre_sort.sort()
      all_indexes_r[i] = statistics.median(pre_sort)
    
    
    str_length_c = singlecell['col'].unique().tolist()
    all_indexes_c = {}
    for i in str_length_c: 
      pre_list = np.where(singlecell['col'] == i)
      indexes = pre_list[0].tolist()
      pre_sort = singlecell.iloc[indexes]['passed_filters'].tolist()
      pre_sort.sort()
      all_indexes_c[i] = statistics.median(pre_sort)
    
    apply_medians_c = singlecell.loc[:,'col'].apply(lambda x: all_indexes_c[x])

    row_mean = statistics.mean(list(all_indexes_r.values()))
    row_std = statistics.stdev(list(all_indexes_r.values()))
    col_mean = statistics.mean(list(all_indexes_c.values()))
    col_std = statistics.stdev(list(all_indexes_c.values()))
    # identify limit more than x standard deviations above mean
    rows_limit = row_mean + deviations * row_std
    cols_limit = col_mean + deviations * col_std

    # create table with only diagonal tixels from singlecell table
    diag_sc = singlecell[singlecell['row'] == singlecell['col']]
    diag_mean = diag_sc['passed_filters'].mean()

    # create 'adjust' column with reads to downsample
    if diag_mean > rows_limit:
        diag_sc['adjust'] = np.ceil(diag_sc['passed_filters'] * (row_mean / diag_mean))
    elif diag_mean > cols_limit:
        diag_sc['adjust'] = np.ceil(diag_sc['passed_filters'] * (col_mean / diag_mean))
    else:
        diag_sc['adjust'] = diag_sc['passed_filters']

    return diag_sc[['barcode', 'adjust']]

def combine_tables(singlecell, deviations=1):
    row_singlecell = singlecell.copy()
    col_singlecell = singlecell.copy()
    dia_singlecell = singlecell.copy()
    row_reductions = get_reductions(row_singlecell, 'row', deviations)
    col_reductions = get_reductions(col_singlecell, 'col', deviations)
    diag_reductions = get_diag_reductions(dia_singlecell, deviations)

    # concat rows and columns
    # If a tixel occurs twice, take the average value
    combined_table = averageDuplicates([row_reductions.values.tolist(), col_reductions.values.tolist(), diag_reductions.values.tolist()])
    
    return combined_table

def clean_fragments(fragments_path, r_table):
    # Reduce high tixels by randomly downsampling fragments.tsv according to
    # reduction table
    skip_headers = [i+1 for i in range(50)]
    print("Loading fragments.tsv")
    # fragments = gzip.open(fragments_path)
    fragments = pd.read_csv(fragments_path, sep='\t', header=skip_headers, skiprows=skip_headers)
    fragments.columns = ['V1', 'V2', 'V3', 'barcode', 'V4']
    frag_copy = fragments.copy()
    outlier_barcodes = list(r_table.keys())
    fragments = fragments[fragments['barcode'].isin(outlier_barcodes) == False]
    frag_copy = frag_copy[frag_copy['barcode'].isin(outlier_barcodes) == True]
    # To each df in the list, randomly downsample if in reduction list
    print("Downsampling....")
    barcode_groups = frag_copy.groupby('barcode')
    list_concat = []
    for i in outlier_barcodes:
      unique_barcode = barcode_groups.get_group(i)
      if unique_barcode.shape[0] > int(r_table[i]):
        unique_barcode = unique_barcode.sample(n= int(r_table[i]))
      list_concat.append(unique_barcode)
    
    downsampled_barcode = pd.concat(list_concat)
    fragments_cleaned = pd.concat([downsampled_barcode, fragments])

    return fragments_cleaned

# do stuff ---------------------------------------------------------------------

args = [i for i in sys.argv]

singlecell = filter_sc(args[1], args[2])
reduct_dict = combine_tables(singlecell, deviations=args[4])
cleaned = clean_fragments(args[3], reduct_dict)

cleaned.to_csv("/root/{}_fragments.tsv".format(args[0]), sep='\t', index=False, header=False)
