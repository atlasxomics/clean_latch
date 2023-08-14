import csv
import logging
import math
import numpy as np
import pandas as pd
import statistics
import sys

from typing import Dict, List

logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s",
    level=logging.INFO
)

metrics_output = None
  
def average_duplicates(big_list: List[List[int]]) -> Dict[str, float]:
  """Combine row, col, diag reduction lists; if a barcode occurs in 
  more then one list, returns the average.
  """

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

def filter_sc(singlecell_path: str, position_path: str) -> pd.DataFrame:
  """ Reformat data, remove headers, apply custom column names for
  dataframes, add -1 to positions, remove off tixels.
  """

  singlecell = pd.read_csv(singlecell_path).drop(0, axis=0)
  
  positions = pd.read_csv(position_path, header=None, usecols=[0,1,2,3])
  positions.columns = ['barcode', 'on_off', 'row', 'col']
  positions['barcode'] = positions.loc[:,'barcode'].apply(lambda x: x + "-1")

  merged = pd.merge(positions.astype(object), singlecell.astype(object))
  filtered = merged[merged['on_off'] == 1]

  return filtered

def get_reductions(
    singlecell: pd.DataFrame,
    axis_id: str,
    deviations: int
  ) -> pd.DataFrame:
  """ Return table with barcode|barcode_index|adjust where "adjust"
  is the new value to reduce outlier lanes to; table to be used to
  reduce fragments.tsv  
  """
  global metrics_output
  
  # calculate axis medians
  str_length = singlecell[axis_id].unique().tolist()
  all_indexes = {}
  for i in str_length: 
    pre_list = np.where(singlecell[axis_id] == i)
    indexes = pre_list[0].tolist()
    pre_sort = singlecell.iloc[indexes]['passed_filters'].tolist()
    pre_sort.sort()
    all_indexes[i] = statistics.median(pre_sort)
  
  axisid_info = singlecell[['barcode', axis_id]] 
  apply_medians = singlecell.loc[:,axis_id].apply(lambda x: all_indexes[x])
  singlecell[axis_id] = apply_medians

  mean = statistics.mean(list(all_indexes.values()))
  std = statistics.stdev(list(all_indexes.values()))

  # identify lanes more than x standard deviations above mean
  upper_limit = mean + deviations * std
  
  # Filter singlecell table to only outliers
  singlecell = singlecell[singlecell[axis_id] > upper_limit]
  
  # Store rows/cols being downsampled in global variable
  bad_barcodes = singlecell['barcode'].values.tolist()
  downsampled_elements = set()
  for i in bad_barcodes:
    correct_element = axisid_info.loc[axisid_info['barcode'] == i]
    convert_list = correct_element.values.tolist()
    downsampled_elements.add(str(convert_list[0][1]))
  set_to_string = ', '.join(downsampled_elements)
  metrics_output[axis_id] = set_to_string
  # Add "adjust" column containing value to reduce reads to
  singlecell = singlecell.assign(
    adjust = lambda x: np.ceil(x.passed_filters * (mean / x[axis_id]))
  )
  final = singlecell[['barcode', 'adjust']]
  return final

def get_diag_reductions(
    singlecell: pd.DataFrame,
    deviations: int
  ) -> pd.DataFrame:
  """Return reduction table for diagonal if median of diagonal counts
  an outlier compared to either rows or columns.
  """
  global metrics_output
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
      diag_sc['adjust'] = np.ceil(
        diag_sc['passed_filters'] * (row_mean / diag_mean)
      )
  elif diag_mean > cols_limit:
      diag_sc['adjust'] = np.ceil(
        diag_sc['passed_filters'] * (col_mean / diag_mean)
      )
  else:
      diag_sc['adjust'] = diag_sc['passed_filters']

  if diag_sc.shape[0] > 0: metrics_output['down'] = 'TRUE'
  else: metrics_output['down'] = 'FALSE'
  return diag_sc[['barcode', 'adjust']]

def combine_tables(
    singlecell: pd.DataFrame,
    deviations: int=1
  ) -> pd.DataFrame:

  row_singlecell = singlecell.copy()
  col_singlecell = singlecell.copy()
  dia_singlecell = singlecell.copy()
  row_reductions = get_reductions(row_singlecell, 'row', deviations)
  col_reductions = get_reductions(col_singlecell, 'col', deviations)
  diag_reductions = get_diag_reductions(dia_singlecell, deviations)

  # Concat rows and columns, if a tixel occurs twice, take the average value
  combined_table = average_duplicates([
    row_reductions.values.tolist(),
    col_reductions.values.tolist(),
    diag_reductions.values.tolist()
  ])

  return combined_table

def clean_fragments(
    fragments_path: str,
    r_table: Dict[str, float]
  ) -> pd.DataFrame:
  """Reduce high tixels by randomly downsampling fragments.tsv
  according to reduction table.
  """
  global metrics_output
  logging.info("Loading fragments.tsv")
  fragments = pd.read_csv(
    fragments_path,
    sep='\t',
    header=None,
    comment='#'
  )
  fragments.columns = ['V1', 'V2', 'V3', 'barcode', 'V4']
  metrics_output['og'] = fragments.shape[0]
  frag_copy = fragments.copy()
  outlier_barcodes = list(r_table.keys())

  logging.info("Splitting fragments.tsv")
  normal_frags = fragments[fragments['barcode'].isin(outlier_barcodes) == False]
  outlier_frags = frag_copy[frag_copy['barcode'].isin(outlier_barcodes) == True]

  # To each df in the list, randomly downsample if in reduction list
  logging.info("Downsampling....")
  barcode_groups = outlier_frags.groupby('barcode')
  list_concat = []
  for i in outlier_barcodes:
    outlier = barcode_groups.get_group(i)
    if outlier.shape[0] > int(r_table[i]):
      outlier = outlier.sample(n=math.floor(r_table[i]))
    list_concat.append(outlier)
  
  downsampled_frags = pd.concat(list_concat)
  fragments_cleaned = pd.concat([downsampled_frags, normal_frags])
  metrics_output['final'] = fragments_cleaned.shape[0]
  metrics_output['pct'] = metrics_output['final'] / metrics_output['og']
  return fragments_cleaned

if __name__ == '__main__':

  metrics_output = {}
  run_id = sys.argv[1]
  singlecell_path = sys.argv[2]
  position_path = sys.argv[3]
  fragments_path = sys.argv[4]
  deviations = int(sys.argv[5])

  singlecell = filter_sc(singlecell_path, position_path)
  reduct_dict = combine_tables(singlecell, deviations)
  cleaned = clean_fragments(fragments_path, reduct_dict)

  cleaned.to_csv(
    f'{run_id}_fragments.tsv',
    sep='\t',
    index=False,
    header=False
  )
  
  fields = ['Run_Id', 'Columns downsampled', 'Rows downsampled', 'Diagonal downsampled', 
            'Original fragments', 'Final fragments', 'pct_difng']
  filename = f'{run_id}_cleaning_metrics.csv'
  with open(filename, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(fields)
    writer.writerow(list(metrics_output.values()))
