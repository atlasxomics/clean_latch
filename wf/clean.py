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
bad_elements = []
  
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

  singlecell = pd.read_csv(singlecell_path, usecols=[0,8]).drop(0, axis=0)
  
  positions = pd.read_csv(position_path, header=None, usecols=[0,1,2,3])
  positions.columns = ['barcode', 'on_off', 'row', 'col']
  positions['barcode'] = positions.loc[:,'barcode'].apply(lambda x: x + "-1")

  merged = pd.merge(positions.astype(object), singlecell.astype(object))
  filtered = merged[merged['on_off'] == 1]

  return filtered

def get_neighbors(current_value: int, repeat: List[int]) -> List[int]:
  global bad_elements
  
  all_neighbors = {}
  row = current_value[0]
  col = current_value[1]
  
  #right
  if col + 1 < 50 and [row, col + 1] not in bad_elements:
    all_neighbors['r'] = [row, col + 1]
  #left
  if col - 1 >= 0 and [row, col - 1] not in bad_elements:
    all_neighbors['l'] = [row, col - 1]
  #down
  if row + 1 < 50 and [row + 1, col] not in bad_elements:
    all_neighbors['d'] = [row + 1, col]
  #up
  if row - 1 >= 0 and [row - 1, col] not in bad_elements:
    all_neighbors['u'] = [row - 1, col]
  #leftUp
  if row - 1 >= 0 and col - 1 >= 0 and [row - 1, col - 1] not in bad_elements:
    all_neighbors['lu'] = [row - 1, col - 1]
  #leftDown
  if row + 1 < 50 and col - 1 >= 0 and [row + 1, col - 1] not in bad_elements:
    all_neighbors['ld'] = [row + 1, col - 1]
  #rightUp
  if row - 1 >= 0 and col + 1 < 50 and [row - 1, col + 1] not in bad_elements:
    all_neighbors['ru'] = [row - 1, col + 1]
  #rightDown
  if row + 1 < 50 and col + 1 < 50 and [row + 1, col + 1] not in bad_elements:
    all_neighbors['rd'] = [row + 1, col + 1]

  return all_neighbors

def multiple_degree(first_neighbors: List[int], degree: int, current: int) -> List[int]:
  current_neighbors = first_neighbors.copy()
  actual_degree = degree - 1
  for i in first_neighbors:
    for x in range(actual_degree):
      children = get_neighbors(i, current_neighbors)
      current_neighbors += children
  current_neighbors.remove(current)
  return current_neighbors

def neighbors_reductions(
    singlecell: pd.DataFrame,
    outliers: List[int],
    degree: int,
    global_mean: float,
    axis_id: str
  ) -> pd.DataFrame:
  """ Return table with barcode|barcode_index|adjust where "adjust"
  is the new value to reduce outlier lanes to; table to be used to
  reduce fragments.tsv  
  """
  
  singlecell['adjust'] = 0
  for i in outliers:
    current_tixel = singlecell.iloc[i]
    row = current_tixel['row']
    col = current_tixel['col']
    neighbors = get_neighbors([row, col], [])
    # if degree > 1: neighbors += multiple_degree(neighbors, degree, i)        
    on_tixels = []
    for x,j in neighbors.items():
      try:
        current_neighbor = singlecell.loc[(singlecell['row'] == j[0]) & (singlecell['col'] == j[1])]
        if x not in ['lu', 'ld', 'ru', 'rd']:
          on_tixels.append(current_neighbor['passed_filters'].values[0])
        else:
          on_tixels.append(current_neighbor['passed_filters'].values[0] * .7)
      except Exception as e:
        pass
    if len(on_tixels) > 0: 
      mean = statistics.mean(on_tixels)
    else:
      if axis_id != 'diag':
        normalized_median = (global_mean / singlecell.iloc[[i], [1]].values[0][0])
        mean = singlecell.iloc[[i], [4]].values[0][0] * normalized_median
      else:
        mean = singlecell.iloc[[i], [4]].values[0][0] * global_mean
    singlecell.iloc[[i], [5]] = mean
      
  sliced = singlecell[['barcode', 'adjust']]
  filtered = sliced[sliced['adjust'] != 0]
  return filtered

def get_reductions(
    singlecell: pd.DataFrame,
    axis_id: str,
    deviations: int,
    degree: int
  ) -> pd.DataFrame:
  """ Return table with barcode|barcode_index|adjust where "adjust"
  is the new value to reduce outlier lanes to; table to be used to
  reduce fragments.tsv  
  """
  global metrics_output
  global bad_elements
  
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
  singlecell['on_off'] = apply_medians
  og_singlecell = singlecell.copy()

  mean = statistics.mean(list(all_indexes.values()))
  std = statistics.stdev(list(all_indexes.values()))

  # identify lanes more than x standard deviations above mean
  upper_limit = mean + deviations * std
  
  # Filter singlecell table to only outliers
  singlecell = singlecell[singlecell['on_off'] > upper_limit]
  
  # Store rows/cols being downsampled in global variable
  bad_barcodes = singlecell['barcode'].values.tolist()
  downsampled_elements = set()
  
  for i in bad_barcodes:
    correct_element = axisid_info.loc[axisid_info['barcode'] == i]
    convert_list = [[i, j] for i,j in correct_element.values.tolist()]
    downsampled_elements.add(str(convert_list[0][1]))
  set_to_string = ', '.join(downsampled_elements)
  metrics_output[axis_id] = set_to_string
  
  # Add "adjust" column containing value to reduce reads to
  all_elem_ids = []
  for i in downsampled_elements:
    outlier_ids = np.where(og_singlecell[axis_id] == int(i))
    all_elem_ids += outlier_ids[0].tolist()
  for i in all_elem_ids:
    element = og_singlecell.iloc[i]
    row = element['row']
    col = element['col']
    bad_elements.append([row, col])
  updated_singlecell = neighbors_reductions(og_singlecell, all_elem_ids, degree, mean, axis_id)
  
  final = updated_singlecell
  return final

def get_diag_reductions(
    singlecell: pd.DataFrame,
    deviations: int,
    degree: int
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
  all_elem_ids = np.where(singlecell['row'] == singlecell['col'])[0].tolist()
  diag_mean = diag_sc['passed_filters'].mean()
  final_dataFrame = None
  # create 'adjust' column with reads to downsample
  if diag_mean > rows_limit:
    final_dataFrame = neighbors_reductions(singlecell, all_elem_ids, degree, (row_mean/diag_mean), 'diag')
    metrics_output['down'] = 'TRUE'
  elif diag_mean > cols_limit:
    final_dataFrame = neighbors_reductions(singlecell, all_elem_ids, degree, (col_mean/diag_mean), 'diag')
    metrics_output['down'] = 'TRUE'
  else:
      diag_sc['adjust'] = diag_sc['passed_filters']
      final_dataFrame =  diag_sc[['barcode', 'adjust']]
      metrics_output['down'] = 'FALSE'

  return final_dataFrame

def combine_tables(
    singlecell: pd.DataFrame,
    deviations: int=1,
    degree: int=1
  ) -> pd.DataFrame:

  row_singlecell = singlecell.copy()
  col_singlecell = singlecell.copy()
  dia_singlecell = singlecell.copy()
  row_reductions = get_reductions(row_singlecell, 'row', deviations, degree)
  col_reductions = get_reductions(col_singlecell, 'col', deviations, degree)
  diag_reductions = get_diag_reductions(dia_singlecell, deviations, degree)

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
  metrics_output['run_id'] = run_id
  singlecell_path = sys.argv[2]
  position_path = sys.argv[3]
  fragments_path = sys.argv[4]
  deviations = int(sys.argv[5])
  degree = 1

  singlecell = filter_sc(singlecell_path, position_path)
  reduct_dict = combine_tables(singlecell, deviations, degree)
  cleaned = clean_fragments(fragments_path, reduct_dict)

  cleaned.to_csv(
    f'{run_id}_fragments.tsv',
    sep='\t',
    index=False,
    header=False
  )
  
  fields = [
    'Run_Id',
    'Columns downsampled',
    'Rows downsampled',
    'Diagonal downsampled',
    'Original fragments',
    'Final fragments',
    'pct_diff'
  ]
  filename = f'{run_id}_cleaning_metrics.csv'
  with open(filename, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(fields)
    writer.writerow(list(metrics_output.values()))