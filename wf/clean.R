library("dplyr")

# inputs -----------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

run_id          <- args[1]
singlecell_path <- args[2]
position_path   <- args[3]
fragments_path  <- args[4]
deviations      <- as.integer(args[5])

# define functions -------------------------------------------------------------

filter_sc <- function(
  singlecell_path,
  position_path
  ) {
  # Return filtered (on-tissue) singlecell (sc) table with row/column indices

  singlecell <- read.csv(singlecell_path)
  positions  <- read.csv(position_path, header = FALSE)

  positions$V1 <- paste0(positions$V1, "-1") # Make barcode columns match

  merged <- merge(positions, singlecell, by.x = "V1", by.y = "barcode")
  singlecell <- merged[which(merged$V2 == 1), ]
}

get_reductions <- function(
  singlecell,
  axis_id,
  deviations
  ) {
  # Return table with barcode|barcode_index|adjust where "adjust" is the new
  # value to reduce outlier lanes to; table to be used to reduce fragments.tsv

  # calculate axis medians
  medians <- aggregate(
    singlecell$passed_filters,
    by = list(index = singlecell[, axis_id]),
    FUN = median
  )

  # add axis median to singlecell table
  singlecell <- merge(
    singlecell,
    medians,
    by.x = axis_id,
    by.y = "index"
  )

  mean <- mean(medians$x)
  sd <- sd(medians$x)

  # identify lanes more than x standard deviations above mean
  upper_limit <- mean + deviations * sd
  outliers <- subset(medians, x > upper_limit)

  # Filter singlecell table to only outliers
  singlecell <- singlecell[singlecell[, axis_id] %in% outliers$index, ]

  # Add "adjust" column containing value to reduce reads to
  singlecell$adjust <- ceiling(
    singlecell$passed_filters * (mean / singlecell$x)
  )

  # return only barcode|adjust_value
  r_table <- singlecell[, c("V1", "adjust")]
}

get_diag_reductions <- function(singlecell, deviations) {
  # Return reduction table for diagonal if median of diagonal counts an outlier
  # compared to either rows or columns.

  row_medians <- aggregate(
    singlecell$passed_filters,
    by = list(Category = singlecell[, "V3"]),
    FUN = median
  )

  col_medians <- aggregate(
    singlecell$passed_filters,
    by = list(Category = singlecell[, "V4"]),
    FUN = median
  )

  # calculate the mean and standard deviation of the medians
  rows_mean <- mean(row_medians$x)
  cols_mean <- mean(col_medians$x) 
  
  rows_sd <- sd(row_medians$x)
  cols_sd <- sd(col_medians$x)
  
  # identify limit more than x standard deviations above mean 
  rows_limit <- rows_mean + deviations * rows_sd
  cols_limit <- cols_mean + deviations * cols_sd 
  
  # create table with only diagonal tixels from singlecell table
  diag_sc <- subset(singlecell, singlecell$V3 == singlecell$V4)
  diag_mean <- mean(diag_sc$passed_filters)
  
  # create 'adjust' column with reads to downsample
  if (diag_mean > rows_limit) {
    diag_sc$adjust <- ceiling(diag_sc$passed_filters * (rows_mean / diag_mean))
  }
  else if (diagonal > cols_limit) {
    diag_sc$adjust <- ceiling(diag_sc$passed_filters * (cols_mean / diag_mean))
  } 
  else {
    diag_sc$adjust <- diag_sc$passed_filters
  }
  
  return(diag_sc[, c('V1', 'adjust')])
}

combine_tables <- function(
  singlecell,
  deviations = 1,
  row_id = "V3",
  col_id = "V4"
) {

  row_reductions  <- get_reductions(singlecell, row_id, deviations)
  col_reductions  <- get_reductions(singlecell, col_id, deviations)
  diag_reductions <- get_diag_reductions(singlecell, deviations)

  # concat rows and columns
  combined_table <- bind_rows(row_reductions, col_reductions, diag_reductions)

  # If a tixel occurs twice, take the average value
  combined_table <- aggregate(
    combined_table$adjust,
    by = list(combined_table$V1),
    FUN = mean
  )

  colnames(combined_table) <- c("barcode", "adjust")
  return(combined_table)
}

clean_fragments <- function(fragments_path, r_table) {
  # Reduce high tixels by randomly downsampling fragments.tsv according to
  # reduction table

  print("Loading fragments.tsv")
  fragments <- read.table(gzfile(fragments_path), sep = "\t")
  print("fragments.tsv loaded")

  # Split the data frame into a list of data frames by barcode
  print("Splitting fragments")
  df_list <- split(fragments, fragments$V4)
  bc_list <- c(r_table$barcode)

  # To each df in the list, randomly downsample if in reduction list
  print("Downsampling....")
  df_list_cleaned <- lapply(df_list, function(df) {

    barcode <- unique(df$V4)

    if (barcode %in% bc_list) {
      n <- r_table[r_table$barcode == barcode, ]$adjust
      df <- sample_n(df, n)
    } else {
      df <- df
    }
    return(df)
  }
  )

  # Combine the filtered data frames into a single data frame
  fragments_cleaned <- do.call(rbind, df_list_cleaned)
  rownames(fragments_cleaned) <- NULL

  return(fragments_cleaned)
}

# do stuff ---------------------------------------------------------------------

singlecell   <- filter_sc(singlecell_path, position_path)
reduct_table <- combine_tables(singlecell, deviations = deviations)
cleaned <- clean_fragments(fragments_path, reduct_table)

write.table(
  cleaned,
  paste0("/root/", run_id, "_fragments.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
