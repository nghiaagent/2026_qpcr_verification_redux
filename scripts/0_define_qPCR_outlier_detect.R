# Load packages
library("rlang")
library("tidyverse")

# Functions defined here:
## Strict
### qPCR_clean_strict:
### Leaves no replicates behind in worst case
### qPCR_outlier_strict:
### Marks all replicates as outliers for removal in worse case

# Create function that removes outliers
# in each biological rep using the following method
# If there is 1 observation:
# Remove biological rep
# If there are 2 observations:
# Remove both reps if
# values are too different from each other (selected threshold)
# If there are > 2 observations:
# Remove reps that are too far away from the median Ct (selected threshold)
# If there are > 2 observations, all values far away from each other:
# Remove biological rep

qPCR_clean_strict <- function(.data, Ct, threshold, ...) {
  # convert ... to character vector.
  # This allows for group_by to be passed to this function
  dots <- names(enquos(..., .named = TRUE))

  # calculate the median and count for every sample
  median_Ct <- .data %>%
    group_by(...) %>%
    summarise(
      median_Ct = median({{ Ct }}, na.rm = TRUE),
      count = sum(!is.na({{ Ct }}))
    ) %>%
    ungroup()

  # remove biological reps with only 1 tech reps,
  # test if distance from the median is < threshold
  dev_test <- .data %>%
    full_join(median_Ct, by = dots) %>%
    filter(count > 1) %>%
    mutate(distance_med = abs({{ Ct }} - median_Ct)) %>%
    mutate(
      keep = if_else(
        count == 2,
        if_else(distance_med * 2 < threshold, TRUE, FALSE),
        if_else(distance_med < threshold, TRUE, FALSE)
      )
    )

  # count the # of TRUEs per sample
  count_true <- dev_test %>%
    group_by(...) %>%
    summarise(count_keep = sum(keep, na.rm = TRUE)) %>%
    ungroup()

  # remove all outliers
  clean_data <- dev_test %>%
    full_join(count_true, by = dots) %>%
    filter(
      count_keep > 1,
      keep == TRUE
    ) %>%
    select(-(median_Ct:count_keep))
  clean_data
}

qPCR_outlier_strict <- function(.data, Ct, threshold, ...) {
  # convert ... to character vector.
  # This allows for group_by to be passed to this function
  dots <- names(enquos(..., .named = TRUE))

  # calculate the median and count for every sample
  median_Ct <- .data %>%
    group_by(...) %>%
    summarise(
      median_Ct = median({{ Ct }}, na.rm = TRUE),
      count = sum(!is.na({{ Ct }}))
    ) %>%
    ungroup()

  # remove biological reps with only 1 tech reps,
  # test if distance from the median is < threshold
  dev_test <- .data %>%
    full_join(median_Ct, by = dots) %>%
    filter(count > 1) %>%
    mutate(distance_med = abs({{ Ct }} - median_Ct)) %>%
    mutate(
      keep = if_else(
        count == 2,
        if_else(distance_med * 2 < threshold, TRUE, FALSE),
        if_else(distance_med < threshold, TRUE, FALSE)
      )
    )

  # count the # of TRUEs per sample
  count_true <- dev_test %>%
    group_by(...) %>%
    summarise(count_keep = sum(keep, na.rm = TRUE)) %>%
    ungroup()

  # remove all outliers
  clean_data <- dev_test %>%
    full_join(count_true, by = dots) %>%
    filter(
      count_keep > 1,
      keep == TRUE
    ) %>%
    select(-(median_Ct:count_keep)) %>%
    mutate(is_outlier = FALSE)

  # Join the original data with cleaned data
  # Split dataset into 2 types: outliers and non-outliers, then merge.

  data_outlier <- anti_join(
    .data,
    clean_data
  ) %>%
    mutate(is_outlier = TRUE)

  data_not_outlier <- clean_data %>%
    mutate(is_outlier = FALSE)

  clean_data <- rbind(data_outlier, data_not_outlier)

  clean_data
}
