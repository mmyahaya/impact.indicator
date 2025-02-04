library(dplyr)
library(sf)

source("R/impact_cat.R")
source("R/taxa_cube.R")

taxa_Acacia<-readRDS("Data/taxa_Acacia.rds")


countries_sf<-readRDS("Data/countries_shapefile.rds")
SA.sf<-dplyr::filter(countries_sf,name=="South Africa") %>% select(name,geometry)
acacia_cube<-taxa_cube(taxa=taxa_Acacia,
                       region=SA.sf,
                       res=0.25,
                       first_year=2015)
impact_data<-readRDS("Data/eicat_data.rds")




data_cube_df =  acacia_cube$cube$data

sbs.fun<-function(y){
  sbs.taxon <- data_cube_df %>%
    dplyr::filter(year == y) %>%
    dplyr::mutate(obs=1) %>%
    dplyr::select(scientificName, cellCode, obs) %>%
    # remove duplicates of a species per site
    dplyr::distinct(scientificName, cellCode, .keep_all = TRUE) %>%
    tidyr::pivot_wider(names_from = scientificName, values_from = obs) %>%
    dplyr::arrange(cellCode) %>%
    tibble::column_to_rownames(var = "cellCode")
  return(sbs.taxon)
}


# Resize a data frame to the size of another
resize <- function(A, n) {
  if (nrow(A) > n) {
    A %>%
      slice(1:n)

  } else if (nrow(A) < n) {
    A %>%
      bind_rows(tibble(id = NA, value = NA) %>% slice(rep(1, n - nrow(A))))

  } else {

    A
  }
}


full_species_list <- sort(unique(data_cube_df$scientificName))
period <- data_cube_df %>%
  pull("year") %>%
  unique()

data_list<-purrr::map(period,sbs.fun)

boot_fun<-function(x){

  sbs.taxon<-x

  type="mean cumulative"

  species_list<-colnames(sbs.taxon)

  if (!exists("eicat_score_list")){
    eicat_score_list=impact_cat(impact_data = impact_data,
                                species_list = full_species_list,
                                col_category="impact_category",
                                col_species="scientific_name",
                                col_mechanism="impact_mechanism",
                                trans = 1)

  }

  if (type %in% c("precautionary", "precautionary cumulative")) {
    eicat_score <- eicat_score_list[species_list, "max"]

    # impact score multiply by species by site
    impactScore <- sweep(sbs.taxon, 2, eicat_score, FUN = "*")

    if (type == "precautionary") {
      siteScore <- apply(impactScore, 1, function(x) {
        max(x,
            na.rm = TRUE
        )
      }) %>%
        # suppress warning when -Inf produced  by max() due to site with no impact
        suppressWarnings()

      #d rop -Inf
      siteScore <- siteScore[siteScore!=-Inf]

      num_cells <- length(unique(data_cube_df$cellCode))

      impact<-sum(siteScore,na.rm = TRUE)/num_cells

      return(impact)
    } else {
      # Precautionary cumulative
      siteScore <- apply(impactScore, 1, function(x) {
        sum(x,
            na.rm = TRUE
        )
      })

      num_cells <- length(unique(data_cube_df$cellCode))

      impact<-sum(siteScore,na.rm = TRUE)/num_cells
      return(impact)
    }
  } else if (type %in% c("mean cumulative", "mean")) {
    eicat_score <- eicat_score_list[species_list, "mean"]

    # impact score multiply by species by site
    impactScore <- sweep(sbs.taxon, 2, eicat_score, FUN = "*")


    if (type == "mean cumulative") {
      siteScore <- apply(impactScore, 1, function(x) {
        sum(x,
            na.rm = TRUE
        )
      })
      num_cells <- length(unique(data_cube_df$cellCode))

      impact<-sum(siteScore,na.rm = TRUE)/num_cells
      return(impact)
    } else {
      # mean
      siteScore <- apply(impactScore, 1, function(x) {
        mean(x,
             na.rm = TRUE
        )
      })

      num_cells <- length(unique(data_cube_df$cellCode))

      impact<-sum(siteScore,na.rm = TRUE)/num_cells
      return(impact)
    }
  } else if (type == "cumulative") {
    eicat_score <- eicat_score_list[species_list, "max_mech"]

    # impact score multiply by species by site
    impactScore <- sweep(sbs.taxon, 2, eicat_score, FUN = "*")

    siteScore <- apply(impactScore, 1, function(x) {
      sum(x,
          na.rm = TRUE
      )
    })

    num_cells <- length(unique(data_cube_df$cellCode))

    impact<-sum(siteScore,na.rm = TRUE)/num_cells
    return(impact)
  } else {
    cli::cli_abort(c(
      "{.var type} is not valid",
      "x" = "{.var type} must be from the options provided",
      "See the function desciption or double check the spelling"
    ))
  }



  # eicat_score <- eicat_score_list[species_list,"max"]
  #
  # #impact score multiply by species by site
  # impactScore <- sweep(sbs.taxon,2,eicat_score,FUN = "*")
  #
  #
  #
  # siteScore<-apply(impactScore,1, function(x) sum(x,
  #                                                       na.rm = TRUE))
  #
  # num_cells <- length(unique(data_cube_df$cellCode))
  #
  # impact<-sum(siteScore,na.rm = TRUE)/num_cells
  #
  # return(impact)
}

#' Perform bootstrapping for a calculated statistic over time
#'
#' This function generate `samples` bootstrap replicates of a statistic applied
#' to a data cube per time point (e.g., year, month ...).
#'
#' @param data_cube_df A dataframe containing data in biodiversity data cube
#' format. See `b3gbi::process_cube()`.
#' @param fun A function which when applied to data returns the statistic(s) of
#' interest.
#' @param samples The number of bootstrap replicates. A single positive integer.
#' @param ref_group A string indicating the reference time point to compare the
#' statistic. Default `NA`, no reference time point is used.
#' @param temporal_col_name The temporal column name of `data_cube_df`
#' (e.g., year, month ...) containing time point values. Default `year`.
#' @param seed A positive numeric value setting the seed for random number
#' generation to ensure reproducibility. If `NA` (default), then `set.seed()`
#' is not called at all. If not `NA`, then the random number generator state is
#' reset (to the state before calling this function) upon exiting this function.
#'
#' @returns The returned value is a list of objects of class `"boot"` per time
#' point. See `boot::boot()`.

perform_bootstrap_ts <- function(
    data_cube_df,
    fun,
    samples = 1000,
    ref_group = NA,
    temporal_col_name = "year",
    seed = NA) {
  require("dplyr")
  require("rlang")

  # Check if seed is NA or a number
  stopifnot("`seed` must be a numeric vector of length 1 or NA." =
              (is.numeric(seed) | is.na(seed)) &
              length(seed) == 1)


  # Set seed if provided
  if (!is.na(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
    }
    set.seed(seed)
  }



  if (is.na(ref_group)) {
    # Define bootstrapping for a calculated statistic
    boot_statistic <- function(data, indices, fun) {

      ifelse(is.vector(data),
             d <- data[indices],
             d <- data[indices,])
      return(fun(d))
    }

    bootstrap_list <- data_list %>%
      # Perform bootstrapping
      lapply(function(x) {
        boot::boot(
          data = x,
          statistic = boot_statistic,
          R = samples,
          fun = fun)
      })

    names(bootstrap_list) <- period
  } else {
    # Define bootstrapping for a difference in a calculated statistic
    boot_statistic_diff <- function(data, ref_data, indices, fun) {

      if(is.vector(data)){
        stat <- fun(data[indices])
        ref_stat <- fun(ref_data[indices])

        return(stat - ref_stat)
      } else{
        stat <- fun(data[indices,])
        ref_data <- resize(ref_data,nrow(data))
        ref_stat <- fun(ref_data[indices,])

        return(stat - ref_stat)
      }

    }


    sum_data_list <- purrr::map(period,sbs.fun)
    names(sum_data_list) <- period
    # Perform bootstrapping
    bootstrap_list <- sum_data_list[
      setdiff(names(sum_data_list), as.character(ref_group))
    ] %>%
      lapply(function(x) {
        boot::boot(
          data = x,
          statistic = boot_statistic_diff,
          R = samples,
          fun = fun,
          ref_data = sum_data_list[[as.character(ref_group)]])
      })
  }
  #names(bootstrap_list) <- period
  return(bootstrap_list)
}


boot::boot(
  data = data_list[[1]],
  statistic = boot_statistic,
  R = 1000,
  fun = boot_fun)

A<-perform_bootstrap_ts(data_cube_df =  acacia_cube$cube$data,
                        fun = boot_fun,
                        ref_group = NA,
                        samples = 500,
                        seed = 123)


#' Convert list of `boot` objects to dataframe
#'
#' This function converts a list of objects of class `"boot"` per time point
#' into a dataframe containing all required summaries.
#'
#' @param bootstrap_list A list of objects of class `"boot"` per time point.
#' @param temporal_list_name The temporal list names of `bootstrap_list`
#' (e.g., year, month ...) containing time point values. Default `year`.
#'
#' @returns The returned value is a dataframe containing the bootstrap sample
#' index (`sample`), the time point column (e.g. `year`), the bootstrap estimate
#' of the statistic (`est_boot`), the original sample estimate of the statistic
#' (`est_original`), the standard deviation of the bootstrap replications
#' (`se_boot`), and the bootstrap bias (`bias_boot`).

bootstrap_list_to_df <- function(bootstrap_list, temporal_list_name = "year") {
  require("dplyr")
  require("rlang")

  bootstrap_data_df <- sapply(bootstrap_list, function(df) df$t) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample") %>%
    tidyr::pivot_longer(cols = -sample,
                        values_to = "est_boot",
                        names_to = temporal_list_name) %>%
    dplyr::mutate(sample = as.numeric(sample))

  bootstrap_summaries <- data.frame(
    temp_col = names(bootstrap_list),
    est_original = sapply(bootstrap_list, function(df) df$t0),
    se_boot = sapply(bootstrap_list, function(df) stats::sd(df$t))
  )

  bootstrap_data_full <- bootstrap_data_df %>%
    dplyr::full_join(
      bootstrap_summaries,
      by = dplyr::join_by(!!temporal_list_name == "temp_col")
    ) %>%
    dplyr::arrange(.data$sample, .data[[temporal_list_name]]) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      diff = .data$est_boot - .data$est_original,
      {{ temporal_list_name }} := as.numeric(.data[[temporal_list_name]])) %>%
    dplyr::group_by(.data[[temporal_list_name]]) %>%
    dplyr::mutate(bias_boot = mean(.data$diff)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-"diff")

  return(bootstrap_data_full)
}


bootstrap_data_full <- bootstrap_list_to_df(A)


#' Calculate confidence intervals for list of `boot` objects
#'
#' This function calculates confidence intervals for a list of objects of class
#' `"boot"` per year into a dataframe containing all required summaries.
#'
#' @param bootstrap_list A list of objects of class `"boot"` per year.
#' @param ... Additional argument to be passed to the `boot::boot.ci()`
#' function.
#' @param temporal_list_name The temporal list names of `bootstrap_list`
#' (e.g., year, month ...) containing time point values. Default `year`.
#'
#' @returns The returned value is a dataframe containing the time point,
#' the type of interval (`int_type`), the lower limit of the confidence
#' interval (`ll`), the upper limit of the confidence interval (`ul`), and the
#' confidence level of the intervals (`conf_level`).

get_bootstrap_ci_old <- function(
    bootstrap_list,
    ...,
    temporal_list_name = "year") {
  require("dplyr")
  require("rlang")

  # Calculate nonparametric confidence intervals
  conf_ints <- lapply(bootstrap_list, boot::boot.ci, ...)

  # Remove null values
  conf_ints[sapply(conf_ints, is.null)] <- NULL

  # Exit if there are no values
  if (length(conf_ints) == 0) {
    return(conf_ints)
  }

  # Get interval names
  indices_to_remove <- match(c("R", "t0", "call"), names(conf_ints[[1]]))
  interval_types <- names(conf_ints[[1]])[-indices_to_remove]

  # Get confidence level
  conf_level <- conf_ints[[1]][[interval_types[1]]][1]

  # Summarise for each confidence interval upper and lower limits in dataframes
  out_list <- vector(mode = "list", length = length(interval_types))
  for (i in seq_along(interval_types)) {
    type <- interval_types[i]

    ll <- sapply(conf_ints, function(list) {
      vec <- list[[type]]
      vec[length(vec) - 1]
    })
    ul <- sapply(conf_ints, function(list) {
      vec <- list[[type]]
      vec[length(vec)]
    })

    out_list[[i]] <- data.frame(time_point = as.numeric(names(conf_ints)),
                                int_type = type,
                                ll = ll,
                                ul = ul)
  }

  # Create combined dataframe
  conf_df_out <- do.call(rbind.data.frame, out_list) %>%
    tidyr::complete("time_point" = as.numeric(names(bootstrap_list)),
                    .data$int_type) %>%
    dplyr::arrange(.data$time_point, .data$int_type) %>%
    dplyr::mutate(conf_level = conf_level) %>%
    dplyr::rename({{ temporal_list_name }} := "time_point")
  rownames(conf_df_out) <- NULL

  return(conf_df_out)
}


# Calculate confidence intervals
ci_df <- get_bootstrap_ci_old(
  A)

bootstrap_data_final <- bootstrap_data_full %>%
  full_join(ci_df, by = join_by(year), relationship = "many-to-many")



numbers_data <- acacia_cube$cube$data %>%
  group_by(year) %>%
  summarize(num_occ = sum(.data$obs),
            num_spec = n_distinct(taxonKey),
            .groups = "drop") %>%
  mutate(label_occ = paste("n_occ", num_occ, sep = "="),
         label_spec = paste("n_spec", num_spec, sep = "="),
         label_tot = paste(label_occ, label_spec, sep = "\n"))




bootstrap_data_final %>%
  ggplot(aes(x = year)) +
  geom_boxplot(aes(y = est_boot, group = year)) +
  geom_point(aes(y = est_original), colour = "firebrick", size = 3) +
  geom_label(data = numbers_data, aes(y = 0.81, label = label_tot),
             size = 3, label.padding = unit(0.35, "lines")) +
  labs(y = "impact") +
  scale_y_continuous( breaks = seq(-10, 10, 0.25)) +
  scale_x_continuous(breaks = sort(unique(bootstrap_data_final$year)))




bootstrap_data_final %>%
  ggplot(aes(x = year)) +
  geom_violin(aes(y = est_boot, group = year)) +
  geom_point(aes(y = est_original), colour = "firebrick", size = 3) +
  geom_label(data = numbers_data, aes(y = 2.5, label = label_tot),
             size = 3, label.padding = unit(0.35, "lines")) +
  labs(y = "impact") +
  scale_y_continuous( breaks = seq(-10, 10, 0.25)) +
  scale_x_continuous(breaks = sort(unique(bootstrap_data_final$year)))



bootstrap_data_final_perc<-bootstrap_data_final %>%
  filter(int_type=="percent") %>%
  distinct(year,est_original,ll,ul)





ggplot2::ggplot(data = bootstrap_data_final_perc, aes(x=year)) +
  ggplot2::geom_line(ggplot2::aes(y = est_original, x = year),
                     colour = "red",
                     stat = "identity",
                     linewidth = 1
  ) +
  ggplot2::labs(
    title = "precautionary impact indicator for acacia",
    y = "impact value"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 14))+

  geom_errorbar(aes(ymin = ll, ymax = ul),
                colour = "green",
                alpha = 0.3,
                width = 1,
                linewidth = 1)+
  geom_ribbon(aes(ymin = predict(loess(ll ~ year)),
                  ymax = predict(loess(ul ~ year))),
              alpha=0.3,
              fill="red")


#' Perform bootstrapping over a data cube for a calculated statistic
#'
#' This function generate `samples` bootstrap replicates of a statistic applied
#' to a data cube.
#'
#' @param data_cube A data cube object (class 'processed_cube', see
#' `b3gbi::process_cube()`) or a dataframe (from $data slot of
#' 'processed_cube').
#' @param fun A function which when applied to `data` returns the statistic(s)
#' of interest.
#' @param grouping_var ...
#' @param samples The number of bootstrap replicates. A single positive integer.
#' @param ref_group A string indicating the reference time point to compare the
#' statistic. Default `NA`, no reference time point is used.
#' @param seed A positive numeric value setting the seed for random number
#' generation to ensure reproducibility. If `NA` (default), then `set.seed()`
#' is not called at all. If not `NA`, then the random number generator state is
#' reset (to the state before calling this function) upon exiting this function.
#'
#' @returns The returned value is a list of objects of class `"boot"` per time
#' point. See `boot::boot()`.

bootstrap_cube <- function(
    data_cube,
    fun,
    grouping_var,
    samples = 1000,
    ref_group = NA,
    seed = NA) {
  require("dplyr")
  require("rlang")

  # Check if seed is NA or a number
  stopifnot("`seed` must be a numeric vector of length 1 or NA." =
              (is.numeric(seed) | is.na(seed)) &
              length(seed) == 1)

  # Set seed if provided
  if (!is.na(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
    }
    set.seed(seed)
  }

  if (inherits(data_cube, "processed_cube")) {
    # Generate bootstrap replicates
    resample_df <- modelr::bootstrap(data_cube$data, samples, id = "id")

    # Function for bootstrapping
    bootstrap_resample <- function(x, fun) {
      resample_obj <- x$strap[[1]]
      indices <- as.integer(resample_obj)
      data <- resample_obj$data[indices, ]

      data_cube_copy <- data_cube
      data_cube_copy$data <- data

      fun(data_cube_copy)$data %>%
        mutate(sample = as.integer(x$id))
    }
  } else {
    # Generate bootstrap replicates
    resample_df <- modelr::bootstrap(data_cube, samples, id = "id")

    # Function for bootstrapping
    bootstrap_resample <- function(x, fun) {
      resample_obj <- x$strap[[1]]
      indices <- as.integer(resample_obj)
      data <- resample_obj$data[indices, ]

      fun(data) %>%
        mutate(sample = as.integer(x$id))
    }
  }

  # Perform bootstrapping
  bootstrap_samples_list_raw <- resample_df %>%
    split(seq_len(nrow(resample_df))) %>%
    purrr::map(bootstrap_resample, fun = fun, .progress = TRUE)

  if (!is.na(ref_group)) {
    # Calculate true statistic
    if (inherits(data_cube, "processed_cube")) {
      t0_full <- fun(data_cube)$data
    } else {
      t0_full <- fun(data_cube)
    }

    ref_val <- t0_full %>%
      filter(.data[[grouping_var]] == !!ref_group) %>%
      pull(.data$diversity_val)

    t0 <- t0_full %>%
      filter(.data[[grouping_var]] != !!ref_group) %>%
      mutate(diversity_val = .data$diversity_val - ref_val)

    # Get bootstrap samples as a list
    bootstrap_samples_list <- lapply(bootstrap_samples_list_raw, function(df) {
      ref_val <- df %>%
        filter(.data[[grouping_var]] == !!ref_group) %>%
        pull(.data$diversity_val)

      df %>%
        filter(.data[[grouping_var]] != !!ref_group) %>%
        mutate(diversity_val = .data$diversity_val - ref_val)
    })
  } else {
    # Calculate true statistic
    if (inherits(data_cube, "processed_cube")) {
      t0 <- fun(data_cube)$data
    } else {
      t0 <- fun(data_cube)
    }

    # Get bootstrap samples as a list
    bootstrap_samples_list <- bootstrap_samples_list_raw
  }

  # Summarise in dataframe
  bootstrap_samples_df <- bootstrap_samples_list %>%
    dplyr::bind_rows() %>%
    dplyr::rename("rep_boot" = "diversity_val") %>%
    dplyr::left_join(t0, by = grouping_var) %>%
    dplyr::rename("est_original" = "diversity_val") %>%
    dplyr::mutate(
      est_boot = mean(.data$rep_boot),
      se_boot = stats::sd(.data$rep_boot),
      .by = all_of(grouping_var)) %>%
    dplyr::mutate(bias_boot = .data$est_boot - .data$est_original) %>%
    dplyr::arrange(.data[[grouping_var]]) %>%
    dplyr::select("sample", all_of(grouping_var), "est_original",
                  dplyr::everything())

  return(bootstrap_samples_df)
}



bootstrap_insect_data <- bootstrap_cube(
  data_cube = acacia_cube$cube$data,
  fun = boot_fun,
  grouping_var = "year",
  samples = 1000,
  seed = 123)
