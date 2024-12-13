
taxa_Acacia<-readRDS("Data/taxa_Acacia.rds")


countries_sf<-readRDS("Data/countries_shapefile.rds")
SA.sf<-filter(countries_sf,name=="South Africa") %>% select(name,geometry)
acacia_cube<-taxa_cube(taxa=taxa_Acacia,
                       region=SA.sf,
                       res=0.25,
                       first_year=2015)
impact_data<-readRDS("Data/eicat_data.rds")

cube<-acacia_cube$cube
sbs.fun<-function(y){
  sbs.taxon<-cube$data %>%
    filter(year==y) %>%
    dplyr::select(scientificName,cellCode,obs) %>%
    group_by(scientificName,cellCode) %>%
    summarise(across(obs, sum), .groups = "drop") %>%
    pivot_wider(names_from = scientificName, values_from = obs) %>%
    arrange(cellCode) %>%
    column_to_rownames(var = "cellCode")
  sbs.taxon<-as.matrix(sbs.taxon)
  return(sbs.taxon)
}


boot_statistic <- function(data, indices, fun) {
  d <- data[indices,]
  return(fun(d))
}

boot_fun<-function(x){

  sbs.taxon<-x

  species_list<-colnames(sbs.taxon)

  if (!exists("eicat_score_list")){
    eicat_score_list=impact_cat(impact_data = impact_data,
                                species_list = species_list,
                                col_category="impact_category",
                                col_species="scientific_name",
                                col_mechanism="impact_mechanism",
                                trans = 1)

  }


  eicat_score<-eicat_score_list[species_list,"max"]

  #impact score multiply by species by site
  impactScore = sweep(sbs.taxon,2,eicat_score,FUN = "*")

  # Remove rows with all NAs
  impactScore_clean <- impactScore[rowSums(is.na(impactScore)) !=
                                     ncol(impactScore)
                                   , ]

  # Remove columns with all NAs
  if(length(impactScore_clean)!=0){
    impactScore_clean <- impactScore_clean[,
           colSums(is.na(impactScore_clean)) != nrow(impactScore_clean)]

  }

  siteScore<-apply(impactScore_clean,1, function(x) sum(x,
                                                        na.rm = TRUE))

  impact<-sum(siteScore,na.rm = TRUE)/cube$num_cells

  return(impact)
}






fun=boot_fun
samples<-100

period<-acacia_cube$cube$data$year %>% unique()
bootstrap_list <- purrr::map(period,sbs.fun) %>%
  # Perform bootstrapping
  purrr::map(~boot::boot(
    data = .,
    statistic = boot_statistic,
    R = samples,
    fun = fun))
bootstrap_list




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
  
  period <- data_cube_df %>% 
    pull(temporal_col_name) %>% 
    unique()
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
    
    bootstrap_list <- purrr::map(period,sbs.fun) %>%
      # Perform bootstrapping
      lapply(function(x) {
        boot::boot(
          data = x,
          statistic = boot_statistic,
          R = samples,
          fun = fun)
      })
  } else {
    # Define bootstrapping for a difference in a calculated statistic
    boot_statistic_diff <- function(data, ref_data, indices, fun) {
      stat <- fun(data[indices,])
      ref_data <- resize(ref_data,nrow(data))
      ref_stat <- fun(ref_data[indices,])
      
      return(stat - ref_stat)
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
  
  return(bootstrap_list)
}
A<-perform_bootstrap_ts(data_cube_df =  acacia_cube$cube$data,
                        fun = boot_fun,
                        ref_group = NA,
                        seed = 123)



#' Calculate confidence intervals for list of `boot` objects
#'
#' This function calculates confidence intervals for a list of objects of class
#' `"boot"` per year into a dataframe containing all required summaries.
#'
#' @param bootstrap_samples_df A dataframe containing the bootstrap samples.
#' @param grouping_var ...
#' @param type A vector of character strings representing the type of intervals
#' required. The value should be any subset of the values `c("perc", "bca")` or
#' simply `"all"` (default) which will compute both types of intervals.
#' @param conf A scalar or vector containing the confidence level(s) of the
#' required interval(s). Default 0.95.
#' @param aggregate ...
#'
#' @returns The returned value is a dataframe containing the time point,
#' the type of interval (`int_type`), the lower limit of the confidence
#' interval (`ll`), the upper limit of the confidence interval (`ul`), and the
#' confidence level of the intervals (`conf_level`).

get_bootstrap_ci <- function(
    bootstrap_samples_df,
    grouping_var,
    type = c("perc", "bca"),
    conf = 0.95,
    aggregate = TRUE) {
  require("dplyr")
  require("rlang")
  
  # Check if crossv_method is loo or kfold
  type <- tryCatch({
    match.arg(type, c("perc", "bca"))
  }, error = function(e) {
    stop("`type` must be one of 'perc', 'bca'.",
         call. = FALSE)
  })
  
  alpha <- (1 - conf) / 2
  
  if (type == "perc") {
    conf_df <- bootstrap_samples_df %>%
      mutate(
        int_type = type,
        ll = stats::quantile(.data$rep_boot, probs = alpha),
        ul = stats::quantile(.data$rep_boot, probs = 1 - alpha),
        conf_level = conf,
        .by = all_of(grouping_var))
  } else {
    stop("bca not implemented yet")
  }
  
  if (aggregate) {
    conf_df_out <- conf_df %>%
      select(-c("sample", "rep_boot")) %>%
      distinct()
  } else {
    conf_df_out <- conf_df
  }
  
  return(conf_df_out)
}


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
get_bootstrap_ci_old(A)

