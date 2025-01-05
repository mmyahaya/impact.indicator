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


  eicat_score <- eicat_score_list[species_list,"max"]

  #impact score multiply by species by site
  impactScore <- sweep(sbs.taxon,2,eicat_score,FUN = "*")



  siteScore<-apply(impactScore,1, function(x) sum(x,
                                                        na.rm = TRUE))
  
  num_cells <- length(unique(data_cube_df$cellCode))

  impact<-sum(siteScore,na.rm = TRUE)/num_cells

  return(impact)
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
B<-perform_bootstrap_ts(data_cube_df =  acacia_cube$cube$data,
                        fun = boot_fun,
                        ref_group = NA,
                        samples = 100,
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



bootstrap_list_to_df(A)




