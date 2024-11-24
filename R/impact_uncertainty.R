sbs.fun<-function(y){
  sbs.taxon<-taxa_cube$data %>%
    filter(year==y) %>%
    dplyr::select(scientificName,cellCode,obs) %>%
    group_by(scientificName,cellCode) %>%
    summarise(across(obs, sum), .groups = "drop") %>%
    pivot_wider(names_from = scientificName, values_from = obs) %>%
    arrange(cellCode) %>%
    column_to_rownames(var = "cellCode")  #%>%
  #mutate(across(all_of(everything()), ~ ifelse(is.na(.),0,.)))

  #colnames(sbs.taxon)<-NULL
  sbs.taxon<-as.matrix(sbs.taxon)
  return(sbs.taxon)
}


my_boot_statistic <- function(data, indices, fun) {
  d <- data[indices]
  #y<-indices
  return(fun(d))
}

my_fun<-function(x){

  species_list<-colnames(x)
  sbs.taxon<-x

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
  
  siteScore<-apply(impactScore_clean,1, function(x) max(x,
                                                        na.rm = TRUE))
  
  impact<-sum(siteScore,na.rm = TRUE)/cube$num_cells

  return(impact)
}

fun=my_fun
samples<-10

period<-acacia_cube$cube$data$year %>% unique()
bootstrap_list <- map(period,sbs.fun) %>%
  # Perform bootstrapping
  purrr::map(~boot::boot(
    data = .,
    statistic = my_boot_statistic,
    R = samples,
    fun = fun))
bootstrap_list


