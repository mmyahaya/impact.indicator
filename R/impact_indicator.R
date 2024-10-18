# install.packages("devtools")
#devtools::install_github("b-cubed-eu/b3gbi")


sbs.fun<-function(y){
  sbs.taxon<-taxa_cube$data %>%
    dplyr::filter(year==y) %>%
    dplyr::select(scientificName,cellCode,obs) %>%
    dplyr::group_by(scientificName,cellCode) %>%
    dplyr::summarise(across(obs, sum), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = scientificName, values_from = obs) %>%
    dplyr::arrange(cellCode) %>%
    tibble::column_to_rownames(var = "cellCode")
  return(sbs.taxon)
}

full_species_list<-sort(unique(taxa_cube$data$scientificName))
period=unique(taxa_cube$data$year)

sbs.taxon_list<-map(period,sbs.fun)

impact_values<-c()
for(y in period){
  #sbs.taxon_list[[paste0("year","_",y)]]<-sbs.taxon
  sbs.taxon<-sbs.fun(y)

  species_list<-unique(names(sbs.taxon))

  if(!exists("taxa_status_list")){
    full_species_list<-sort(unique(taxa_cube$data$scientificName))
    taxa_status_list<-taxa_status(species_list = full_species_list,
                                    source = "WCVP",
                                    region = "South Africa")
  }

  intro.sf<-taxa_cube$data %>%
    dplyr::filter(year==y) %>%
    dplyr::left_join(taxa_status_list,
              by = c("scientificName" = "taxon"))


  status.sf <- intro.sf %>%
    dplyr::group_by(cellCode) %>%
    dplyr::summarise(
      total_intro_obs = sum(obs[introduction_status == "introduced"], na.rm = TRUE),
      total_native_obs = sum(obs[introduction_status == "native"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(dplyr::across(c(total_intro_obs, total_native_obs),
                                ~ ifelse(.==0,NA,.))) %>%
    dplyr::mutate(intro_native=total_intro_obs/total_native_obs) %>%
    dplyr::arrange(cellCode)
  if (!exists("eicat_score_list")){
    eicat_score_list=impact_cat(data = eicat_data,species_list = species_list,
                                  fun="max")
  }

  eicat_score<-eicat_score_list[species_list,]

  siteScore<-status.sf$intro_native

  abdundance_impact = sweep(sbs.taxon,2,eicat_score,FUN = "*")
  impactScore = siteScore*abdundance_impact
  impact<-sum(impactScore,na.rm = TRUE)
  impact_values<-c(impact_values,impact)
}








