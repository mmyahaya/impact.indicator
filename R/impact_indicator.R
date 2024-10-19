# install.packages("devtools")
#devtools::install_github("b-cubed-eu/b3gbi")


impact_indicator<-function(cube,
                           status_source,
                           status_data=NULL,
                           region="South Africa",
                           col_scientificName=NULL,
                           col_introductionStatus=NULL,
                           impact_data = eicat_data,
                           col_impact=NULL,
                           col_name=NULL,
                           fun="max"){


  full_species_list<-sort(unique(cube$data$scientificName))
  period=unique(cube$data$year)

  impact_values<-c()
  for(y in period){
    sbs.taxon<-cube$data %>%
      dplyr::filter(year==y) %>%
      dplyr::select(scientificName,cellCode,obs) %>%
      dplyr::group_by(scientificName,cellCode) %>%
      dplyr::summarise(across(obs, sum), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = scientificName, values_from = obs) %>%
      dplyr::arrange(cellCode) %>%
      tibble::column_to_rownames(var = "cellCode")

    species_list<-unique(names(sbs.taxon))

    if(!exists("taxa_status_list")){
      full_species_list<-sort(unique(cube$data$scientificName))
      taxa_status_list<-taxa_status(species_list = full_species_list,
                                    status_source = "WCVP",
                                    status_data=status_data,
                                    region=region,
                                    col_scientificName=col_scientificName,
                                    col_introductionStatus=col_introductionStatus)
    }

    intro.sf<-cube$data %>%
      dplyr::filter(year==y) %>%
      dplyr::left_join(taxa_status_list,
                       by = c("scientificName" = "taxon"))


    status.sf <- intro.sf %>%
      dplyr::group_by(cellCode) %>%
      dplyr::summarise(
        total_intro_obs = sum(obs[introduction_status == "introduced"],
                              na.rm = TRUE),
        total_native_obs = sum(obs[introduction_status == "native"],
                               na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(dplyr::across(c(total_intro_obs, total_native_obs),
                                  ~ ifelse(.==0,NA,.))) %>%
      dplyr::mutate(intro_native=total_intro_obs/total_native_obs) %>%
      dplyr::arrange(cellCode)
    if (!exists("eicat_score_list")){
      eicat_score_list=impact_cat(impact_data = impact_data,
                                  species_list = full_species_list,
                                  col_impact=col_impact,
                                  col_name=col_name,
                                  fun=fun)
    }

    eicat_score<-eicat_score_list[species_list,]

    siteScore<-status.sf$intro_native

    abdundance_impact = sweep(sbs.taxon,2,eicat_score,FUN = "*")
    impactScore = siteScore*abdundance_impact
    impact<-sum(impactScore,na.rm = TRUE)
    impact_values<-rbind(impact_values,c(y,impact))
  }

  impact_values<-as.data.frame(impact_values)
  names(impact_values)<-c("year","value")
  return(impact_values)
}



impact_value<-impact_indicator(cube=taxa_cube,
                           status_source="WCVP",
                           status_data=NULL,
                           region="South Africa",
                           col_scientificName=NULL,
                           col_introductionStatus=NULL,
                           impact_data = eicat_data,
                           col_impact=NULL,
                           col_name=NULL,
                           fun="max")

system.time(
  {
    full_species_list<-sort(unique(cube$data$scientificName))
    period=unique(cube$data$year)

    impact_values<-c()
    species_values<-data.frame()
    for(y in period){
      sbs.taxon<-cube$data %>%
        dplyr::filter(year==y) %>%
        dplyr::select(scientificName,cellCode,obs) %>%
        dplyr::group_by(scientificName,cellCode) %>%
        dplyr::summarise(across(obs, sum), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = scientificName, values_from = obs) %>%
        dplyr::arrange(cellCode) %>%
        tibble::column_to_rownames(var = "cellCode")

      species_list<-unique(names(sbs.taxon))

      if(!exists("taxa_status_list")){
        full_species_list<-sort(unique(cube$data$scientificName))
        taxa_status_list<-taxa_status(species_list = full_species_list,
                                      status_source = "WCVP",
                                      status_data=NULL,
                                      region="South Africa",
                                      col_scientificName=NULL,
                                      col_introductionStatus=NULL)
      }

      intro.sf<-cube$data %>%
        dplyr::filter(year==y) %>%
        dplyr::left_join(taxa_status_list,
                         by = c("scientificName" = "taxon"))


      status.sf <- intro.sf %>%
        dplyr::group_by(cellCode) %>%
        dplyr::summarise(
          total_intro_obs = sum(obs[introduction_status == "introduced"],
                                na.rm = TRUE),
          total_native_obs = sum(obs[introduction_status == "native"],
                                 na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::mutate(dplyr::across(c(total_intro_obs, total_native_obs),
                                    ~ ifelse(.==0,NA,.))) %>%
        dplyr::mutate(intro_native=total_intro_obs/total_native_obs) %>%
        dplyr::arrange(cellCode)
      if (!exists("eicat_score_list")){
        eicat_score_list=impact_cat(impact_data = eicat_data,
                                    species_list = full_species_list,
                                    col_impact=NULL,
                                    col_name=NULL,
                                    fun="max")
      }

      eicat_score<-eicat_score_list[species_list,]

      siteScore<-status.sf$intro_native

      abdundance_impact = sweep(sbs.taxon,2,eicat_score,FUN = "*")
      impactScore = siteScore*abdundance_impact
      impact<-sum(impactScore,na.rm = TRUE)
      impact_values<-rbind(impact_values,c(y,impact))
      
     impact_species <- eicat_score_list %>% 
        na.omit() %>% 
        rownames()
     
     speciesScore<-colSums(impactScore,na.rm = TRUE) %>% 
       as.data.frame() %>% 
       t() %>% 
       as.data.frame() %>% 
       select(any_of(impact_species))
     
     species_values<-bind_rows(species_values,speciesScore)
    
    
      
    }
  }
)



df1 <- tibble(x = 1:2, y = letters[1:2])
df2 <- tibble(x = 4:5, z = 1:2)
bind_rows(df1, df2,as.data.frame(A))

impact_values<-as.data.frame(impact_values)
ggplot() + geom_line(aes(y = value, x = year,colour="red"),
                     data = impact_value, stat="identity")


test %>% as.data.frame() %>% View()
A<-data.frame()
A
names(A)<-names(test)
