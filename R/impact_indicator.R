


sbs.fun<-function(y){
  sbs.taxon<-taxon_cube$data %>%
    filter(year==y) %>% 
    dplyr::select(scientificName,cellCode,obs) %>%
    group_by(scientificName,cellCode) %>%
    summarise(across(obs, sum), .groups = "drop") %>%
    pivot_wider(names_from = scientificName, values_from = obs) %>%
    arrange(cellCode) %>% 
    column_to_rownames(var = "cellCode") 
  return(sbs.taxon)
}

full_species_list<-sort(unique(taxon_cube$data$scientificName))
period=

sbs.taxon_list<-map(period,sbs.fun)
for(y in period[-c(1:12)]){
  
  
  sbs.taxon_list[[paste0("year","_",y)]]<-sbs.taxon
  
  species_list<-unique(names(sbs.taxon))
  
  if(!exists("taxon_status_list")){
    full_species_list<-sort(unique(taxon_cube$data$scientificName))
    taxon_status_list<-taxon_status(species_list = full_species_list,
                                    source = "WCVP",
                                    region = "South Africa")
  }
  
  intro.sf<-taxon_cube$data %>%
    filter(year==y) %>%
    left_join(taxa_list_status,
              by = c("scientificName" = "taxon"))
  
  
  status.sf <- intro.sf %>%
    group_by(cellCode) %>%
    summarise(
      total_intro_obs = sum(obs[introduction_status == "introduced"], na.rm = TRUE),
      total_native_obs = sum(obs[introduction_status == "native"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(c(total_intro_obs, total_native_obs), ~ ifelse(.==0,NA,.))) %>%
    mutate(intro_native=total_intro_obs/total_native_obs) %>%
    arrange(cellCode)
  if (!exists("eicat_score_list")){
    eicat_score_list=eicat_impact(eicat_data = eicat_data,species_list = species_list,
                                  fun="max")
  }
  
  eicat_score<-eicat_score_list[species_list,]
  
  siteScore<-status.sf$intro_native
  
  impact_metrics<-list(year=y,sbs.taxon=sbs.taxon,eicat_score=eicat_score,siteScore)
  impact_list[[as.character(y)]]<-impact_metrics
}




species_list<-unique(names(sbs.taxon))


taxon_status_list<-taxon_status(species_list = full_species_list,
                                source = "WCVP",
                                region = "South Africa")


test_data<-data.frame(name=species_list[1:100],
                      status=sample(c("introduced","native"),100,replace = T))

taxon_status_list<-taxon_status(species_list = species_list,
                                source = "manual",
                                status_data = test_data)



intro.sf<-taxon_cube$data %>% 
  filter(year==y) %>% 
  left_join(taxa_list_status,
            by = c("scientificName" = "taxon"))


status.sf <- intro.sf %>%
  group_by(cellCode) %>%
  summarise(
    total_intro_obs = sum(obs[introduction_status == "introduced"], na.rm = TRUE),
    total_native_obs = sum(obs[introduction_status == "native"], na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(across(c(total_intro_obs, total_native_obs), ~ ifelse(.==0,NA,.))) %>% 
  mutate(intro_native=total_intro_obs/total_native_obs) %>% 
  arrange(cellCode)


#intro.sf<-left_join(intro.sf,status.sf, by="cellCode")



eicat_score=eicat_impact(eicat_data = eicat_data,species_list = species_list,
                         fun="max")
eicat<-eicat_score[species_list,]

siteScore<-status.sf$intro_native

abdundance_impact = sweep(sbs.taxon,2,eicat,FUN = "*")
impactScore = siteScore*abdundance_impact
impact<-sum(impactScore,na.rm = TRUE)

# specieImpact<-colSums(impactScore,na.rm = TRUE)
# specieImpact<-specieImpact[specieImpact>0]
# 
# siteImpact<-rowSums(impactScore,na.rm = TRUE)
# siteImpact<-siteImpact[siteImpact>0]





