# install.packages("devtools")
#devtools::install_github("b-cubed-eu/b3gbi")


impact_indicator<-function(cube,
                           impact_data = NULL,
                           col_impact=NULL,
                           col_name=NULL,
                           type=NULL,
                           coords=NULL){


  full_species_list<-sort(unique(cube$data$scientificName))
  period<-unique(cube$data$year)
  impact_values<-c()
  species_values<-data.frame()
  for(y in period){
    sbs.taxon<-cube$data %>%
      dplyr::filter(year==y) %>%
      dplyr::select(scientificName,cellCode,obs) %>%
      dplyr::group_by(scientificName,cellCode) %>%
      dplyr::summarise(across(obs, first), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = scientificName, values_from = obs) %>%
      dplyr::arrange(cellCode) %>%
      tibble::column_to_rownames(var = "cellCode")

    species_list<-unique(names(sbs.taxon))

    if (!exists("eicat_score_list")){
      eicat_score_list=impact_cat(impact_data = impact_data,
                                  species_list = full_species_list,
                                  col_impact=col_impact,
                                  col_name=col_name)

      impact_species <- eicat_score_list %>%
        na.omit() %>%
        rownames()
    }

    if(type %in% c("precautionary","precautionary cumulative")){

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
      # if any species has impact
      if(!is.null(dim(impactScore_clean))){
        if(type=="precautionary"){
          siteScore<-apply(impactScore_clean,1, function(x) max(x,
                                                                na.rm = TRUE))

          impact<-sum(siteScore,na.rm = TRUE)/cube$num_cells
          impact_values<-rbind(impact_values,c(y,impact))
          
          siteScore<- siteScore%>%
            as.data.frame() %>% 
            tibble::rownames_to_column(var = "siteID")
          
          names(siteScore)[2]<-as.character(y)
          
        } else {

          #Precautionary cumulative
          siteScore<-apply(impactScore_clean,1, function(x) sum(x,
                                                                na.rm = TRUE))


          impact<-sum(siteScore,na.rm = TRUE)/cube$num_cells
          impact_values<-rbind(impact_values,c(y,impact))
          
          siteScore<- siteScore%>%
            as.data.frame() %>% 
            tibble::rownames_to_column(var = "siteID")
          
          names(siteScore)[2]<-as.character(y)
          
          
        }
      }
      else { # return NA if no species has impact
        siteScore<-data.frame("siteID"=NA,"year"=NA)
        names(siteScore)[2]<-as.character(y)
        impact<-NA
        impact_values<-rbind(impact_values,c(y,impact))
      }


    } else if (type %in% c("mean cumulative","mean")){
      eicat_score<-eicat_score_list[species_list,"mean"]

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
      # if any species has impact
      if(!is.null(dim(impactScore_clean))){
        if(type=="mean cumulative"){
          siteScore<-apply(impactScore_clean,1, function(x) sum(x,
                                                                na.rm = TRUE))

          impact<-sum(siteScore,na.rm = TRUE)/cube$num_cells
          impact_values<-rbind(impact_values,c(y,impact))
          
          siteScore<- siteScore%>%
            as.data.frame() %>% 
            tibble::rownames_to_column(var = "siteID")
          
          names(siteScore)[2]<-as.character(y)
          
        } else {

          #mean
          siteScore<-apply(impactScore_clean,1, function(x) mean(x,
                                                                na.rm = TRUE))


          impact<-sum(siteScore,na.rm = TRUE)/cube$num_cells
          impact_values<-rbind(impact_values,c(y,impact))
          siteScore<- siteScore%>%
            as.data.frame() %>% 
            tibble::rownames_to_column(var = "siteID")
          
          names(siteScore)[2]<-as.character(y)
        }
      }
      else { # return NA is no species has impact
        siteScore<-data.frame("siteID"=NA,"year"=NA)
        names(siteScore)[2]<-as.character(y)
        impact<-NA
        impact_values<-rbind(impact_values,c(y,impact))
      }




    } else {
      eicat_score<-eicat_score_list[species_list,"max_mech"]

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

      # if any species has impact
      if(!is.null(dim(impactScore_clean))){
        if(type=="cumulative"){
          siteScore<-apply(impactScore_clean,1, function(x) sum(x,
                                                                na.rm = TRUE))

          impact<-sum(siteScore,na.rm = TRUE)/cube$num_cells
          impact_values<-rbind(impact_values,c(y,impact))
          siteScore<- siteScore%>%
            as.data.frame() %>% 
            tibble::rownames_to_column(var = "siteID")
          
          names(siteScore)[2]<-as.character(y)
          
        }
      }
      else { # return NA is no species has impact
        siteScore<-data.frame("siteID"=NA,"year"=NA)
        names(siteScore)[2]<-as.character(y)
        impact<-NA
        impact_values<-rbind(impact_values,c(y,impact))
      }


    }
    
    coords<-left_join(coords,siteScore,by="siteID")
    # Species impact
    speciesScore<-colSums(impactScore,na.rm = TRUE) %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame() %>%
      select(any_of(impact_species))

    species_values<-bind_rows(species_values,speciesScore)



  }

  impact_values<-as.data.frame(impact_values)
  names(impact_values)<-c("year","value")
  rownames(species_values)<-as.character(period)
  return(list("impact_values"=impact_values,"species_values"=species_values,
              "sitedf"=coords))
}



impact_value<-impact_indicator(cube=acacia_cube$cube,
                           impact_data = eicat_data,
                           col_impact=NULL,
                           col_name=NULL,
                           type = "mean",
                           coords=acacia_cube$coords)

ggplot() + geom_line(aes(y = value, x = year),colour="red",
                     data = impact_value$impact_values, stat="identity")+
  labs(
    title = "Impact risk",
    y = "sum of risk map value"
  )+theme(text=element_text(size=20))

df<-impact_value$species_values %>%
  rownames_to_column("year") %>%
  mutate(year=as.numeric(year)) %>%
  gather(-year,key = "Alien_species", value = "impact_score")


ggplot(df, aes(x = year, y = impact_score)) +
  geom_line(aes(color = Alien_species))
length(unique(taxa_Acacia$species))

sitedf<-impact_value$sitedf
sitedf %>% 
  gather(-c(siteID,X,Y),key="year",value="impact") %>% 
  na.omit() %>% 
  #filter(year>=2021) %>% 
  ggplot() +
  geom_tile(
            aes(x=X,y=Y,fill=impact),color="black")+
  geom_sf(data = SA.sf, fill = NA, color = "black", alpha = 0.5)+
  scale_fill_gradient2(low = "forestgreen",
                       mid = "yellow",
                       high = "red")+
  theme_minimal() +
  facet_wrap(~year)
  
