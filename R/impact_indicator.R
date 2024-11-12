#' Compute impact indicator
#'
#' @param cube The list containing data cube of class `sim_cube` from 
#' `b3gbi::process_cube()`.
#' @param impact_data The dataframe of species impact which contains columns of category,
#'  species and mechanism.
#' @param col_category The name of the column containing the impact categories.
#' The first two letters each categories must be an EICAT short names 
#' (e.g "MC - Minimal concern")
#' @param col_species The name of the column containing species names 
#' @param col_mechanism The name of the column containing mechanisms of impact
#' @param trans Numeric. The type of transformation to convert the EICAT categories to 
#' numerical values. 1 converts ("MC", "MN", "MO", "MR", "MV") to (0,1,2,3,4)
#' 2 converts ("MC", "MN", "MO", "MR", "MV") to (1,2,3,4,5) and 
#' 3 converts ("MC", "MN", "MO", "MR", "MV") to (1,10,100,1000,10000)
#' @param type The type indicators based on the aggregation of within and 
#' across species in a site. The type can be precautionary, precautionary cumulative,
#' mean, mean cumulative or cumulative. 
#' @param coords The dataframe containing coordinates of the sites of the region. 
#'
#' @return A list containing dataframes impact indicator, species impact indicator
#' and sites impact indicator
#' @export
#'
#' @examples
#' taxa_Acacia<-readRDS("Data/taxa_Acacia.rds")
#'  acacia_cube<-taxa_cube(taxa=taxa_Acacia,
#'                       region=SA.sf,
#'                       res=0.25,
#'                       first_year=2010)
#' eicat_data<-readRDS("Data/eicat_data.rds")
#' impact_value<-impact_indicator(cube=acacia_cube$cube,
#'                               impact_data = eicat_data,
#'                               col_category="impact_category",
#'                               col_species="scientific_name",
#'                                col_mechanism="impact_mechanism",
#'                               trans=1,
#'                               type = "mean cumulative")

impact_indicator<-function(cube,
                           impact_data = NULL,
                           col_category=NULL,
                           col_species=NULL,
                           col_mechanism=NULL,
                           trans=1,
                           type=NULL){


  full_species_list<-sort(unique(cube$data$scientificName))
  
  period<-unique(cube$data$year)
  
  #create empty vector and dataframe for impact indicators and species impact
  #indicator
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
                                  col_category=col_category,
                                  col_species=col_species,
                                  col_mechanism = col_mechanism,
                                  trans = trans)

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

        } else {

          #Precautionary cumulative
          siteScore<-apply(impactScore_clean,1, function(x) sum(x,
                                                                na.rm = TRUE))


          impact<-sum(siteScore,na.rm = TRUE)/cube$num_cells
          impact_values<-rbind(impact_values,c(y,impact))

        }
      }
      else { # return NA if no species has impact

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

        } else {

          #mean
          siteScore<-apply(impactScore_clean,1, function(x) mean(x,
                                                                na.rm = TRUE))


          impact<-sum(siteScore,na.rm = TRUE)/cube$num_cells
          impact_values<-rbind(impact_values,c(y,impact))
        }
      }
      else { # return NA is no species has impact
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

        }
      }
      else { # return NA is no species has impact
        impact<-NA
        impact_values<-rbind(impact_values,c(y,impact))
      }

    }

  }

  impact_values<-as.data.frame(impact_values)
  names(impact_values)<-c("year","value")
  return(impact_values)
}



#' Compute site impact indicator
#'
#' @param cube The list containing data cube of class `sim_cube` from 
#' `b3gbi::process_cube()`.
#' @param impact_data The dataframe of species impact which contains columns of category,
#'  species and mechanism.
#' @param col_category The name of the column containing the impact categories.
#' The first two letters each categories must be an EICAT short names 
#' (e.g "MC - Minimal concern")
#' @param col_species The name of the column containing species names 
#' @param col_mechanism The name of the column containing mechanisms of impact
#' @param trans Numeric. The type of transformation to convert the EICAT categories to 
#' numerical values. 1 converts ("MC", "MN", "MO", "MR", "MV") to (0,1,2,3,4)
#' 2 converts ("MC", "MN", "MO", "MR", "MV") to (1,2,3,4,5) and 
#' 3 converts ("MC", "MN", "MO", "MR", "MV") to (1,10,100,1000,10000)
#' @param type The type indicators based on the aggregation of within and 
#' across species in a site. The type can be precautionary, precautionary cumulative,
#' mean, mean cumulative or cumulative. 
#' @param coords The dataframe containing coordinates of the sites of the region. 
#'
#' @return The dataframe of impact indicator per sites
#' @export
#'
#' @examples
#' 
#' taxa_Acacia<-readRDS("Data/taxa_Acacia.rds")
#'  acacia_cube<-taxa_cube(taxa=taxa_Acacia,
#'                       region=SA.sf,
#'                       res=0.25,
#'                       first_year=2010)
#' eicat_data<-readRDS("Data/eicat_data.rds")
#' siteImpact<-site_impact(cube=acacia_cube$cube,
#'                        impact_data = eicat_data,
#'                        col_category="impact_category",
#'                        col_species="scientific_name",
#'                        col_mechanism="impact_mechanism",
#'                        trans=1,
#'                       type = "precautionary cumulative",
#'                        coords=acacia_cube$coords)
#' 
#' 
#'
site_impact<-function(cube,
                      impact_data = NULL,
                      col_category=NULL,
                      col_species=NULL,
                      col_mechanism=NULL,
                      trans=1,
                      type=NULL,
                      coords=NULL){
  
  
  full_species_list<-sort(unique(cube$data$scientificName))
  
  period<-unique(cube$data$year)
  
  
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
                                  col_category=col_category,
                                  col_species=col_species,
                                  col_mechanism = col_mechanism,
                                  trans = trans)
      
    }
    
    if(type %in% c("precautionary","precautionary cumulative")){
      
      eicat_score<-eicat_score_list[species_list,"max"]
      
      #impact score multiply by species by site
      impactScore = sweep(sbs.taxon,2,eicat_score,FUN = "*")
      
      if(type=="precautionary"){
        siteScore<-apply(impactScore,1, function(x) max(x,
                                                        na.rm = TRUE)) %>% 
          suppressWarnings()
      } else { # elsecompute precautionary cumulative
        siteScore<-apply(impactScore,1, function(x) sum(x,
                                                        na.rm = TRUE))
      }
      
    } else if (type %in% c("mean cumulative","mean")){
      eicat_score<-eicat_score_list[species_list,"mean"]
      
      #impact score multiply by species by site
      impactScore = sweep(sbs.taxon,2,eicat_score,FUN = "*")
      
      if(type=="mean"){
        siteScore<-apply(impactScore,1, function(x) mean(x,
                                                         na.rm = TRUE))
      } else { # else compute mean cumulative
        siteScore<-apply(impactScore,1, function(x) sum(x,
                                                        na.rm = TRUE))
      }
      
      
    } else if(type=="cumulative") {
      eicat_score<-eicat_score_list[species_list,"max_mech"]
      
      #impact score multiply by species by site
      impactScore = sweep(sbs.taxon,2,eicat_score,FUN = "*")
      
      siteScore<-apply(impactScore,1, function(x) sum(x,
                                                      na.rm = TRUE))
      
    }
    else(stop("'type' is not valid. Make sure it is from the options provided"))
    
    # convert siteScore to dataframe
    siteScore<- siteScore%>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "siteID")
    
    names(siteScore)[2]<-as.character(y)
    
    coords<-left_join(coords,siteScore,by="siteID")
    
  }
  
  # remove -Inf produced by max(.,na.rm=TRUE)
  # remove NaN produced by mean(.,na.rm=TRUE)
  # remove 0 produced by sum(.,na.rm=TRUE)
  
  coords <- coords %>% 
    mutate(across(all_of(everything()), ~ ifelse(.==-Inf|is.nan(.)|.==0,NA,.)))
  
  return(coords)
}



# siteImpact%>% 
#   gather(-c(siteID,X,Y),key="year",value="impact") %>% 
#   na.omit() %>% 
#   filter(year>=2021) %>% 
#   ggplot() +
#   geom_tile(
#     aes(x=X,y=Y,fill=impact),color="black")+
#   geom_sf(data = SA.sf, fill = NA, color = "black", alpha = 0.5)+
#   scale_fill_gradient(low = "yellow",
#                       high = "red")+
#   theme_minimal() +
#   labs(
#     title = "impact risk map",
#     y = "Latitude", x="Longitude"
#   )+
#   theme(text=element_text(size=14))+
#   facet_wrap(~year)



#' Compute species impact indicator
#'
#' @param cube The list containing data cube of class `sim_cube` from 
#' `b3gbi::process_cube()`.
#' @param impact_data The dataframe of species impact which contains columns of category,
#'  species and mechanism.
#' @param col_category The name of the column containing the impact categories.
#' The first two letters each categories must be an EICAT short names 
#' (e.g "MC - Minimal concern")
#' @param col_species The name of the column containing species names 
#' @param col_mechanism The name of the column containing mechanisms of impact
#' @param trans Numeric. The type of transformation to convert the EICAT categories to 
#' numerical values. 1 converts ("MC", "MN", "MO", "MR", "MV") to (0,1,2,3,4)
#' 2 converts ("MC", "MN", "MO", "MR", "MV") to (1,2,3,4,5) and 
#' 3 converts ("MC", "MN", "MO", "MR", "MV") to (1,10,100,1000,10000)
#' @param type The type indicators based on the aggregation of within and 
#' across species in a site. The type can be precautionary, precautionary cumulative,
#' mean, mean cumulative or cumulative. 
#'
#' @return A dataframe of impact indicator per species
#' @export
#'
#' @examples
#' 
#' taxa_Acacia<-readRDS("Data/taxa_Acacia.rds")
#'  acacia_cube<-taxa_cube(taxa=taxa_Acacia,
#'                       region=SA.sf,
#'                       res=0.25,
#'                       first_year=2010)
#' eicat_data<-readRDS("Data/eicat_data.rds")
#' 
#' speciesImpact<-species_impact(cube=acacia_cube$cube,
#'                         impact_data = eicat_data,
#'                         col_category="impact_category",
#'                         col_species="scientific_name",
#'                         col_mechanism="impact_mechanism",
#'                         trans=1,
#'                         type = "mean")
species_impact<-function(cube,
                         impact_data = NULL,
                         col_category=NULL,
                         col_species=NULL,
                         col_mechanism=NULL,
                         trans=1,
                         type=NULL){
  
  
  full_species_list<-sort(unique(cube$data$scientificName))
  
  period<-unique(cube$data$year)
  
  #create empty vector for species impact
  speciesImpact<-c()
  
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
                                  col_category=col_category,
                                  col_species=col_species,
                                  col_mechanism = col_mechanism,
                                  trans = trans)
      
      impact_species <- eicat_score_list %>%
        na.omit() %>%
        rownames()
    }
    
    if(type=="max"){
      eicat_score<-eicat_score_list[species_list,type]
      
    } else if(type=="mean"){
      eicat_score<-eicat_score_list[species_list,type]
      
      
    } else if(type=="max_mech"){
      eicat_score<-eicat_score_list[species_list,type]
    } else (stop("type should be one of max, mean or max_mech"))
    
    # site by species impact
    impactScore = sweep(sbs.taxon,2,eicat_score,FUN = "*")
    # Species impact sport
    speciesScore<-colSums(impactScore,na.rm = TRUE)/cube$num_cells
    
    speciesImpact<-bind_rows(speciesImpact,speciesScore)
  }
  
  
  speciesImpact<- speciesImpact%>%
    select(any_of(impact_species)) %>% 
    as.data.frame()
  
  rownames(speciesImpact)<-as.character(period)
  return(speciesImpact)
}



# speciesImpact %>%
#   rownames_to_column("year") %>%
#   mutate(year=as.numeric(year)) %>%
#   gather(-year,key = "Alien_species", value = "impact_score") %>%
#   ggplot(aes(x = year, y = impact_score)) +
#   geom_line(aes(color = Alien_species),linewidth=1.5)+
#   theme_minimal() +
#   labs(
#     title = "species impact",
#     y = "impact score"
#   )+
#   theme(text=element_text(size=14))
