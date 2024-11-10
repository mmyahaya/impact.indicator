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
impact_indicator<-function(cube,
                           impact_data = NULL,
                           col_category=NULL,
                           col_species=NULL,
                           col_mechanism=NULL,
                           trans=1,
                           type=NULL,
                           coords=NULL){


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
    # Site impact sport
    speciesScore<-colSums(impactScore,na.rm = TRUE)/cube$num_cells

    speciesScore%>%
      as.data.frame() %>%
      t() %>%
      as.data.frame()

    species_values<-bind_rows(species_values,speciesScore)

  }

  species_values<- species_values%>%
    select(any_of(impact_species))

  impact_values<-as.data.frame(impact_values)
  names(impact_values)<-c("year","value")
  rownames(species_values)<-as.character(period)
  return(list("impact_values"=impact_values,"species_values"=species_values,
              "sitedf"=coords))
}

