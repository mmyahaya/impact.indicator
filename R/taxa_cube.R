taxaFun <- function(taxa,country.sf,limit=500,country='ZA',res=0.25){
  
  grid <- country.sf %>%
    sf::st_make_grid(cellsize = c(res,res),
                     offset = c(sf::st_bbox(country.sf)$xmin,
                                sf::st_bbox(country.sf)$ymin)) %>%
    sf::st_sf() %>%
    dplyr::mutate(cellid = dplyr::row_number())
  
  # download taxa (target taxa) if the scientific name is given as character
  if("character" %in% class(taxa)){
    taxa.gbif_download = rgbif::occ_data(scientificName=taxa, # download data from gbif
                                         country=country,
                                         hasCoordinate=TRUE,
                                         hasGeospatialIssue=FALSE,
                                         limit = limit)
    
    taxa.df = as.data.frame(taxa.gbif_download$data) #extract data from the downloaded file
  } else if("data.frame" %in% class(taxa)){ #check if data fame contains the required columns
    if(any(!c("decimalLatitude","decimalLongitude",
              "species","speciesKey","coordinateUncertaintyInMeters","dateIdentified","year") %in% colnames(taxa))){
      requiredcol<-c("decimalLatitude","decimalLongitude","species","coordinateUncertaintyInMeters","dateIdentified","year")
      missingcol<-requiredcol[!c("decimalLatitude","decimalLongitude","species","speciesKey","coordinateUncertaintyInMeters","dateIdentified","year") %in% colnames(taxa)]
      cli::cli_abort(c("{missingcol} is/are not in the {.var taxa} column ", 
                       "x" = "{.var taxa} should be a data of GBIF format "))
    }
    # take taxa data frame if accurate
    taxa.df<-taxa
  } else { # stop and report if taxa is not a scientific name or dataframe
    cli::cli_abort(c("{.var taxa} is not a character or dataframe"))
  }
  
  
  taxa.sf = taxa.df %>%
    dplyr::select(decimalLatitude,decimalLongitude,
                  species,speciesKey,coordinateUncertaintyInMeters,dateIdentified,year) %>% #select occurrence data
    dplyr::filter_all(all_vars(!is.na(.))) %>% # remove rows with missing data
    dplyr::filter(coordinateUncertaintyInMeters<=res*1000) %>% 
    #dplyr::mutate(coordinateUncertaintyInMeters = coordinateUncertaintyInMeters/(res*1000)^2) %>% 
    dplyr::mutate(dateIdentified = as.Date(dateIdentified)) %>%  # convert date to date format
    sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
                 crs = 4326) %>% 
    sf::st_join(grid) %>% 
    as.data.frame() %>% 
    dplyr::select(-geometry) %>% 
    dplyr::mutate(occurrences=1) %>% 
    dplyr::group_by(species,coordinateUncertaintyInMeters,year,speciesKey,
                    iucnRedListCategory,
                    cellid) %>% 
    dplyr::summarise(across(occurrences, sum), .groups = "drop")
  
  
  taxa_cube<-b3gbi::process_cube(taxa.sf,grid_type = "custom",
                                 cols_cellCode = "cellid", cols_year = "year",
                                 cols_occurrences = "occurrences",
                                 cols_species = "species",cols_speciesKey = "speciesKey",
                                 cols_minCoordinateUncertaintyInMeters = "coordinateUncertaintyInMeters")
 
  
  return(taxa_cube)
}
