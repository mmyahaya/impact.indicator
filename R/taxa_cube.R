#' @title Prepare Data Cubes
#'
#'@description Prepare data cube to calculate species impact . The function `taxa_cube`
#' can take in the scientific name of the taxa of interest as in character or
#' GBIF occurrences data containing necessary columns. The GBIF occurrences is
#' downloaded if scientific names is given. The function returns data cubes of
#' focal taxa and reference taxa. If the reference taxa is not specified, the
#' function returns the focal taxa as reference taxa.
#'@param taxa Character or dataframe. The character should be the scientific
#'name of the focal taxa while the dataframe is the GBIF occurrences data which must
#'contain "decimalLatitude","decimalLongitude","species","speciesKey",
#'"coordinateUncertaintyInMeters","dateIdentified", and "year".
#'@param country.sf sf object. The shapefile of the region of study
#'@param limit Integer. Number of records to return from GBIF download.
#'Default is set to 500

#' @param country Character. Country for which the GBIP occurrences data should
#' be downloaded. Country should be 2-letter country code (ISO-3166-1).
#' Default is ZA
#' @param res Numeric. The resolution of grid cells to be used. Default is 0.25
#'
#' @return A dataframe containing the `sim_cubes` of taxa.
#'

taxaFun <- function(taxa,country.sf,limit=500,
                    country='ZA',
                    res=0.25,
                    first_year=NULL){

  grid <- country.sf %>%
    sf::st_make_grid(cellsize = c(res,res),
                     offset = c(sf::st_bbox(country.sf)$xmin,
                                sf::st_bbox(country.sf)$ymin)) %>%
    sf::st_sf() %>%
    dplyr::mutate(cellid = dplyr::row_number())

  # download taxaif the scientific name is given as character
  if("character" %in% class(taxa)){
    taxa.gbif_download = rgbif::occ_data(scientificName=taxa, # download data from gbif
                                         country=country,
                                         hasCoordinate=TRUE,
                                         hasGeospatialIssue=FALSE,
                                         limit = limit)

    taxa.df = as.data.frame(taxa.gbif_download$data) #extract data from the downloaded file
  } else if("data.frame" %in% class(taxa)){ #check if data fame contains the required columns
    if(any(!c("decimalLatitude","decimalLongitude",
              "species","speciesKey","coordinateUncertaintyInMeters",
              "dateIdentified","year") %in% colnames(taxa))){
      requiredcol<-c("decimalLatitude","decimalLongitude","species",
                     "speciesKey","iucnRedListCategory","coordinateUncertaintyInMeters","dateIdentified","year")
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
                  species,speciesKey,iucnRedListCategory,
                  coordinateUncertaintyInMeters,year) %>% #select occurrence data
    dplyr::filter_all(all_vars(!is.na(.))) %>% # remove rows with missing data
    dplyr::filter(coordinateUncertaintyInMeters<=res*1000) %>%
    #dplyr::mutate(coordinateUncertaintyInMeters = coordinateUncertaintyInMeters/(res*1000)^2) %>%
    #dplyr::mutate(dateIdentified = as.Date(dateIdentified)) %>%  # convert date to date format
    sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
                 crs = 4326) %>%
    sf::st_join(grid) %>%
    as.data.frame() %>%
    dplyr::select(-geometry) %>%
    dplyr::mutate(occurrences=1)
  # %>%
  #   dplyr::group_by(species,coordinateUncertaintyInMeters,year,speciesKey,
  #                   iucnRedListCategory,
  #                   cellid) %>%
  #   dplyr::summarise(across(occurrences, sum), .groups = "drop")


  taxa_cube<-b3gbi::process_cube(taxa.sf,grid_type = "custom",
                                 cols_cellCode = "cellid",
                                 cols_year = "year",
                                 cols_occurrences = "occurrences",
                                 cols_species = "species",
                                 cols_speciesKey = "speciesKey",
      cols_minCoordinateUncertaintyInMeters = "coordinateUncertaintyInMeters",
      first_year = first_year)


  return(taxa_cube)
}


acacia_cube$data
# countries_sf<-readRDS(paste0(getwd(),"/Data/countries_shapefile.rds"))
# SA.sf<-filter(countries_sf,name=="South Africa") %>% select(name,geometry)
# taxa_cube<-taxaFun(taxa = occ_Fabacae, country.sf = SA.sf, res=1)
# full_species_list<-sort(unique(taxa_cube$data$scientificName))
SA.sf<-sf::st_read("C:/Users/26485613/OneDrive - Stellenbosch University/Downloads/Code_Data/Code_Data/boundary_south_africa_land_geo.shp")


acacia_cube<-taxaFun(taxa = taxa_Acacia, country.sf = SA.sf, res=0.25,first_year=2010)


calc_ts.obs_richness(acacia_cube)
occ_density_ts(taxa_cube)




