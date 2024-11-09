#' Title
#'
#' @param impact_data The data frame of species impact with columns of category,
#'  species
#' @param species_list 
#' @param col_category 
#' @param col_species 
#' @param trans 
#' @param col_mechanism 
#'
#' @return
#' @export
#'
#' @examples
impact_cat<-function(impact_data,
                     species_list,
                     col_category=NULL,
                     col_species=NULL,
                     col_mechanism=NULL,
                     trans=1){


  if(all(c("category",
           "species",
           "mechanism")%in%names(impact_data))){
    impact_data<-impact_data
  } else if(all(c(!is.null(col_category),
                  !is.null(col_species),
                  !is.null(col_mechanism)))){
    impact_data <- impact_data %>%
      rename(all_of(c(category=col_category,
                      species=col_species,
                      mechanism=col_mechanism)))

  } else{ stop("category, species and mechanism are not found in  impact_data. col_category and col_species must be given")}

  category_max_mean <- impact_data %>%
    mutate(category = substr(category, 1, 2)) %>%
    filter(category %in% c("MC", "MN", "MO", "MR", "MV")) %>%
    select(species, mechanism, category) %>% 
    mutate(category_value=cat_num(category,trans)) %>%
    distinct(species, mechanism, category, 
             .keep_all = TRUE) %>%
    group_by(species) %>%
    summarise(
      max_value = max(category_value, na.rm = TRUE),
      mean_value = mean(category_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(species %in% species_list)


  category_sum_mech = impact_data %>%
    dplyr::mutate(category=substr(category,1,2)) %>%
    dplyr::filter(category %in% c("MC","MN","MO","MR","MV")) %>%
    dplyr::select(species,mechanism,category) %>%
    mutate(category_value=cat_num(category,trans)) %>% 
    distinct(species,mechanism,category, 
             .keep_all = TRUE) %>%

    dplyr::group_by(species,mechanism) %>%
    dplyr::summarise(across(category_value,max),
                     .groups = "drop") %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(across(category_value,sum),
                     .groups = "drop") %>%
    dplyr::filter(species %in% species_list)

    category_M<-dplyr::left_join(category_max_mean,category_sum_mech,
                        by=join_by(species)) %>%
    tibble::column_to_rownames(var = "species")

    names(category_M)<-c("max","mean","max_mech")

  impact_matrix<-category_M

  na.df<-as.data.frame(matrix(NA,
                              nrow = length(setdiff(species_list,
                                                    rownames(impact_matrix))),
                              ncol = ncol(impact_matrix)))
  row.names(na.df)<-setdiff(species_list,rownames(impact_matrix))
  names(na.df)<-names(impact_matrix) # column names
  impact_matrix<-rbind(impact_matrix,na.df)
  
  
  impact_matrix <- impact_matrix %>%
    dplyr::mutate(rowname = row.names(.)) %>%
    dplyr::arrange(rowname) %>%
    dplyr::select(-rowname)

  return(impact_matrix)
}




cat_num<-function(cat,trans){
  name<-c("MC", "MN", "MO", "MR", "MV")
  if(trans==1){
    x<-c(0,1,2,3,4)
    names(x)<-name
  } else if (trans==2){
    x<-c(1,2,3,4,5)
    names(x)<-name
  }
  else if (trans==3){
    x<-c(0,1,10,100,1000,10000)
    names(x)<-name
  } else { stop("`trans` must be 1, 2, or 3")}
  
  return(x[cat])
}
