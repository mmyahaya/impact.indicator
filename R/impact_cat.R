impact_cat<-function(impact_data,
                     species_list,
                     col_impact=NULL,
                     col_name=NULL){


  if(all(c("impact_category","scientific_name")%in%names(impact_data))){
    impact_data<-impact_data
  } else if((!is.null(col_impact)&!is.null(col_name))){
    impact_data <- impact_data %>%
      rename(all_of(c(impact_category=col_impact, scientific_name=col_name)))

  } else{ stop("impact_category and scientific_name are not found in  impact_data. col_impact and col_name must be given")}

  category_max_mean <- impact_data %>%
    mutate(impact_category = substr(impact_category, 1, 2)) %>%
    filter(impact_category %in% c("MC", "MN", "MO", "MR", "MV")) %>%
    select(scientific_name, impact_mechanism, impact_category) %>% 
    mutate(category_value=cat_num(impact_category,3)) %>%
    # mutate(category_value = case_when(
    #   impact_category == "MC" ~ 0,
    #   impact_category == "MN" ~ 1,
    #   impact_category == "MO" ~ 2,
    #   impact_category == "MR" ~ 3,
    #   impact_category == "MV" ~ 4,
    #   TRUE ~ 0  # Default case
    # )) %>%
    distinct(scientific_name, impact_mechanism, impact_category, 
             .keep_all = TRUE) %>%
    group_by(scientific_name) %>%
    summarise(
      max_value = max(category_value, na.rm = TRUE),
      mean_value = mean(category_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(scientific_name %in% species_list)


  category_sum_mech = impact_data %>%
    dplyr::mutate(impact_category=substr(impact_category,1,2)) %>%
    dplyr::filter(impact_category %in% c("MC","MN","MO","MR","MV")) %>%
    dplyr::select(scientific_name,impact_mechanism,impact_category) %>%
    mutate(category_value=cat_num(impact_category,3)) %>% 
    # dplyr::mutate(category_value = case_when(
    #   impact_category == "MC" ~ 0,
    #   impact_category == "MN" ~1,
    #   impact_category == "MO" ~2,
    #   impact_category == "MR" ~3,
    #   impact_category == "MV" ~4,
    #   TRUE ~ 0  # Default case, if any value falls outside the specified ranges
    # )) %>%
    distinct(scientific_name,impact_mechanism,impact_category, 
             .keep_all = TRUE) %>%

    dplyr::group_by(scientific_name,impact_mechanism) %>%
    dplyr::summarise(across(category_value,max),
                     .groups = "drop") %>%
    dplyr::group_by(scientific_name) %>%
    dplyr::summarise(across(category_value,sum),
                     .groups = "drop") %>%
    dplyr::filter(scientific_name %in% species_list)


   


    category_M<-dplyr::left_join(category_max_mean,category_sum_mech,
                        by=join_by(scientific_name)) %>%
    tibble::column_to_rownames(var = "scientific_name")

    names(category_M)<-c("max","mean","max_mech")


  #impact_matrix<-category_M %>% dplyr::select(all_of(fun))
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