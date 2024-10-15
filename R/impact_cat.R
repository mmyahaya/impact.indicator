impact_cat<-function(data,
                     species_list,
                     col_impact=NULL,
                     col_name=NULL,
                     fun="max"){
  
  
  if(all(c("impact_category","scientific_name")%in%names(data))){
    data<-data
  } else if((!is.null(col_impact)&!is.null(col_name))){
    data <- data %>% 
      rename(all_of(c(impact_category=col_impact, scientific_name=col_name)))
    
  } else{ stop("required column is not given")}
  
  
  if(fun=="max"){
    f<-function(x) max(x,na.rm = T)
  } else if(fun=="min"){
    f<-function(x) min(x,na.rm = T)
  } else if(fun=="mean"){
    f<-function(x) mean(x,na.rm = T)
  } else {stop("`fun` should be max, min or mean character")}
  
  category_M = data %>% 
    dplyr::mutate(impact_category=substr(impact_category,1,2)) %>% 
    dplyr::filter(impact_category %in% c("MC","MN","MO","MR","MV")) %>% 
    dplyr::select(scientific_name,
           impact_category) %>% 
    dplyr::mutate(category_value = case_when(
      impact_category == "MC" ~ 0,
      impact_category == "MN" ~1,
      impact_category == "MO" ~2,
      impact_category == "MR" ~3,
      impact_category == "MV" ~4,
      TRUE ~ 0  # Default case, if any value falls outside the specified ranges
    )) %>% 
    dplyr::group_by(scientific_name,impact_category) %>%
    
    dplyr::summarise("impact"=across(category_value, first),
                     "frequency"= across(category_value, length),
                     .groups = "drop") %>% 
   
    tidyr::pivot_wider(names_from = impact_category, 
                       values_from = c(impact,frequency)) %>% 
   
    dplyr::filter(scientific_name %in% species_list) %>% 
   
    tibble::column_to_rownames(var = "scientific_name") 
  
  
  df_impact <- category_M %>% 
    select(starts_with("impact"))
  
  impact_matrix<-as.data.frame(apply(df_impact,1,f))
  
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

  
  # df_frquency <- category_M %>% 
  #   select(starts_with("frequency"))
  
  return(impact_matrix)
  
  
}


#impact_cat(Combined_eicat_data,species_list)
