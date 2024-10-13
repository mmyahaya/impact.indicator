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
    mutate(impact_category=substr(impact_category,1,2)) %>% 
    filter(impact_category %in% c("MC","MN","MO","MR","MV")) %>% 
    select(scientific_name,
           impact_category) %>% 
    mutate(category_value = case_when(
      impact_category == "MC" ~ 0,
      impact_category == "MN" ~1,
      impact_category == "MO" ~2,
      impact_category == "MR" ~3,
      impact_category == "MV" ~4,
      TRUE ~ 0  # Default case, if any value falls outside the specified ranges
    )) %>% 
    group_by(scientific_name,impact_category) %>%
    #choose the first trait value if there are multiples trait for a species
    summarise(across(category_value, first), .groups = "drop") %>% 
    # reshape to wide format to have specie by trait dataframe
    pivot_wider(names_from = impact_category, values_from = category_value) %>% 
    # select species that are only present in gbif data
    filter(scientific_name %in% species_list) %>% 
    #convert species names to row names
    column_to_rownames(var = "scientific_name") 
  
  
  category_M<-data.frame("fun"=apply(category_M,1,f))
  #names(category_M)<-fun
  na.df<-as.data.frame(matrix(NA,
                              nrow = length(setdiff(species_list,rownames(category_M))),
                              ncol = ncol(category_M)))
  row.names(na.df)<-setdiff(species_list,rownames(category_M))
  names(na.df)<-names(category_M) # column names
  category_M<-rbind(category_M,na.df)
  category_M <- category_M %>%
    dplyr::mutate(rowname = row.names(.)) %>%  
    dplyr::arrange(rowname) %>%               
    dplyr::select(-rowname)                    
  return(category_M)
  
  
}