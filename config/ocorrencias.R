

#o_file = spp_data
#sp_names = spp_names
#shp = grid_study_area
#proj = CRS(crs_occ_data)
#

occurrences_to_shapefile <- function(o_file, sp_names, shp){
  #if(!is.null(sp_names)){
  #  sp_names2 <- sp_names %>% 
  #    to_snake_case() %>% 
  #    abbreviate(minlength = 10)
  #}
  
  suppressMessages(
    sp_occurrences <- o_file %>%
      vroom() %>%
      clean_names() %>%
      remove_empty(c("rows", "cols"))
  )
  
  if(is.null(sp_names)){
    sp_names <- sp_occurrences$species %>% unique()
  }
  
  if (!is.na(sp_names) && !is.null(sp_names) && length(sp_names)>0){
    sp_occurrences <-  sp_occurrences %>%
      filter(species %in% sp_names)
  }

  
  if (is.character(shp)){
    shp <- shp %>% 
      st_read(verbose = F)
  }

  sp_occurrences <- sp_occurrences %>%
    mutate(species=species %>% to_snake_case() %>% abbreviate(minlength = 10))
  
  sp_occurrences$decimal_longitude <- sp_occurrences$decimal_longitude %>% 
    as.character() %>% 
    as.numeric()
  
  sp_occurrences$decimal_latitude <- sp_occurrences$decimal_latitude %>% 
    as.character() %>% 
    as.numeric()
  
  coordinates(sp_occurrences) <- ~decimal_longitude+decimal_latitude
  
  crs(sp_occurrences) <-  CRS("+init=epsg:4326")
  
  sp_occurrences %>%  
    st_as_sf() %>%
    st_transform(st_crs(shp)) %>%
    return()
}



#shp_occ=shape_ocorrencias
#shp_area=shape_grid_estudo
#sp_names=nomes_especies

occurrences_to_pa_shapefile <- function(shp_occ, shp_area, sp_names){
  if(nchar(sp_names)>10){sp_names <- sp_names %>% 
    to_snake_case() %>% 
    abbreviate(minlength = 10)}
  
  shp_area <- as(shp_area, "Spatial")
  shp_occ <- as(shp_occ, "Spatial")
  
  pa_matrix <- shp_area@data[,F]
  
  for (spp in sp_names){
    pa_matrix <- pa_matrix %>% 
      bind_cols(
        shp_area %>%
          over(shp_occ[shp_occ@data$species==spp, ])
      )
  }
  names(pa_matrix) <- sp_names
  
  pa_matrix[!is.na(pa_matrix)] <- 1
  pa_matrix[is.na(pa_matrix)] <- 0
  
  pa_matrix <- pa_matrix %>% 
    mutate_all(~ as.integer(.))
  
  rownames(pa_matrix) <- rownames(shp_area@data)
  
  #grid_pa_matrix <- shp_area 
  #grid_pa_matrix@data <-  pa_matrix
  
  grid_pa_matrix <- cbind(shp_area, pa_matrix)
  
  grid_pa_matrix %>%
    st_as_sf()%>%
    return()
}

map_of_occurrences <- function(shp_occ, shp_area, title="", crs_subtitle=T){
  map_tmp <- ggplot(
    data =  fortify(shp_area),
    aes(
      x = long, 
      y = lat,
      group=group
    )
  )

  map_tmp <- map_tmp +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = comma) +
    coord_equal()
  
  if (title!=""){
    map_tmp <- map_tmp + 
      labs(title=title)
  } 
  
  if (crs_subtitle){
    map_tmp <- map_tmp + 
      labs(subtitle = paste0(crs(shp_area)))
  }

  map_tmp <- map_tmp + 
    geom_polygon(colour = "black", fill=NA)
  
  map_tmp <- map_tmp + 
    geom_point(data = fortify(shp_occ), aes(x = long, y = lat, group=species, color=species)) 

  return(map_tmp)  
}

map_of_pa <- function(shp_pa, shp_area, sp_names){
  sp_names <- sp_names %>% 
    to_snake_case() %>% 
    abbreviate(minlength = 10)
  
  df_temp <- shp_pa@data %>% 
    rownames_to_column("id") %>%
    pivot_longer(-id, names_to = "especies", values_to="presenca") 
  
  legend_title  <- ifelse(length(sp_names)==1, sp_names[[1]], "Presença/Ausência")
  
  ggplot(data = fortify(shp_area) %>% left_join(df_temp)) +
       aes(x = long, y = lat, group = group) +
       geom_polygon(aes(fill = as.factor(presenca))) +
       scale_fill_manual(legend_title, values=c("darkseagreen2", "tomato")) +
       scale_x_continuous(labels = comma) +
       scale_y_continuous(labels = comma) +
       coord_equal()
}