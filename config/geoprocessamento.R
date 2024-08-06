library(stars)
library(here)
library(tidyverse)

#nome_shape_area_estudo <- here("/home/reginaldo/UTFPR/projetos/modelagem_de_especies/castela_tweedii/entrada/shape_area_estudo/lim_pais_a[6933].shp")  


#shape_area_estudo <- sf::read_sf(nome_shape_area_estudo, quiet = T)

repair_shp <- function(shp) {
  #if (sf::st_is(x, c("POLYGON", "MULTIPOLYGON")))
  #if (!(sf::st_geometry(shp) %>% inherits(c("sfc_POLYGON","sfc")))){
  if (!(sf::st_geometry(shp) %>% class() %in% "sfc" %>% any())){
    stop("Invalid study area file!")
  }
  if (shp %>% sf::st_crs() %>% is.na()){
    shp <- shp %>% 
      sf::st_set_crs(4326)
  } 
  shp <- shp %>% sf::st_make_valid() 
  return(shp)
}

drop_areas <- function(shp, inf_limit=0) {
  if (!(sf::st_geometry(shp) %>% class() %in% "sfc" %>% any())){
    stop("Invalid study area file!")
  }
  if (inf_limit<=0){
    stop("Invalid area limit!")
  }
  
  if (shp %>% sf::st_is(c("POLYGON", "MULTIPOLYGON")) %>% all()) {
    kept_areas <- shp[shp %>% sf::st_area() >= units::set_units(inf_limit,"m^2"),]

    if (shp %>% sf::st_area() %>% sum() != kept_areas %>% sf::st_area() %>% sum()) {
      return(kept_areas)
    }
  }
  return(shp)      
}

make_grid <- function(shp, cell_width=0, cell_height=0, var_names=NULL, centroid=T, epsg=NULL){
  sfc_as_cols <- function(x, geometry, names = c("x","y")) {
    if (missing(geometry)) {
      geometry <- sf::st_geometry(x)
    } else {
      geometry <- rlang::eval_tidy(enquo(geometry), x)
    }
    stopifnot(inherits(x,"sf") && inherits(geometry,"sfc_POINT"))
    ret <- sf::st_coordinates(geometry)
    ret <- tibble::as_tibble(ret)
    stopifnot(length(names) == ncol(ret))
    x <- x[ , !names(x) %in% names]
    ret <- setNames(ret,names)
    dplyr::bind_cols(x,ret)
  }
  
  if (!(sf::st_geometry(shp) %>% class() %in% "sfc" %>% any())){
    stop("Invalid study area file!")
  }
  if (cell_width<=0 || cell_height<=0){
    stop("Invalid cell size!")
  }
  shp <- shp %>%
    dplyr::rename_all(tolower)
  
  if (!is.null(var_names)){
    shp <- shp %>%
      dplyr::select(var_names %>% tolower() %>% all_of()) 
  }
  
  if(!is.null(epsg)){
    shp <- st_transform(shp, st_crs(paste0("+init=epsg:",epsg)))
  }
  
  bbox <- shp %>% 
    sf::st_bbox() %>% 
    as.vector()
  
  xmin <- bbox[1]
  ymin <- bbox[2]
  xmax <- bbox[3]
  ymax <- bbox[4]
  
  n_x_cell <- (abs(xmax - xmin) / cell_width) %>%
    ceiling()
  n_y_cell <- (abs(ymax - ymin) / cell_height) %>%
    ceiling()
  
  new_xmin <- ((xmin / cell_width) %>% trunc()) * cell_width
  new_ymin <- ((ymin / cell_height) %>% trunc()) * cell_height
  new_xmax <- new_xmin + (n_x_cell * cell_width)
  new_ymax <- new_ymin + (n_y_cell * cell_height)
  
  attr(sf::st_geometry(shp), "bbox") <- sf::st_bbox(
    c(xmin=new_xmin, xmax=new_xmax, ymax=new_ymax, ymin=new_ymin) %>% 
      sf::st_bbox()
  )
  
  shp_grid <- shp %>% 
    sf::st_make_grid(what = "polygons", n=c(n_x_cell, n_y_cell))
  
  shp_grid <- st_sf(geometry=shp_grid) %>%
    rowid_to_column("cell_id")
  
  
  grid_cells <- shp_grid %>% 
    sf::st_intersects(shp) %>% 
    as.data.frame() %>%
    rename(cell_id=row.id, area_id=col.id) %>%
    dplyr::left_join(shp_grid, by = "cell_id")
  
  shp_grid <- grid_cells %>% 
    dplyr::left_join(
      shp %>% 
        as.data.frame() %>% 
        tibble::rowid_to_column("area_id") %>% 
        dplyr::select(-geometry),
      by = "area_id"
    ) %>%
    dplyr::group_by(cell_id) %>% 
    dplyr::summarise_all(~ ifelse(is.numeric(.), mean(.), .), na.rm = TRUE) %>%
    dplyr::select(-area_id) %>%
    #tibble::rowid_to_column("cell_id") %>%
    st_sf(crs=sf::st_crs(shp), sf_column_name = "geometry")
  
  if (centroid){
    centroids <-  shp_grid %>%
      sf::st_geometry() %>%
      sf::st_centroid() %>% 
      sf::st_coordinates() %>%
      as.data.frame()
    
    shp_grid <- shp_grid %>%
      dplyr::bind_cols(centroids) %>%
      dplyr::rename(x_centroid = X, y_centroid = Y)
  }
  return(shp_grid)
}

area_map <- function(shp, title="", crs_subtitle=T, lat="decimalLatitude", long="decimalLongitude", group="species", colour="black", fill=NA) {
  number_format <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE)
  
  map_data <- NULL
  ifelse(is(shp, "Spatial"), map_data <- fortify(shp), map_data <- shp)
  
  map_tmp <- ggplot(
    data =  map_data,
    aes_string(
      x = long, 
      y = lat, 
      group = group
    )
  ) + 
    scale_x_continuous(labels = number_format) +
    scale_y_continuous(labels = number_format) +
    coord_equal()
  
  if (title!=""){
    map_tmp <- map_tmp + 
      labs(title=title)
  } 
  
  if (crs_subtitle && is(shp, "Spatial")){
    map_tmp <- map_tmp + 
      labs(subtitle = paste0(crs(shp)))
  }
  
  if (is.na(fill)){
    map_tmp <- map_tmp + 
      geom_polygon(colour = colour, fill = NA)
  } else {
    map_tmp <- map_tmp + 
      geom_polygon(colour = NA, aes_string(fill = fill))
  }
  
  return(map_tmp)  
}


#shp = grid_study_area
#raster_folder= folder_current_rasters
#var_names= raster_vars
#scenario = 'ssp370'



add_raster <- function(shp, raster_folder=NULL, var_names=NULL, scenario=NULL){
  if (is.null(raster_folder)){
    stop("Invalid raster folder!")
  }
  if (is.null(var_names)){
    var_names <- raster_folder %>%
      list.files()
  } 
  if(is.null(scenario)){ # if scenario == NULL, we are using current data
    l <- raster_folder %>% list.files(full.names = T)
    var_names2 <- paste0(var_names, '\\D')
    raster_list <- l[grep(paste(var_names2,collapse="|"),l)]
    
    if (length(raster_list)!=length(var_names)){
      stop("It is not possible to determine exactly which rasters will be used from the given list of names!")
    }
    
    raster_stack <- raster_list %>%
      stack() %>%
      tryCatch(error = function(e) e)
    
  } else { # if scenario != NULL, we are using future data
    l <- raster_folder %>% list.files(full.names = T, rec=T)
    l <- l[grep(scenario,l)]
    #if (!is.null(selected_gcms)){
    #  l <- l[grep(paste(selected_gcms,collapse="|"),l)]
    #}
    raster_stack <-  l %>%
      stack() %>%
      tryCatch(error = function(e) e)
  }
  
  # if tryCatch returns error, than we need to work on raster data:
  # 1. arrumar todos os rasters (extensão/resolução/...)
  if(any(class(raster_stack)=='error')){
    stop("Rasters are not compatible so they can't be stacked.")
  }
  
  if(!is.null(scenario)){
    names(raster_stack) <- paste0('bio_', 1:19)
    raster_stack <- raster_stack[[unlist(var_names)]]
  } else {
    raster_stack <- raster_stack[[str_sort(names(raster_stack), numeric=T)]]
    names(raster_stack) <- var_names
  }
  
  if(!as.character(crs(raster_stack))==as.character(CRS(crs(shp)))){
    raster_stack <- projectRaster(raster_stack, crs=CRS(crs(shp)))
  }
    
  # 2. Obter bbox do shape
  # 3. crop mask stack
  raster_stack <- raster_stack %>%
    crop(shp) %>%
    mask(shp) %>%
    stack()
  
  v <- lapply(shp$cell_id, function(i){
    sc <- subset(shp, cell_id==i)
    vals <- as.vector(terra::extract(raster_stack, sc, fun=mean, na.rm=T))
    return(vals)})
  
  #for(i in shp$cell_id){
  #  sc <- subset(shp, cell_id==i)
  #  vals <- as.vector(terra::extract(raster_stack, sc, fun=mean, na.rm=T))
  #  v[[i]] <- vals
  #}
  
  v <- v %>% 
    do.call(rbind,.) %>% 
    as.data.frame()
  
  v <- v %>% mutate(across(everything(), ~replace_na(.x, NA)))
  
  # 7. rasteriza shape (dentro do stars)
  shp_raster <- shp %>%
    cbind(v) %>%
    st_as_sf() 
  
  #n <- shp %>% length()
  
  names(shp_raster) <- c(names(shp)[-which(names(shp)=='geometry')], names(raster_stack), "geometry")
  
  # set classes:
  cols <- !lapply(shp_raster, class) == 'numeric' & !names(shp_raster) %in% c('geometry', 'cell_id')
    
  shp_raster[,cols] <- shp_raster[,cols] %>% 
    mutate(across(where(~ !is.numeric(.)), as.numeric))
  
  return(shp_raster)
}
#add_raster <- function(shp, raster_folder=NULL, cell_width=0, cell_height=0, var_names=NULL){
#  #if (!(shp %>% is("Spatial"))){
#  #  stop("Arquivo da área de estudo inválido!")
#  #}
#  if (is.null(raster_folder)){
#    stop("Pasta de rasters inválida!")
#  }
#  if (cell_width<=0 || cell_height<=0){
#    stop("Tamanho de célula inválido!")
#  }
#  
#  if (is.null(var_names)){
#    var_names <- raster_folder %>%
#      list.files()
#  } 
#  
#  raster_list <- raster_folder %>%
#    list.files(full.names = T) %>% 
#    keep(~ .x %>% str_detect(fixed(var_names, ignore_case = T)) %>% any())
#  
#  if (length(raster_list)!=length(var_names)){
#    stop("Não é possível determinar com exatidão quais rasters serão usados a partir da lista de nomes informada!")
#  }
#  
#  raster_stack <- raster_list %>%
#    stack()
#  
#  raster_tmp_file <- tempfile() %>% paste0(".tif")
#  raster_stack %>% 
#    writeRaster(
#      raster_tmp_file,
#      format = "GTiff",
#      bylayer = F, #bylayer = T,
#      #suffix = lista_rasters_bio %>% names(),
#      options = c("dstnodata =-9999.0"),
#      overwrite=T,
#    )
#  
#  shp_countour_file <- tempfile() %>% paste0(".shp")
#  shp %>%
#    raster::aggregate(dissolve=T) %>%
#    #rgeos::gBuffer(width=-(min(c(cell_width, cell_height)))/10, capStyle = "SQUARE", joinStyle = "BEVEL") %>%
#    #rgeos::gUnionCascaded() %>%
#    as("SpatialPolygonsDataFrame") %>% 
#    writeOGR(
#      dsn = shp_countour_file, 
#      layer=".", 
#      driver="ESRI Shapefile",
#      overwrite=T
#    )
#  
#  shp_area_file <- tempfile() %>% paste0(".shp")
#  shp %>%
#    writeOGR(
#      dsn= shp_area_file, 
#      layer=".", 
#      driver="ESRI Shapefile",
#      overwrite=T
#    )
#  
#  shp_grid_file <- tempfile() %>% paste0(".shp")
#  shp %>%
#    as("SpatialPolygonsDataFrame") %>%
#    writeOGR(
#      dsn= shp_grid_file, 
#      layer=".", 
#      driver="ESRI Shapefile",
#      overwrite=T
#    )
#  
#  raster_file_reescaled_countour <- tempfile() %>% paste0(".tif")
#  raster_reescaled_countour <- gdalwarp( ### sp::st_warp
#    raster_tmp_file,
#    raster_file_reescaled_countour,
#    s_srs = raster::crs(raster_stack),
#    t_srs = raster::crs(shp),
#    cutline = shp_countour_file,
#    crop_to_cutline = T, 
#    r = 'average', # resampling method
#    tr = c(cell_width, cell_height), # target resolution
#    tap = T, # should output extent include minimum extent?
#    te = shp %>% bbox() %>% as("vector"), #extent of output
#    te_srs = raster::crs(shp), # crs of output
#    dstnodata = "-9999.0",
#    ot = 'Float32', # data type
#    co = c("BIGTIFF=YES"), #"COMPRESS=DEFLATE", "PREDICTOR=2","ZLEVEL=9"),
#    #wo = c("CUTLINE_ALL_TOUCHED=TRUE"),
#    multi = T, # multithread
#    output_Raster = T, # sai um rasterBrick
#    overwrite = T,
#    verbose = T
#  ) %>% 
#    terra::crop(shp)
#  
#  # Alternative:
#  #          raster_reescaled_countour <- terra::extract(terra::rast(raster_tmp_file),
#  #                                                      terra::vect(shp), 
#  #                                                      fun='mean',
#  #                                                      cells=T, 
#  #                                                      na.rm=T,
#  #                                                      )
#  #          shp <- cbind(shp, raster_reescaled_countour)
#  #          if (!all(shp$cell_id == shp$ID)){
#  #            stop("Reescaling failed: cell IDs do not match.")
#  #          }
#  # End of alternative.
#  
#  raster_reescaled_countour_masked <- raster_reescaled_countour %>% 
#    terra::rast() 
#  
#  raster_reescaled_countour_masked <- raster_reescaled_countour_masked %>%
#    terra::mask(shp %>% terra::vect(), touches=F) %>% 
#    raster::stack()
#  
#  raster_grid <- gdal_rasterize(
#    shp_area_file,
#    tempfile() %>% paste0(".tif"),
#    #burn = 0,
#    a = "cell_id",
#    #at = T, 
#    co = c("BIGTIFF=YES"), 
#    a_nodata = "-9999.0",    
#    tr = c(cell_width, cell_height),
#    #tap= T,
#    ot = 'Float32',
#    output_Raster = T,
#    te = shp %>% bbox() %>% as("vector"),
#    verbose = T
#  )
#  
#  raster_grid <- raster_grid %>% 
#    terra::rast() %>%
#    terra::mask((raster_reescaled_countour_masked %>% terra::rast())[[1]])
#  
#  
#  grid_cells <- raster_grid %>%
#    as.vector() %>%
#    discard(is.na)
#  
#  shp@data <- shp@data %>%
#    as.data.frame()
#  
#  shp_grid <- shp[grid_cells %>% as.integer(),] 
#  
#  shp_grid@data$cell_id <- 1:length(grid_cells)
#  shp_grid %>% 
#    spChFIDs(as.character(1:length(grid_cells)))
#  
#  
#  shp_grid@data <- shp_grid@data %>% bind_cols(
#    raster_reescaled_countour_masked %>%
#      as.list() %>%
#      map_dfc(~ .x %>% values() %>% discard(is.na) %>% as.data.frame()) %>%
#      rename_all(~ (var_names %>% unlist()))
#  )
#  
#  return(shp_grid)
#}










