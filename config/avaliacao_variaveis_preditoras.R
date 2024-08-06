add_var_shp <- function(df, shp){
  new_df <- shp %>% 
    get_var_shp()
  
  if (is.null(df) || missing(df)){
    return(new_df)
  }
  df %>% 
    bind_cols(new_df) %>%
    return()
}

get_var_shp <- function(shp){
  if (is.character(shp)){
    shp <- shp %>% 
      readOGR() 
  }
  shp@data <- shp@data %>%
    select(-(c("cell_id") %>% any_of()))
  
    return(shp@data)
}

get_predictors_as_df <- function(shp){
  df_var_preditoras <- shp %>%
    readRDS() %>%
    st_as_sf() %>% 
    as.data.frame() %>%
    select(-c('geometry','x_centroid','y_centroid'))
  return(df_var_preditoras)
}

remove_var <- function(df, var_names){
  if (!missing(var_names) && (is.list(var_names) || is.vector(var_names))){
    df %>%
      select(-(var_names %>% 
          to_snake_case() %>% 
          abbreviate(minlength = 10) %>%
          any_of()
          )
      )
  } else {
    df
  }
}

center_scale <- function(df, var_names){
  if (!missing(var_names) && (is.list(var_names) || is.vector(var_names))){
    if(!all(var_names %in% colnames(df))){
      print('Not every variable chosen is in data.frame. Keeping only with available variables.')
      var_names <- var_names[which(var_names %in% colnames(df))]
      df %>%
        mutate_at(var_names, ~ as.numeric(scale(.))) %>%
        as.data.frame()
    } else {
      df %>%
        mutate_at(var_names, ~ as.numeric(scale(.))) %>%
        as.data.frame()
    }
    
  } else {
    df
  }
}

corr_plot <- function(df){
  my_ggcorrplot(
    df %>% na.omit() %>% cor() %>% round(2),
    method = "circle",
    type = "upper",
    tl.col = "black",
    tl.srt = 45,
    lab = T,
    lab_size = 3,
    p.mat = df %>% cor_pmat()
  )
}



#data=df_var_pca_bio
#nrEixos=4
#pcaName="bio"

calc_pca <- function(data, nrEixos=5, pcaName=""){
  a_pca <- data %>% PCA(scale.unit = T, graph = F, ncp = nrEixos)
  
  # Loadings (i.e. standard coordinates) are not given by FactoMineR's methods. They return principal coordinates.
  # You can calculate them by dividing variables' coordinates on a dimension by this dimension's eigenvalue's square root.
  a_pca$var$loadings <- a_pca$var$coord %>% 
    sweep(2, a_pca$eig[1:ncol(a_pca$var$coord),1] %>% sqrt(), FUN="/") %>%
    as.data.frame()
  if (!is.na(pcaName)){
    colnames(a_pca$var$coord) <- paste(colnames(a_pca$var$coord), pcaName) %>% make_clean_names()
    colnames(a_pca$var$cos2) <- paste(colnames(a_pca$var$cos2), pcaName) %>% make_clean_names()
    colnames(a_pca$var$cor) <- paste(colnames(a_pca$var$cor), pcaName) %>% make_clean_names()
    colnames(a_pca$var$contrib) <- paste(colnames(a_pca$var$contrib), pcaName) %>% make_clean_names()
    colnames(a_pca$var$loadings) <- paste(colnames(a_pca$var$loadings), pcaName) %>% make_clean_names()
    
    colnames(a_pca$ind$coord) <- paste(colnames(a_pca$ind$coord), pcaName) %>% make_clean_names()
    colnames(a_pca$ind$cos2) <- paste(colnames(a_pca$ind$cos2), pcaName) %>% make_clean_names()
    colnames(a_pca$ind$contrib) <- paste(colnames(a_pca$ind$contrib), pcaName) %>% make_clean_names()
  }
  
  return(a_pca)
}

proj_pca <- function(a_pca, list_scenarios, vars=NULL){
  if(is.null(vars)){
    vars <- rownames(a_pca$var$coord)
  }
  s <- lapply(list_scenarios, function(x){
    p <- x %>% 
      as.data.frame() %>% 
      select(all_of(vars)) %>% 
      center_scale(vars) %>% 
      predict(a_pca,.)
    p <- cbind(x, p$coord)
    ids <- grep(pattern='Dim',colnames(p))
    colnames(p)[ids] <- colnames(a_pca$var$coord)
    return(p)
  })
  return(s)
}

dt_pca_summ <- function(a_pca){
  a_pca %>% 
    dimdesc(axes = 1:.$call$ncp) %>% 
    compact(~ .$quanti) %>% 
    map_df(~ as.data.frame(.$quanti) %>% rownames_to_column("variables"), .id="dim") %>%
    select(dim, variables, correlation) %>% 
    pivot_wider(names_from = dim, values_from = correlation) %>%
    mutate_if(is.numeric, ~ round(.,4)) %>%
    datatable(options = list(pageLength = 10, scrollX=T))
} 

pca_screeplot <- function(aPca){
  retainComp <- comp_pca_retain(aPca)
  aPlot <- fviz_screeplot(aPca, addlabels = TRUE) + 
    geom_vline(xintercept=unlist(retainComp), linetype="dashed", color = c("red", "blue", "green"), show.legend = T) +
    annotate("text", label="Broken Stick", x=retainComp$broken.stick, y=20, angle=45) + 
    annotate("text", label="Kaiser Mean", x=retainComp$kaiser.mean, y=25, angle=45) +
    annotate("text", label="Horn P", x=retainComp$horn.p, y=30, angle=45) 
  return(aPlot)
}


contrib_scree <- function(a_pca){
  aContribPlot <- fviz_pca_var(a_pca, col.var="contrib", 
                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                               legend.title = "Contrib",
                               repel = TRUE # 
  ) 
  
  aScreePlot <- pca_screeplot(a_pca) 
  
  aScreePlot + aContribPlot
}

contrib_dims <- function(a_pca){
  plotList <- list()
  nrEixos <- ncol(a_pca$var$coord)
  for (i in 1:nrEixos){
    plotList[[i]] <- a_pca %>% fviz_contrib(
      title = colnames(a_pca$var$contrib)[i],
      choice = "var",
      axes = i
    )
  }
  
  plotList[[i+1]] <- a_pca %>% fviz_contrib(
    title= paste("All", nrEixos, "dims."), 
    choice = "var", 
    axes=1:nrEixos
  )
  wrap_plots(plotList)
}


contrib_corr <- function(a_pca){
  my_ggcorrplot(a_pca$var$contrib, 
                method = "circle", 
                tl.col="black", 
                tl.srt=45, 
                lab=T,
                legend.title = "Contrib",
                lab_size = 2) 
}

cos2_corr <- function(a_pca){
  my_ggcorrplot(
    a_pca$var$cos2,
    method = "circle",
    tl.col = "black",
    tl.srt = 45,
    lab = T,
    legend.title = "Cos2",
    lab_size = 3
  )
}

pca_cos2 <- function(a_pca){
  fviz_pca_var(a_pca, gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE, col.var="cos2")
}

cos2_dims <- function(a_pca){
  plotList <- list()
  nrEixos <- ncol(a_pca$var$cos2)
  for (i in 1:nrEixos){
    plotList[[i]] <- a_pca %>% fviz_cos2(
      title = colnames(a_pca$var$cos2)[i], 
      choice = "var", axes=i
    ) + 
      ylab("Cos2")
  }
  
  plotList[[i+1]] <- a_pca %>% fviz_cos2(
    title= paste("All", nrEixos, "dims."), 
    choice = "var", 
    axes=1:nrEixos
  ) + 
    ylab("Cos2")
  
  wrap_plots(plotList)
}

pca_bi_plot <- function(a_pca, pres_aus){
  a_pca %>% fviz_pca_biplot(
    geom.ind = "point",
    repel = T,
    pointshape = 21,
    pointsize = 2.5,
    col.ind = as.factor(ifelse(pres_aus == 1, "Presence", "Absence")),
    fill.ind = as.factor(ifelse(pres_aus == 1, "Presence", "Absence")),
    col.var = "black",
    palette = c("#00AFBB", "#E7B800", "#FC4E07"),
    addEllipses = TRUE,
    legend.title = "Groups"
  )  # +
  #ggpubr::fill_palette("jco")  #+ # Indiviual fill color
  #ggpubr::color_palette("npg")      # Variable colors
}

comp_pca_retain <- function(factoMinerObject){
  eigs <- factoMinerObject$eig[,1]
  
  broken.stick.distribution <- lapply(
    X=1:length(eigs),
    FUN=function(x,n){return(tail((cumsum(1/x:n))/n,n=1))},
    n=length(eigs)
  ) %>% 
    unlist() * 100
  n.broken.stick <- (eigs/sum(eigs)*100 > broken.stick.distribution) %>% 
    discard(. != T) %>% 
    length()
  n.kaiser.mean <- (eigs > mean(eigs)) %>% 
    discard(. != T) %>% 
    length()
  invisible(capture.output(
    n.horn <- paran(factoMinerObject$call$X, quietly = T, status = F)$Retained
  ))
  return(list(broken.stick=n.broken.stick, kaiser.mean=n.kaiser.mean, horn.p = n.horn))
}


get_df_vif <- function(var, pa, var_names=NULL, sp_names=NULL){
  if(is.null(var_names)){
    var_names <- colnames(var)
    var_names <- var_names[!var_names %in% c('cell_id',"x_centroid", "y_centroid", "geometry" )]
    cat('var_names not informed. Variables detected: ',paste(var_names, collapse=', '))
  }
  
  if(is.null(sp_names)){
    sp_names <- colnames(pa)[-ncol(pa)]
    cat('sp_names not informed. Detected species: ',paste(sp_names, collapse=', '))
  }
  
  df_vif <- cbind(as.data.frame(var) %>% select(all_of(var_names)), 
                  as.data.frame(pa) %>% select(!'geometry'))
  l <- lapply(sp_names, function(x){df <- filter(df_vif,!!sym(x)== 1)
                                    df <- df %>% select(!all_of(sp_names))})
  names(l) <- sp_names
  return(l)
}











