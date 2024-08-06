# Aqui são carregadas as bibliotecas e configurações necessárias para a execução dos documentos
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE, fig.width = 12, out.width="100%")

# Geoprocessamento
library(raster)                              # manipulação de rasters
library(rgdal)                               # readOGR, ler arquivos raster e shapefile
library(rasterDT)                            # rasterização
library(sf)                                  # manipulação de formas geométricas com interface com gdal, rgeos e proj4g
library(rgeos)                               # habilita operações geométricas
library(gdalUtils)

# Gráficos
library(ggplot2)                             # plotagem de gráficos
library(patchwork)                           # combinar gráficos ggplot independentes
library(cowplot)                             # combinar gráficos ggplot independentes
library(plotly)                              # adiciona interatividade nos gráficos
library(mapview)                             # ferramenta de visualização de mapas geográficos
library(ggfortify)                           # permite que o ggplot plote gráficos rasters e shapefiles
library(scales)                              # habilita o uso de diferentes projeções gegráficas nos gráficos de mapas
library(factoextra)                          # Visualização de PCA
library(paran)                               # Análise paralela de Horn
library(Rtsne)                               # Redutor de dimensionalidade T-SNE  
library(ggcorrplot)

# Tabelas
library(DT)                                  # renderiza tabelas interativas

# Processamento e computação paralela 
library(parallel)                            # habilita processamento paralelo
library(snow)                                # habilita computação paralela

# SDM
library(sdm)                                 # modelagem de distribuição de espécies
library(FactoMineR)                          # PCA
library(rdist)                               # distancia euclideana
library(usdm)                                # vif
#library(clusternor)                          # xmeans clustering
library(ade4)                                # Enfa

# Miscelânea
library(tidyverse)                           # funções de manipulação de dataframes e  listas
library(purrrlyr)                            # funções que auxiliam no uso conjunto do dplyr com o purrr
library(here)                                # encontra o caminho (pasta) do projeto de maneira inteligente  e segura
library(fs)                                  # funções de manipulação de pastas e arquivos
library(stringr)                             # manipulação de strings
library(data.table)                          # le arquivos csv
library(janitor)                             # faz limpeza de datasets e nomes de variáveis
library(snakecase)                           # transformação de strings em snake case
library(lubridate)                           # manipulação de datas
library(vroom)

select <- dplyr::select

options(java.parameters = "-Xmx1g", java.awt.headless="false")
source(here("config/ggcorrplot.R"))
source(here("config/xmeans.R"))
source(here("config/algoritmos_predicao.R"))
source(here("config/geoprocessamento.R"))
source(here("config/ocorrencias.R"))
source(here("config/avaliacao_variaveis_preditoras.R"))
source(here("config/treinamento_avaliacao.R"))

pastaRaiz <- here()

if (interactive()){
  #pastaRaiz <- params$title %>% make_clean_names() %>% here()
} else {
  knitr::opts_chunk$set(echo=TRUE, message=TRUE, warning=TRUE, error=TRUE, include=TRUE)
  myKnit = function(inputFile, encoding, out=NULL){
    #   if (is.null(out)){
    #     out <- inputFile %>% fs::path_file() %>% fs::path_ext_remove()
    #   } else{
    #     out <- out %>% fs::path(inputFile %>% fs::path_file() %>% fs::path_ext_remove())
    #   }
    # pastaRaiz <<- out
    # rmarkdown::render(inputFile, encoding = encoding, output_dir = out)
    rmarkdown::render(inputFile, encoding = encoding, envir = new.env())
  }
  pastaRaiz <- ""
}


#---------------------------------------------------------------------------------------------
# Executa a predição a partir dos modelos treinados e transforma dados de suitability em dados de frequencia e, posteriormente, em presença/ausência

# Essa função recebe como entrada: os modelos treinados; um dataframe em que linhas são as células da grid e as colunas as varáveis 
# preditoras; e uma lista com os nomes dos algoritmos de predição, disponíveis dentro dos modelos treinados,  para fazer a predição.
# A predição é feita para cada célula da grid com cada um dos diversos modelos dos diversos algoritmos de predição desejados. 
# Um limiar é usado para determinar se, em uma dada célula, existe presença ou ausência. O resultado é que, haverá uma matriz
# M x N, onde M são as células e N são as predições dos diversos modelos dos diversos algoritmos. Cada célula terá valor 0 ou 1,
# indicando presença ou ausência na célula. Como podem (e muito provavelmente) haverão diversos modelos de um mesmo algoritmo 
# (em função do uso de replicações como cross-validation ou bootstraping) que indicará ausencia para algumas execuções e presença
# para outras execuções, é necessário fazer um consenso para cada um dos algoritmos em ca célula. O consenso é a média dos 
# valores da célula para cada um dos algoritmos. Assim, o resultado é uma matriz M x N, onde M são as células da grid e N são os 
# consensos (médias para aquele algoritmo). O valor de cada célula será um valor entre 0 e 1. 
DRE_predict <- function(m_treinados, df_p, algoritmo_predicao, tipo_thresh=2, lista_thresh=NULL){
  if (is.null(lista_thresh)){
    lista_thresh <- getEvaluation(m_treinados, stat=c('threshold'), wtest="dep.test", opt=tipo_thresh)
  }
  predi <- predict(m_treinados, df_p, method=algoritmo_predicao) %>% 
                as.data.frame() %>%
                replace(., is.na(.), 0)
  
  predi[is.na(predi[])] <- 0
  
  colnames(predi) <- colnames(predi) %>% 
                        modify(., ~ str_remove(unlist(str_split(., "[.-]"))[1], "id_"))
  
  for (nome_coluna in colnames(predi)){
    thresh <- lista_thresh %>% 
                    filter(modelID==as.integer(nome_coluna))
    
    thresh <- thresh$threshold
    if (length(thresh)>0){
      predi[predi[nome_coluna] >= thresh, nome_coluna] <- 1
      predi[predi[nome_coluna] < thresh, nome_coluna] <- 0
    }
  }
  
  detalhes_modelos <- getModelInfo(m_treinados) %>% 
                          filter(method %in% algoritmo_predicao) %>%
                          select(modelID, method)
  
  for (i in 1:ncol(predi)){
    colnames(predi)[i] <- as.character(detalhes_modelos$method[detalhes_modelos$modelID==as.integer(colnames(predi[i]))])
  }
  
  colunas_por_metodo <- ncol(predi) / length(unique(colnames(predi)))
  predi_temp = data.frame(row.names = 1:nrow(predi))
  for (i in 1:length(unique(colnames(predi)))){
    nome_coluna <- colnames(predi)[(i-1)*colunas_por_metodo+1]
    predi_temp[[nome_coluna]] <- predi[,((i-1)*colunas_por_metodo+1):(i*colunas_por_metodo)] %>% rowSums(.)
    predi_temp[[nome_coluna]] <- predi_temp[[nome_coluna]] / colunas_por_metodo
  }
  
  return(predi_temp)
}
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
gerarPseudoAusenciasDRE <- function (grid, df_var, metodo="random"){
  df_var <- df_var %>% 
    select(-id_celula)
  
  df_tmp_p <- df_var %>% 
    filter(species == 1) 
  
  if (metodo == "k-means") {
    df_tmp <- df_var
    numero_grupos <- ceiling(nrow(df_tmp_p) / 10)
    
    set.seed(numero_grupos)
    grupos <- df_tmp %>% kmeans(numero_grupos)
    

    df_tmp$grupo <- as.factor(grupos$cluster)
    
    grupos <- df_tmp %>%
      group_by(grupo) %>%
      summarise(presencas = sum(species), quant_celulas= n(), dens_presenca=sum(species)/n(), presenca_rel=sum(species)/nrow(df_tmp_p))
    
    grupos <- grupos %>%
      filter(dens_presenca <= median(grupos$dens_presenca))
    
    df_bg <- df_tmp %>%
      filter(species==0) %>%
      filter(grupo %in% grupos$grupo) %>% 
      group_by(grupo) %>% 
      sample_n(nrow(df_tmp_p)/nrow(grupos)) %>% 
      as.data.frame() %>%
      select(-grupo, -species)
    
    # # visualização dos pontos de background
    # grid_temp <- grid
    # grid_temp@data <- df_tmp
    # grid_temp@data$id <- rownames(grid_temp@data)
    # grid_temp <- fortify(grid_temp, region = "id") %>% left_join(grid_temp@data)
    # 
    # ggplot(data = fortify(grid_temp)) +
    #   aes(x = long,
    #       y = lat,
    #       group = group,
    #       text = grupo) +
    #   geom_polygon(aes(fill = as.factor(grupo))) +
    #   #geom_point(
    #   #  data = df_bg %>% filter(species == 1),
    #   #  aes(x = long, y = lat, colour = species),
    #   #  size = 0.5
    #   #) +
    #   theme(
    #     axis.text.x = element_blank(),
    #     axis.text.y = element_blank(),
    #     axis.ticks.x = element_blank(),
    #     axis.ticks.y = element_blank()
    #   ) +
    #   coord_equal()
    
  } else if (metodo=="dre_area" || metodo=="dre_dens"){
    df_tmp <- df_var %>%
      select(-species) 
    
    quant_max <- 0
    num_grupo_quant_max <- 2
    dens_acum <- 0
    num_grupo_dens_acum <- 2
    
    numero_grupos <- trunc(nrow(df_tmp_p))
    set.seed(numero_grupos)
    for (i in 2:numero_grupos){
      suppressWarnings(grupos <- df_tmp %>% kmeans(i))
      
      df_tmp <- df_var
      df_tmp$grupo <- as.factor(grupos$cluster)
      
      grupos <- df_tmp %>%
        group_by(grupo) %>%
        summarise(presencas = sum(species), quant_celulas= n(), dens_presenca=sum(species)/n(), presenca_rel=sum(species)/nrow(df_tmp_p))
      
      grupos_area <- grupos %>%
        filter(dens_presenca <= median(grupos$dens_presenca))
      
      grupos_dens <- grupos %>%
        filter(dens_presenca > median(grupos$dens_presenca))
      
      if (sum(grupos_area$quant_celulas) > quant_max){
        quant_max <- sum(grupos_area$quant_celulas)
        num_grupo_quant_max <- i
      }
      
      if (sum(grupos_dens$presencas) / sum(grupos_dens$quant_celulas) > dens_acum){
        dens_acum <- sum(grupos_dens$presencas) / sum(grupos_dens$quant_celulas)
        num_grupo_dens_acum <- i
      }
    }
    
    if (metodo=="dre_area"){
      numero_grupos <- num_grupo_quant_max  
    } else {
      numero_grupos <- num_grupo_dens_acum
    }
    
    suppressWarnings(grupos <- df_tmp %>% kmeans(numero_grupos))
    
    df_tmp <- df_var
    df_tmp$grupo <- as.factor(grupos$cluster)
    
    grupos <- df_tmp %>%
      group_by(grupo) %>%
      summarise(
        presencas = sum(species), 
        quant_celulas= n(), 
        dens_presenca=sum(species)/n(), 
        presenca_rel=sum(species)/nrow(df_tmp_p)
      )
    
    if (metodo=="dre_area"){
      grupos <- grupos %>%
        filter(dens_presenca <= median(grupos$dens_presenca))
    }
    else {
      grupos <- grupos %>%
        filter(dens_presenca > median(grupos$dens_presenca))
    }
    
    
    df_bg <- df_tmp %>% 
      filter(species == 0) %>%
      filter(grupo %in% grupos$grupo) %>% 
      group_by(grupo) %>% 
      sample_n(ceiling(n() * nrow(df_tmp_p) / grupos$quant_celulas %>% sum())) %>% 
      as.data.frame() %>%
      select(-grupo, -species)
    
  } else {
    df_bg <- df_var %>%
      filter(species==0) %>%
      select(-species) %>%
      unique() %>%
      sample_n(nrow(df_tmp_p))
  }
  return (df_bg)
}



min_max <- function(x){
  (x - min(x)) / (max(x) - min(x))
}

mapaDistrConsenso <- function(df_predicao_distr_pa, shapefile_grid_estudo){
  shp_estudo <- shapefile_grid_estudo
  shp_estudo@data <- df_predicao_distr_pa
  df_temp <- shp_estudo %>%
    fortify() %>%
    merge(shp_estudo %>% as.data.frame(), by.x="id", by.y=0) %>%
    mutate(id_celula=id)
  
  if (nrow(df_temp %>% filter(!(consenso==0 | consenso==1)))==0){
    mapa_temp <- ggplot(data = df_temp) + 
      aes(x = long, y = lat, group = group, text=consenso) +
      geom_polygon(aes(fill = as.factor(consenso))) +
      scale_fill_manual(values=c("darkseagreen2", "tomato")) +
      theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) +
      coord_equal() 
  } else {
    mapa_temp <- ggplot(data = df_temp) + 
      aes(x = long, y = lat, group = group, fill=consenso, text=consenso) +
      geom_polygon() +
      scale_fill_continuous(low="darkseagreen2", high="tomato", limits=c(0,1), name = "Adequabilidade") +
      theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) +
      coord_equal()
  }
  return (mapa_temp)
}

mapaDistrMetodos <- function(df_predicao_distr_pa, shapefile_grid_estudo, metodos){
  shp_estudo <- shapefile_grid_estudo
  shp_estudo <- df_predicao_distr_pa
  nome_metodos <- metodos
  df_temp <- shp_estudo %>%
    fortify() %>%
    merge(shp_estudo %>% as.data.frame(), by.x="id_celula", by.y=0) %>%
    mutate(id_celula=id)
  
  df_temp <-  df_temp %>%
    gather(modelo, consenso, all_of(nome_metodos))
  
  if (nrow(df_temp %>% filter(!(consenso==0 | consenso==1)))==0){
    mapa_temp <- ggplot(data = df_temp) + 
      aes(x = long, y = lat, group = group) +
      geom_polygon(aes(fill = as.factor(consenso))) +
      scale_fill_manual(values=c("darkseagreen2", "tomato")) +
      theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) +
      coord_equal() +
      facet_wrap(~modelo)
  } else {
    mapa_temp <- ggplot(data = df_temp) + 
      aes(x = long, y = lat, group = group, fill=consenso) +
      geom_polygon() +
      scale_fill_continuous(low="darkseagreen2", high="tomato", limits=c(0,1), name = "Adequabilidade") +
      theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) +
      coord_equal() +
      facet_wrap(~modelo) 
  }

  return (mapa_temp)
}

mapaDistrModelos <- function(df_distr, grid_estudo, nome_modelos){
  return (mapaDistrMetodos(df_distr, grid_estudo, nome_modelos))
}

matriz_confusao <- function(obs, pre) {
  cmx<-matrix(nrow=2,ncol=2)
  colnames(cmx) <- rownames(cmx) <- c('P','A')
  cmx[1,1] <- length(which(obs == 1 & pre == 1))
  cmx[2,2] <- length(which(obs == 0 & pre == 0))
  cmx[1,2] <- length(which(obs == 0 & pre == 1))
  cmx[2,1] <- length(which(obs == 1 & pre == 0))
  cmx[] <- as.numeric(cmx)
  return (cmx)
}


calcular_metricas_matriz <- function(mc) {
  TP<-mc[1,1]
  FP<-mc[1,2]
  TN<-mc[2,2]
  FN<-mc[2,1]
  valor_mcc <- ((TP*TN)-(FP-FN)) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  valor_recall <- (TP / (TP+FN))
  
  return (round(c(tp=TP, fp=FP, tn=TN, fn=FN, mcc=valor_mcc, recall=valor_recall), 3))
}


calcular_limiares_modelo <- function(obs, pre){
  th <- sort(unique(pre))
  aval <- matrix(nrow=length(th),ncol=7)
  colnames(aval) <- c("pre", "tp","fp","tn","fn","mcc", "recall")
  aval[,1] <- sort(unique(pre))
  
  for (i in seq_along(th)) {
    w <- which(pre >= th[i])
    pt <- rep(0,length(pre))
    pt[w] <- 1
    aval[i,2:7] <- calcular_metricas_matriz(matriz_confusao(obs,pt))
  }
  return (as.data.frame(aval))
}

calcular_limiar_metrica<- function(limiares, tipo="mcc"){
  return ((limiares %>%
             filter(get(tipo) == max(get(tipo))))$pre %>%
            max())
}


calcular_limiar_todos_modelos <- function(modelos){
  aval_list <- list()
  modelos <- modelos@models$especie
  for (alg in modelos){
    for (modelo in alg){
      if (length(modelo@evaluation)<=0){
        aval_list[modelo@mID] <- 0
      }
      else {
        aval_list[modelo@mID] <- calcular_limiares_modelo(modelo@evaluation$test.dep@observed,
                                                          modelo@evaluation$test.dep@predicted) %>%
          calcular_limiar_metrica("mcc")
      }
    }
  }
  
  return (aval_list %>% map_df(~ as.data.frame(.), .id="modelID") %>% rename(threshold_MCC="."))
}

select_comp_pca_retain <- function(factoMinerObject){
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

my_fviz_screeplot <- function(aPca){
  retainComp <- select_comp_pca_retain(aPca)
  aPlot <- fviz_screeplot(aPca, addlabels = TRUE) + 
    geom_vline(xintercept=unlist(retainComp), linetype="dashed", color = c("red", "blue", "green"), show.legend = T) +
    annotate("text", label="Broken Stick", x=retainComp$broken.stick, y=20, angle=45) + 
    annotate("text", label="Kaiser Mean", x=retainComp$kaiser.mean, y=25, angle=45) +
    annotate("text", label="Horn P", x=retainComp$horn.p, y=30, angle=45) 
  return(aPlot)
}

my_PCA <- function(data, nrEixos, pcaName){
  aPca <- data %>% PCA(scale.unit = T, graph = F, ncp = nrEixos)
  
  # Loadings (i.e. standard coordinates) are not given by FactoMineR's methods. They return principal coordinates.
  # You can calculate them by dividing variables' coordinates on a dimension by this dimension's eigenvalue's square root.
  aPca$var$loadings <- aPca$var$coord %>% 
    sweep(2, aPca$eig[1:ncol(aPca$var$coord),1] %>% sqrt(), FUN="/") %>%
    as.data.frame()
  if (!is.na(pcaName)){
    colnames(aPca$var$coord) <- paste(colnames(aPca$var$coord), pcaName) %>% make_clean_names()
    colnames(aPca$var$cos2) <- paste(colnames(aPca$var$cos2), pcaName) %>% make_clean_names()
    colnames(aPca$var$cor) <- paste(colnames(aPca$var$cor), pcaName) %>% make_clean_names()
    colnames(aPca$var$contrib) <- paste(colnames(aPca$var$contrib), pcaName) %>% make_clean_names()
    colnames(aPca$var$loadings) <- paste(colnames(aPca$var$loadings), pcaName) %>% make_clean_names()

    colnames(aPca$ind$coord) <- paste(colnames(aPca$ind$coord), pcaName) %>% make_clean_names()
    colnames(aPca$ind$cos2) <- paste(colnames(aPca$ind$cos2), pcaName) %>% make_clean_names()
    colnames(aPca$ind$contrib) <- paste(colnames(aPca$ind$contrib), pcaName) %>% make_clean_names()
  }
  
  return(aPca)
}

my_fiz_contrib_scree <- function(aPca){
  aContribPlot <- fviz_pca_var(aPca, col.var="contrib", 
                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                  legend.title = "Contrib",
                                  repel = TRUE # 
  ) 
  
  aScreePlot <- my_fviz_screeplot(aPca) 
  
  aScreePlot + aContribPlot
}

my_fviz_corr <- function(aPca){
  my_ggcorrplot(aPca$var$contrib, 
                method = "circle", 
                tl.col="black", 
                tl.srt=45, 
                lab=T,
                legend.title = "Contrib",
                lab_size = 2) 
}

my_fiz_contrib_dims <- function(aPca){
  plotList <- list()
  nrEixos <- ncol(aPca$var$coord)
  for (i in 1:nrEixos){
    plotList[[i]] <- aPca %>% fviz_contrib(
      title = colnames(aPca$var$contrib)[i],
      choice = "var",
      axes = i
    )
  }
  
  plotList[[i+1]] <- aPca %>% fviz_contrib(
    title= paste("All", nrEixos, "dims."), 
    choice = "var", 
    axes=1:nrEixos
  )
  wrap_plots(plotList)
}


my_fiz_cos2_dims <- function(aPca){
  plotList <- list()
  nrEixos <- ncol(aPca$var$cos2)
  for (i in 1:nrEixos){
    plotList[[i]] <- aPca %>% fviz_cos2(
        title = colnames(aPca$var$cos2)[i], 
        choice = "var", axes=i
      ) + 
      ylab("Cos2")
  }
  
  plotList[[i+1]] <- aPca %>% fviz_cos2(
      title= paste("All", nrEixos, "dims."), 
      choice = "var", 
      axes=1:nrEixos
      ) + 
    ylab("Cos2")
  
  wrap_plots(plotList)
}

my_fviz_pca_var <- function(aPca, ...){
  fviz_pca_var(aPca, gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE, ...)
}


knit_chunk_as_child <- function(){
  if (is.null(knitr::all_labels())){
    cat("You must knit this document!", sep="\n")
    return()
  }
  if (knitr::opts_current$get("label") %>% str_detect("chunk") %>% any()){
    cat("This chunk must not be include as a child document!", sep="\n")
    return()
  }
  
  current_folder <- knitr::opts_current$get("label") %>%
    discard(~ (.) %>% str_detect("chunk") %>% any()) %>% 
    keep(~ (.) %>% here() %>% is_dir()) %>% 
    here() %>% 
    path() %>% 
    dir_info(type = "directory") %>%
    filter(path %>% str_detect("output_data")) %>%
    select(path, modification_time)
  
  if (nrow(current_folder)==0){
    return(
      here(knitr::opts_current$get("label"),paste0(knitr::opts_current$get("label"), ".Rmd")) %>%
        knitr::knit_child(quiet = TRUE) %>%
        cat(sep = '\n')
    )
  } 
  
  all_folders_details <- knitr::all_labels() %>% 
    discard(~ (.) %>% str_detect("chunk") %>% any()) %>%
    keep(~ (.) %>% here() %>% is_dir()) %>% 
    here() %>% 
    path() %>% 
    dir_info(type = "directory") %>%
    filter(path %>% str_detect("output_data")) %>%
    select(path, modification_time) %>%
    rownames_to_column(var = "id") %>%
    mutate(time_dif = as_datetime(.$modification_time) - lag(as_datetime(.$modification_time))) %>%
    replace(is.na(.), 0)
  
  
  is_newest <- all_folders_details %>%
     filter(path == current_folder$path) %>%
     select(time_dif) %>%
     pull() %>%
    {(.<0)}
  
  if (is_newest){
    return(
      here(knitr::opts_current$get("label"),paste0(knitr::opts_current$get("label"), ".Rmd")) %>%
        knitr::knit_child(quiet = TRUE) %>%
        cat(sep = '\n')
    )
  }    
  
  is_last_one <- all_folders_details %>%
    filter(path == current_folder$path) %>%
    select(id) %>%
    pull() %>%
    {(. == all_folders_details$id %>% last())}
  
  if (is_last_one)
    all_folders_details %>%
      datatable(options = list(pageLength = 10, scrollX=T))
}


render_chunk_as_html <- function(){
  current_chunk <- knitr::opts_current$get("label")
  all_chunks <- knitr::all_labels()

  if (is.null(all_chunks)){
    cat("You must knit this document!", sep="\n")
    return()
  }
  if (current_chunk %>% str_detect("chunk") %>% any()){
    cat("This chunk must not be include as a child document!", sep="\n")
    return()
  }

  current_folder <- current_chunk %>%
    discard(~ (.) %>% str_detect("chunk") %>% any()) %>% 
    keep(~ (.) %>% here() %>% is_dir()) %>% 
    here() %>% 
    path() %>% 
    dir_info(type = "directory") %>%
    filter(path %>% str_detect("output_data")) %>%
    select(path, modification_time)
  
  if (nrow(current_folder)==0){
    path <- callr::r(
      function(...) rmarkdown::render(...),
      args = list(
        here(current_chunk,paste0(current_chunk, ".Rmd")), 
        quiet=T, 
        envir = new.env()
      )
    )
  }
  current_folder <- current_chunk %>%
    discard(~ (.) %>% str_detect("chunk") %>% any()) %>% 
    keep(~ (.) %>% here() %>% is_dir()) %>% 
    here() %>% 
    path() %>% 
    dir_info(type = "directory") %>%
    filter(path %>% str_detect("output_data")) %>%
    select(path, modification_time)
  
  all_folders_details <- all_chunks %>% 
    discard(~ (.) %>% str_detect("chunk") %>% any()) %>%
    keep(~ (.) %>% here() %>% is_dir()) %>% 
    here() %>% 
    path() %>% 
    dir_info(type = "directory") %>%
    filter(path %>% str_detect("output_data")) %>%
    select(path, modification_time) %>%
    rownames_to_column(var = "id") %>%
    mutate(time_dif = as_datetime(.$modification_time) - lag(as_datetime(.$modification_time))) %>%
    replace(is.na(.), 0)
  
  is_newest <- all_folders_details %>%
    filter(path == current_folder$path) %>%
    select(time_dif) %>%
    pull() %>%
    {(.<0)}
  
  if (is_newest){
    path <- callr::r(
      function(...) rmarkdown::render(...),
      args = list(
        here(current_chunk,paste0(current_chunk, ".Rmd")), 
        quiet=T, 
        envir = new.env()
      )
    )
  }      
  
  html_src <- current_chunk %>%
    discard(~ (.) %>% str_detect("chunk") %>% any()) %>% 
    keep(~ (.) %>% here() %>% is_dir()) %>% 
    here() %>% 
    path() %>% dir_ls(glob = "*.html") %>% 
    path_rel()
 
  is_last_one <- all_folders_details %>%
    filter(path == current_folder$path) %>%
    select(id) %>%
    pull() %>%
    {(. == all_folders_details$id %>% last())}
  
  if (is_last_one){
    aTable <- all_folders_details %>%
    datatable(filter="none", list(dom = 't'))
    return(
      htmltools::tagList(list(
        htmltools::tags$iframe(
          title = current_chunk %>% make_clean_names(case = "title") ,
          src = html_src,
          width="100%",
          height="900"
        ),
        htmltools::tags$div(
          aTable
        )
      ))
    )
  } else {
    return(
      htmltools::tags$iframe(
        title = current_chunk %>% make_clean_names(case = "title") ,
        src = html_src,
        width="100%",
        height="900"
      )
    )
  }
}

fortify_join <- function(shp){
  shp_temp <- fortify(shp)
  shp_temp$id <- (shp_temp$id %>% as.integer()) + 1
  
  shp_temp <- shp_temp %>% 
    left_join(shp@data %>% rowid_to_column("id"))
  
  return(shp_temp)
}


