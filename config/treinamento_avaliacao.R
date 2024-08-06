tsne_plot <- function(df, df_bg, sp_names){
  sp_names <- sp_names %>% 
    to_snake_case() %>% 
    abbreviate(minlength = 10) %>%
    unname()
  
  result <- lapply(sp_names, function(sp){
    df_pred <- df %>%
      filter(.[sp]==1) %>% 
      bind_rows(df_bg[[sp]]) %>%
      replace(is.na(.), 0) %>% 
      select(sp,all_of(names(df_bg[[sp]])))
    
    # Empiricamente um bom valor de perplexity é: P ~ N^1/2
    perp = round((nrow(df_pred) ^ (1/2)), digits=0)
    df_bg_tsne <- df_pred[!duplicated(df_pred), ]
  
    tsne_bg <- Rtsne(as.matrix(df_bg_tsne), perplexity = perp)
    
    esp_colors <- df_pred %>% 
      unique() %>% 
      select(sp) %>% 
      pull()
    
    tsne_result <- qplot(
      tsne_bg$Y[,1], 
      tsne_bg$Y[,2],
      main = sp,
      geom = "point", 
      colour = esp_colors
    )
    return(tsne_result)
  })
  return(result)
}

fit_data <- function(df_pa, df_var, df_bg){
  #if(class(df_var)=='list'){
  #  df_var <- as.data.frame(df_var)
  #}
  #if(class(df_bg)=='list'){
  #  df_bg <- as.data.frame(df_bg)
  #}
  if('geometry' %in% colnames(df_pa)){
    df_pa %<>% select(-c('geometry'))
  }
  if('cell_id' %in% colnames(df_bg[[1]])){
    df_bg %<>% lapply(function(x) { x["cell_id"] <- NULL; x })
  }
  if('cell_id' %in% colnames(df_var[[1]])){
    df_var %<>% lapply(function(x) { x["cell_id"] <- NULL; x })
  }
  
  fitted_data <- list()
  for (sp in colnames(df_pa)){
    if(class(df_var)=='list'){
      df_var2 <- df_var[[sp]] %>% select(colnames(df_bg[[sp]]))
    } else {
      df_var2 <- df_var %>% select(colnames(df_bg[[sp]]))
    }
    
    predict_data <- sdmData(
      sp %>% 
        paste0(collapse = " + ") %>% 
        paste(".", sep=" ~ ") %>% 
        as.formula(), 
      train = df_pa %>% 
        select(all_of(sp)) %>% 
        bind_cols(df_var2) %>% 
        filter(.[sp]==1),
      bg = df_bg[[sp]] %>% 
        bind_cols(rep(0,nrow(df_bg[[sp]])) %>% list() %>% set_names(sp))
    )
    fitted_data <- fitted_data %>% 
      append(list(predict_data) %>% set_names(sp))
  }
  names(fitted_data) <- colnames(df_pa)
  return(fitted_data)
}

train_models <- function(df_pa, fitted_data, pred_methods, n_exec, n_folds){
  trained_models <- list()
  for (sp in colnames(df_pa)){
    suppressWarnings(
      a_model <-
        sdm(
          data = fitted_data[[sp]],
          methods = pred_methods,
          n = n_exec,
          cv.folds = n_folds,
          replication = "cv"
        )
    )
    trained_models <- trained_models %>% 
      append(list(a_model) %>% set_names(sp))
  }
  return(trained_models)
}

train_models_to_folder <- function(df_pa, fitted_data, pred_methods, n_exec, n_folds, output_folder){
  for (sp in colnames(df_pa)){
    output_folder_tmp <- output_folder %>%
      path(sp)
    if (!folder_models %>% path(paste0(sp, ".sdm")) %>% file_exists()){
      dir_create(output_folder_tmp)
    
      suppressWarnings(
        m <-
          sdm(
            data = fitted_data[[sp]],
            methods = pred_methods,
            n = n_exec,
            cv.folds = n_folds,
            replication = "cv"
          )
      )

      #setting outputs
      eval <- getEvaluation(m, wtest = "indep.test",stat = c('threshold', 'AUC', 'COR', 'sensitivity', 'specificity', 'TSS', 'Kappa'), opt="max(se+sp)") # get evaluation for models
      varimp <- getVarImp(m) 
      varimp <- varimp@varImportanceMean # get variables importance
      varimp_auc <- varimp$AUCtest
      varimp_correlation <- varimp$corTest

      #### Save Some Outputs ####
      if(nrow(eval) >= 1) write.csv(eval, paste0(output_folder_tmp,"/",sp,"_EvalMetrics.csv"))
      if(nrow(varimp_auc) >= 1) write.csv(varimp_auc, paste0(output_folder_tmp,"/",sp,"_VariableImportanceAUC.csv"))
      if(nrow(varimp_correlation) >= 1) write.csv(varimp_correlation, paste0(output_folder_tmp,"/",sp,"_VariableImportanceCorrelation.csv"))

      m %>% 
        saveRDS(output_folder_tmp %>% path(paste0(sp, ".sdm")))
    }
  }
}

sp_thresh <- function(t_models, thr_criteria){
  t_models %>% 
    #map_dfr(~ 
    getEvaluation(
      ., 
      stat=c(
        'threshold', 'AUC', 'COR', 'sensitivity', 'specificity', 'TSS', 'Kappa'
      ), 
      opt=thr_criteria
    )
  #, 
  #.id = "species"
  #)
}

#sp_name = spp_names
#folder = folder_models
#thr_criteria = thresh_criteria
sp_thresh_from_folder <- function(sp_name, folder, thr_criteria){
  ths <- list()
  for (sp in sp_name) {
    m <- sp %>%
      sp_model_from_folder(folder)
    if('glmnet' %in% m[[sp]]@setting@methods){
      m <- sp %>%
        sp_model_from_folder(folder)
      ids <- m[[sp]] %>% getModelInfo()
      ids <- ids$modelID[ids$method=='glmnet']
      df <- sp %>%
        sp_model_from_folder(folder) %>%
        sp_thresh(thr_criteria)
      if(any(df[ids,'threshold']<0)){
        df[ids,'threshold'] <- boot::inv.logit(df[ids,'threshold'])
      }
      ths[[sp]] <- df %>%
        mutate(species=sp)
    } else {
      ths[[sp]] <- m[[sp]] %>%
        sp_thresh(thr_criteria) %>%
        mutate(species=sp)
    }
  }
  return(ths)
}

thresh_mean <- function(t_models, df_thr, stats){ #, pred_methods){
  t_models %>% 
    map_dfr(.f=function(x){rename(x@run.info, species="species") %>% 
                           mutate_if(.,is.factor, as.character)}, .id="species") %>%
    #filter(method %in% pred_methods) %>%
    select(species, modelID, method) %>%
    inner_join(df_thr %>% select(species, modelID, (!!sym(stats))), by=c("species","modelID")) %>%
    group_by(species, method) %>% 
    summarize(stats = mean((!!sym(stats))), .groups = 'rowwise') %>%
    return()
}

validation_metrics_from_folder <- function(sp_name, folder, df_thr, stats='threshold', th=NULL){
  thr_mean <- list()
  for(sp in sp_name){
    thr_mean[[sp]] <- sp %>%
      sp_model_from_folder(folder) %>%
      unlist() %>%
      thresh_mean(df_thr[[sp]], stats)
  }
  if(!is.null(th)){
    thr_mean <- lapply(thr_mean, function(x){x[x[,"stats"]>th,]})
  }
  
  return(thr_mean)
} 

sp_thresh_mean_from_folder <- function(sp_name, folder, df_thr){
  thr_mean <- list()
  for(sp in sp_name){
    thr_mean[[sp]] <- sp %>%
        sp_model_from_folder(folder) %>%
        unlist() %>%
        thresh_mean(df_thr[[sp]])
  }
  return(thr_mean)
} 

sp_model_from_folder <- function(sp_name, folder){
  if(nchar(sp_name)>11){
    sp_name <- sp_name %>% 
      to_snake_case() %>% 
      abbreviate(minlength = 10) %>% 
      unname()
  }
  folder %>% 
    #dir_ls(glob = paste0("*", sp_name, "*.sdm"), recurse = T) %>%
    dir_ls(glob = paste0("*", sp_name, ".sdm"), recurse = T) %>%
    read.sdm() %>%
    list() %>%
    set_names(sp_name)
}

sp_predictions_from_folder <- function(sp_name, folder){
  sp_name <- sp_name %>% 
    to_snake_case() %>% 
    abbreviate(minlength = 10) %>% 
    unname()

  suppressMessages(
    folder %>% 
      dir_ls(glob = paste0("*", sp_name, "*.csv"), recurse = T) %>% 
      purrr::map(vroom) %>% 
      set_names(path_file) %>% 
      set_names(path_ext_remove)
  )
}

predict_sp_to_folder <- function(df_pa, df_var, t_models, shp, pred_methods, thr_criteria, output_folder){
  for (sp in colnames(df_pa)){
    df_pred_freq <- df_var$current %>% 
      DRE_predict(t_models[[sp]], pred_methods, thr_criteria)
    # <- DRE_predict(modelos_treinados, df_predicao_var, alg_predicao, criterio_thres, thresholds_modelos_mcc)
    
    colnames(df_pred_freq) <- pred_methods
    df_pred_freq <- df_pred_freq %>% 
      mutate(consensus=rowMeans(.))
    
    # Consenso usando presenca/ausencia
    df_pred_pres_aus <- df_pred_freq %>%
      select(-consensus) %>%
      mutate_all(~ ifelse(. <= 0.5, 0, 1)) %>%
      mutate(consensus = rowMeans(.) %>% round(., digits = 0))
    
    df_pred_freq <- df_pa %>% 
      select(sp %>% all_of()) %>% 
      cbind(df_pred_freq)
    
    df_pred_pres_aus <- df_pa %>% 
      select(sp %>% all_of()) %>%
      cbind(df_pred_pres_aus)
    
    shp_tmp <- shp
    
    shp_tmp@data <- df_pred_freq 
    shp_tmp %>%
      writeOGR(
        dsn = output_folder %>% fs::path(sp), 
        layer=paste0(sp, "_freq"), 
        driver="ESRI Shapefile"
      )
    
    shp_tmp@data <- df_pred_pres_aus
    shp_tmp %>% 
      writeOGR(
        dsn = output_folder %>% fs::path(sp), 
        layer=paste0(sp, "_pa"), 
        driver="ESRI Shapefile"
      )
  }
}

# TODO: Qual o comportamento da predição para a pasta, sobreescrever sempre ou continuar de onde parou? Qual o comportamento para as outras funções sobreescrever ou continuar?

# Projetar modelos nos cenários



#df_pa = df_pa
#scenarios_list =   projection_data
#models_folder =   folder_models
#pred_methods =   model_selection
#thr_criteria =   thresh_criteria
#output_folder =   directory_projections
#
#thresholds_models_means
#sp = colnames(df_pa)[1]
#scenario_name = names(scenarios_list)[1]

predict_to_folder <- function(df_pa, scenarios_list, models_folder, pred_methods, thr_criteria, output_folder, thresholds_models_means){
  if ('geometry' %in% colnames(df_pa)){
    spp_names <- colnames(df_pa)[-grep("geometry", colnames(df_pa))]
  } else {
    spp_names <- colnames(df_pa)
  }
  
  names(thresholds_models_means) <- lapply(thresholds_models_means, function(x){unique(x$species)})
  selected_models <- lapply(pred_methods, function(x) x$method)
  
  for (sp in spp_names){
    if(selected_models[[sp]] %>% length() == 0){
      print(paste0("No model passing the threshold for ",sp))
    } else {
      t_models <- sp %>%
        sp_model_from_folder(models_folder)
      if (scenarios_list %>% is.data.frame()){
        scenarios_list <- scenarios_list %>%
          list() %>%
          set_names("scenarios")
      }
      
      th <- thresholds_models_means[[sp]] %>% arrange(method)
      th <- th$stats
      for (scenario_name in names(scenarios_list)){
        folder_tmp <- output_folder %>% path(sp) %>% path(scenario_name)
        file_tmp <- scenario_name
        
        if (!dir_exists(folder_tmp)){
          print(paste0('Predicting ',scenario_name))
          df_pred_freq <- scenarios_list %>% 
            pluck(scenario_name) %>% 
            #pluck(sp) %>%
            DRE_predict(t_models %>% pluck(sp), selected_models[[sp]], thr_criteria)
          
          # Consenso
          #df_pred_freq <- df_pred_freq %>% 
          #  mutate(consenso=rowMeans(.))
          consensus <- df_pred_freq %>% rowMeans()
          
          # Weighted mean (AUC)
          ids <- getEvaluation(t_models[[sp]], stat=c('TSS','AUC','threshold'), wtest="dep.test", opt=thr_criteria)$modelID
          df <- cbind(algo=getModelInfo(t_models[[sp]])$method[ids],
                      getEvaluation(t_models[[sp]], stat=c('TSS','AUC','threshold'), wtest="dep.test", opt=thr_criteria))
          w <- aggregate(df$AUC, by=list(df$algo), FUN=mean)$x
          wmean_AUC <- apply(df_pred_freq, 1, function(x) sum(x * w)/(df$algo %>% unique() %>% length()))
          
          # Weighted mean (TSS)
          w <- aggregate(df$TSS, by=list(df$algo), FUN=mean)$x
          wmean_TSS <- apply(df_pred_freq, 1, function(x) sum(x * w)/(df$algo %>% unique() %>% length()))
          
          
          # include ensembles in data.frame
          df_pred_freq$consensus <- consensus
          df_pred_freq$wmean_AUC <- wmean_AUC
          df_pred_freq$wmean_TSS <- wmean_TSS
          
          # Consenso usando presenca/ausencia
          df_pred_pres_aus <- df_pred_freq %>%
            select(-c(consensus, wmean_AUC, wmean_TSS)) %>% 
            select(sort(names(.))) %>%
            mutate_all(funs(ifelse(. > th, 1, 0))) %>% 
            mutate(consensus = rowMeans(.) %>% round(., digits = 0))
          
          # Ensembles:
          df_pred_freq <- rep(sp, nrow(df_pred_freq)) %>% 
            as.data.frame() %>% 
            rename("species"=".") %>%
            bind_cols(rep(scenario_name, nrow(df_pred_freq)) %>% as.data.frame()) %>%
            rename("scenario"=".") %>%
            bind_cols(df_pred_freq)
          
          df_pred_pres_aus <- rep(sp, nrow(df_pred_pres_aus)) %>% 
            as.data.frame() %>% 
            rename("species"=".") %>%
            bind_cols(rep(scenario_name, nrow(df_pred_pres_aus)) %>% as.data.frame()) %>%
            rename("scenario"=".") %>%
            bind_cols(df_pred_pres_aus)
          
          
          folder_tmp %>% 
            dir_create(recurse = T)
          
          colnames(df_pred_freq) <- df_pred_freq %>% 
            colnames() %>%
            to_snake_case() %>%
            abbreviate(minlength = 10) %>%
            unname()
          
          suppressMessages(
            df_pred_freq %>%
              vroom_write(
                folder_tmp %>% 
                  path(file_tmp) %>% 
                  paste0("_freq.csv"),
                delim = ";"
              )
          )
          
          #colnames(df_pred_pres_aus) <- df_pred_freq %>% 
          #  colnames() %>%
          #  to_snake_case() %>%
          #  abbreviate(minlength = 10) %>%
          #  unname()
          
          suppressMessages(
            df_pred_pres_aus %>%
              vroom_write(
                folder_tmp %>% 
                  path(file_tmp) %>% 
                  paste0("_pa.csv"),
                delim = ";"  
              )
          )
        }
      }
    }
  }
}

model_failures <- function(fitted_data, models_folder){
  fails <- list()
  for(i in length(names(fitted_data))){
      m <- sp_model_from_folder(names(fitted_data[[1]]), models_folder) %>%
                pluck(1) %>%
                slot("run.info") %>% 
                select(species, method, success) %>%
                filter(success == F) %>%
                unique() %>%
                nrow()
      fails[[i]] <- m
  }
  return(m)
}



# TODO: Esse método precisa ser refatorado: quebrado em diversos outros métodos. Também precisa ser melhorado com novas funções.
#gerarPseudoAusenciasDRE <- function (shp_pa, df_var, especies, metodo="random", cluster_m="k-means"){
#  especies <- especies %>% 
#    to_snake_case() %>% 
#    abbreviate(minlength = 10) %>%
#    unname()
#  
#  df_pa <- shp_pa@data %>% select(all_of(especies))
#  .add_grupos <- function(a_df, a_df_p, percent_grupos=NA, nro_grupos=NA, cluster_m="k-means"){
#    if (cluster_m=="k-means"){
#      if (is.na(nro_grupos)) {
#        nro_grupos <- round(nrow(a_df_p)/percent_grupos)  
#      }
#      set.seed(nro_grupos)
#      a_df$grupo <- a_df %>% 
#        select(-any_of("grupo")) %>%
#        kmeans(nro_grupos) %>%
#        pluck("cluster") %>%
#        as.factor()
#    }
#    else {
#      # Não está funcionando
#      if (is.na(nro_grupos)){
#        nro_grupos <- ceiling(nrow(a_df_p)/percent_grupos)  
#      } 
#      a_df$grupo <- a_df %>% 
#        as.matrix() %>% 
#        Xmeans(nrow(a_df)/100, nthread = 3, min.clust.size = nro_grupos) %>%
#        pluck("cluster")
#    }
#    return(a_df)
#  }
#  
#  .get_summary_from <- function(a_df_var, a_list_pa){
#    p <- a_list_pa %>% keep(~ .==1)
#    grupos <- a_df_var %>% 
#      bind_cols(especie=a_list_pa) %>%
#      group_by(grupo) %>%
#      summarise(
#        presencas = sum(especie), 
#        quant_celulas= n(), 
#        dens_presenca=sum(especie)/n(), 
#        presenca_rel=sum(especie)/nrow(a_df_var),
#        quant_cel_bg= round(n() * length(p) / nrow(a_df_var)) 
#      ) #%>%
#    #filter(dens_presenca <= mean(.$dens_presenca))
#    return(grupos)
#  }
#  
#  if (metodo == "cluster") {
#    backgrounds = list()
#    for (sp in colnames(df_pa)){
#      df_p <- df_pa %>% select(sp %>% all_of()) %>% filter(.[sp]==1)
#      df_var <- df_var %>% .add_grupos(df_p, percent_grupos = 10, cluster_m = cluster_m) # k = 10% do numero de pres.
#      grupos <- df_var %>% .get_summary_from(df_pa[[sp]])
#      
#      df_bg <- df_var %>%
#        bind_cols(especie=df_pa[[sp]]) %>%
#        filter(especie==0) %>%
#        filter(grupo %in% grupos$grupo) %>%
#        group_by(grupo) %>% 
#        select(-especie) %>%
#        nest() %>%
#        ungroup() %>%
#        mutate(n=grupos$quant_cel_bg) 
#      
#      df_bg_tmp <- try(
#        df_bg %>% 
#          mutate(samp = map2(data, n, sample_n)) %>% 
#          select(-data, -n, -grupo) %>%
#          unnest(samp)
#      )
#      if ("try-error" %in% class(df_bg_tmp)){
#        df_bg <- df_bg %>% 
#          mutate(samp = map2(data, n, ~ sample_n(.x, .y, replace = T))) %>% 
#          select(-data, -n, -grupo) %>%
#          unnest(samp)
#      }
#      else {
#        df_bg <- df_bg_tmp
#      }
#      
#      backgrounds <- backgrounds %>% 
#        append(list(df_bg) %>% set_names(sp))
#    }
#    return(backgrounds)
#    # # visualização dos pontos de background
#    # grid_temp <- grid
#    # grid_temp@data <- df_tmp
#    # grid_temp@data$id <- rownames(grid_temp@data)
#    # grid_temp <- fortify(grid_temp, region = "id") %>% left_join(grid_temp@data)
#    # 
#    # ggplot(data = fortify(grid_temp)) +
#    #   aes(x = long,
#    #       y = lat,
#    #       group = group,
#    #       text = grupo) +
#    #   geom_polygon(aes(fill = as.factor(grupo))) +
#    #   #geom_point(
#    #   #  data = df_bg %>% filter(species == 1),
#    #   #  aes(x = long, y = lat, colour = species),
#    #   #  size = 0.5
#    #   #) +
#    #   theme(
#    #     axis.text.x = element_blank(),
#    #     axis.text.y = element_blank(),
#    #     axis.ticks.x = element_blank(),
#    #     axis.ticks.y = element_blank()
#    #   ) +
#    #   coord_equal()
#    
#  } else if (metodo=="dre_area" || metodo=="dre_dens"){
#    quant_max <- 0
#    num_grupo_quant_max <- 2
#    dens_acum <- 0
#    num_grupo_dens_acum <- 2
#    backgrounds = list()
#    for (sp in colnames(df_pa)){
#      df_p <- df_pa %>% select(sp %>% all_of()) %>% filter(.[sp]==1)
#      
#      numero_grupos <- trunc(nrow(df_p))
#      for (i in 2:numero_grupos){
#        df_var <- df_var %>% .add_grupos(df_p, nro_grupos=i) # k = 1 do numero de pres.
#        grupos <- df_var %>% .get_summary_from(df_pa[[sp]])
#        
#        grupos_area <- grupos %>%
#          filter(dens_presenca <=mean(grupos$dens_presenca))
#        
#        grupos_dens <- grupos %>%
#          filter(dens_presenca > mean(grupos$dens_presenca))
#        
#        if (sum(grupos_area$quant_celulas) > quant_max){
#          quant_max <- sum(grupos_area$quant_celulas)
#          num_grupo_quant_max <- i
#        }
#        if (sum(grupos_dens$presencas) / sum(grupos_dens$quant_celulas) > dens_acum){
#          dens_acum <- sum(grupos_dens$presencas) / sum(grupos_dens$quant_celulas)
#          num_grupo_dens_acum <- i
#        }
#      }
#      
#      if (metodo=="dre_area"){
#        numero_grupos <- num_grupo_quant_max  
#      } else {
#        numero_grupos <- num_grupo_dens_acum
#      }
#      
#      df_var <- df_var %>% .add_grupos(df_p, nro_grupos=numero_grupos) # k = 1 do numero de pres.
#      grupos <- df_var %>% .get_summary_from(df_pa[[sp]])
#      
#      if (metodo=="dre_area"){
#        grupos <- grupos %>%
#          filter(dens_presenca <=mean(grupos$dens_presenca))
#      }
#      else {
#        grupos <- grupos %>%
#          filter(dens_presenca > mean(grupos$dens_presenca))
#      }
#      grupos <- grupos %>% 
#        mutate(quant_cel_bg=round(nrow(df_p) * grupos$quant_celulas / sum(grupos$quant_celulas)))
#      
#      
#      df_bg <- df_var %>%
#        bind_cols(especie=df_pa[[sp]]) %>%
#        filter(especie==0) %>%
#        filter(grupo %in% grupos$grupo) %>%
#        group_by(grupo) %>% 
#        select(-especie) %>%
#        nest() %>%
#        ungroup() %>%
#        arrange(grupo) %>%
#        mutate(n=grupos$quant_cel_bg)
#      
#      df_bg_tmp <- try(
#        df_bg %>% 
#          mutate(samp = map2(data, n, sample_n)) %>% 
#          select(-data, -n, -grupo) %>%
#          unnest(samp)
#      )
#      
#      if ("try-error" %in% class(df_bg_tmp)){
#        df_bg <- df_bg %>% 
#          mutate(samp = map2(data, n, ~ sample_n(.x, .y, replace = T))) %>% 
#          select(-data, -n, -grupo) %>%
#          unnest(samp)
#      }
#      else {
#        df_bg <- df_bg_tmp
#      }
#      
#      backgrounds <- backgrounds %>% 
#        append(list(df_bg) %>% set_names(sp))
#    }
#    return(backgrounds)
#  } 
#  else {
#    backgrounds = list()
#    for (sp in colnames(df_pa)){
#      df_p <- df_pa %>% select(sp %>% all_of()) %>% filter(.[sp]==1)
#      
#      df_bg <- df_var %>%
#        bind_cols(especie=df_pa[[sp]]) %>%
#        filter(especie==0) %>%
#        sample_n(nrow(df_p)) %>%
#        select(-especie)
#      
#      backgrounds <- backgrounds %>% 
#        append(list(df_bg) %>% set_names(sp))
#    }
#    return(backgrounds)
#  }
#}

#shp_pa=shp_matrix_pa
#df_var=df_var_preditors
#especies=spp_names
#method="cluster"
#cluster_m="k-means"
#n_pa=1
#output_folder=here("output_data/models")

#sp=colnames(df_pa)[1]

pseudoabsences <- function (shp_pa, df_var, especies, method="envelope", cluster_m="kmeans", n_pa=1, output_folder=here("output_data/models")){
  especies <- especies %>% 
    to_snake_case() %>% 
    abbreviate(minlength = 10) %>%
    unname()
  
  if('geometry' %in% colnames(shp_pa)){
    df_pa <- shp_pa %>% as.data.frame() %>% select(-c('geometry'))
  } else {
    if('geom' %in% colnames(shp_pa)){
      df_pa <- shp_pa %>% as.data.frame() %>% select(-c('geom'))
    } else {
      df_pa <- shp_pa %>% as.data.frame() 
    }
  }
  
  if(class(df_var)=='list'){
    df_var <- df_var[[especies]]
  }
  
  .add_grupos <- function(a_df, a_df_p, percent_grupos=NA, nro_grupos=NA, cluster_m="kmeans"){
    if (cluster_m=="kmeans"){
      if (is.na(nro_grupos)) {
        nro_grupos <- round(nrow(a_df_p)/percent_grupos)  
      }
      set.seed(nro_grupos)
      a_df$grupo <- a_df %>% 
        select(-any_of("grupo")) %>%
        kmeans(nro_grupos) %>%
        pluck("cluster") %>%
        as.factor()
      return(a_df)
    }
    #else {
    #  # Não está funcionando
    #  if (is.na(nro_grupos)){
    #    nro_grupos <- ceiling(nrow(a_df_p)/percent_grupos)  
    #  } 
    #  a_df$grupo <- a_df %>% 
    #    as.matrix() %>% 
    #    Xmeans(nrow(a_df)/100, nthread = 3, min.clust.size = nro_grupos) %>%
    #    pluck("cluster")
    #}
    return(a_df)
  }
  
  .get_summary_from <- function(a_df_var, a_list_pa){
    p <- a_list_pa %>% keep(~ .==1)
    grupos <- a_df_var %>% 
      bind_cols(especie=a_list_pa) %>%
      group_by(grupo) %>%
      summarise(
        presencas = sum(especie), 
        quant_celulas= n(), 
        dens_presenca=sum(especie)/n(), 
        presenca_rel=sum(especie)/nrow(a_df_var),
        quant_cel_bg= round(n() * length(p) / nrow(a_df_var)) 
      ) #%>%
    #filter(dens_presenca <= mean(.$dens_presenca))
    return(grupos)
  }
  
  if (method == "cluster") {
    backgrounds = list()
    for (sp in names(df_pa)){
      if (!file_exists(paste0(output_folder,'/',sp,'/pseudoabsence_',method,'_',sp,'_',n_pa,'.csv'))){
        df_p <- df_pa[[sp]] %>% select(sp %>% all_of()) %>% filter(.[sp]==1)
        print(colnames(df_var))
        df_var2 <- df_var %>% .add_grupos(df_p, percent_grupos = 10, cluster_m = cluster_m) # k = 10% do numero de pres.
        print(colnames(df_var2))
        grupos <- df_var2 %>% .get_summary_from(df_pa[[sp]])
        
        df_bg <- df_var2 %>%
          bind_cols(especie=df_pa[[sp]]) %>%
          filter(especie==0) %>%
          filter(grupo %in% grupos$grupo) %>%
          group_by(grupo) %>% 
          select(-especie) %>%
          nest() %>%
          ungroup() %>%
          mutate(n=grupos$quant_cel_bg) 
        
        df_bg_tmp <- try(
          df_bg %>% 
            mutate(samp = map2(data, n, sample_n)) %>% 
            select(-data, -n, -grupo) %>%
            unnest(samp)
        )
        if ("try-error" %in% class(df_bg_tmp)){
          df_bg <- df_bg %>% 
            mutate(samp = map2(data, n, ~ sample_n(.x, .y, replace = T))) %>% 
            select(-data, -n, -grupo) %>%
            unnest(samp)
        }
        else {
          df_bg <- df_bg_tmp
        }
        
        backgrounds <- backgrounds %>% 
          append(list(df_bg) %>% set_names(sp))
        
        write.csv(df_bg, paste0(output_folder,'/',sp,'/pseudoabsence_',method,'_',sp,'_',n_pa,'.csv'), row.names = F)
      } else {
        backgrounds[[sp]] <- read.csv(paste0(output_folder,'/',sp,'/pseudoabsence_',method,'_',sp,'_',n_pa,'.csv'), sep=',')
      }
    }
    return(backgrounds)
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
    
  }
  
  if (method == "envelope" || method == "bioclim"){
    backgrounds = list()
    for (sp in colnames(df_pa)){
      if(!file.exists(paste0(output_folder,'/',sp,'/pseudoabsence_',method,'_',sp,'_',n_pa,'.csv'))){
        if(!file.exists(paste0(output_folder,'/',sp,'/bioclim_',sp,'.tif'))){
          
          if(!dir.exists(paste0(output_folder,'/',sp))){
            dir.create(paste0(output_folder,'/',sp), recursive=T)
          }
          print(paste0("Building pseudoabsences for ",sp,"..."))
          
          
          
          sp_dat <- cbind(df_pa[,sp],df_var)
          colnames(sp_dat) <- c(sp,colnames(sp_dat)[-1])
          sp_dat <- sp_dat[sp_dat[,sp]==1,]

          d <- sdmData(reformulate(termlabels = colnames(df_var), response = sp),
                       train = sp_dat)

          m <- sdm(~.,
                   data = d,
                   methods='bioclim')

          p <- predict(m, df_var,
                       filename=paste0(output_folder,'/',sp,'/bioclim_',sp,'.tif'))
        } else {
          p <- raster(paste0(output_folder,'/',sp,'/bioclim_',sp,'.tif'))
        }
        p[p[]>0] <- NA
        env2 <- df_var[!is.na(p),]
        
        if(nrow(env2) > nrow(sp_dat)){
          pa <- env2[sample(nrow(env2), nrow(sp_dat)),]
        } else {
          pa <- env2[sample(nrow(env2), nrow(sp_dat), replace = T),]
        }
        
        backgrounds <- backgrounds %>% 
          append(list(pa) %>% set_names(sp))
        
        write.csv(pa, paste0(output_folder,'/',sp,'/pseudoabsence_',method,'_',sp,'_',n_pa,'.csv'), row.names = F)
      } else {
      backgrounds[[sp]] <- read.csv(paste0(output_folder,'/',sp,'/pseudoabsence_',method,'_',sp,'_',n_pa,'.csv'), sep=',')
      }
    }
    return(backgrounds)
  }
  
  if (method == "dre_area" || method == "dre_dens"){
    quant_max <- 0
    num_grupo_quant_max <- 2
    dens_acum <- 0
    num_grupo_dens_acum <- 2
    backgrounds = list()
    for (sp in names(df_pa)){
      df_p <- df_pa[[sp]] %>% select(sp %>% all_of()) %>% filter(.[sp]==1)
      
      numero_grupos <- trunc(nrow(df_p))
      for (i in 2:numero_grupos){
        df_var <- df_var %>% .add_grupos(df_p, nro_grupos=i) # k = 1 do numero de pres.
        grupos <- df_var %>% .get_summary_from(df_pa[[sp]])
        
        grupos_area <- grupos %>%
          filter(dens_presenca <=mean(grupos$dens_presenca))
        
        grupos_dens <- grupos %>%
          filter(dens_presenca > mean(grupos$dens_presenca))
        
        if (sum(grupos_area$quant_celulas) > quant_max){
          quant_max <- sum(grupos_area$quant_celulas)
          num_grupo_quant_max <- i
        }
        if (sum(grupos_dens$presencas) / sum(grupos_dens$quant_celulas) > dens_acum){
          dens_acum <- sum(grupos_dens$presencas) / sum(grupos_dens$quant_celulas)
          num_grupo_dens_acum <- i
        }
      }
      
      if (method=="dre_area"){
        numero_grupos <- num_grupo_quant_max  
      } else {
        numero_grupos <- num_grupo_dens_acum
      }
      
      df_var <- df_var %>% .add_grupos(df_p, nro_grupos=numero_grupos) # k = 1 do numero de pres.
      grupos <- df_var %>% .get_summary_from(df_pa[[sp]])
      
      if (method=="dre_area"){
        grupos <- grupos %>%
          filter(dens_presenca <=mean(grupos$dens_presenca))
      }
      else {
        grupos <- grupos %>%
          filter(dens_presenca > mean(grupos$dens_presenca))
      }
      grupos <- grupos %>% 
        mutate(quant_cel_bg=round(nrow(df_p) * grupos$quant_celulas / sum(grupos$quant_celulas)))
      
      
      df_bg <- df_var %>%
        bind_cols(especie=df_pa[[sp]]) %>%
        filter(especie==0) %>%
        filter(grupo %in% grupos$grupo) %>%
        group_by(grupo) %>% 
        select(-especie) %>%
        nest() %>%
        ungroup() %>%
        arrange(grupo) %>%
        mutate(n=grupos$quant_cel_bg)
      
      df_bg_tmp <- try(
        df_bg %>% 
          mutate(samp = map2(data, n, sample_n)) %>% 
          select(-data, -n, -grupo) %>%
          unnest(samp)
      )
      
      if ("try-error" %in% class(df_bg_tmp)){
        df_bg <- df_bg %>% 
          mutate(samp = map2(data, n, ~ sample_n(.x, .y, replace = T))) %>% 
          select(-data, -n, -grupo) %>%
          unnest(samp)
      }
      else {
        df_bg <- df_bg_tmp
      }
      
      backgrounds <- backgrounds %>% 
        append(list(df_bg) %>% set_names(sp))
    }
    return(backgrounds)
  } 
  
  if (method == "random") {
    
    backgrounds = list()
    
    for (sp in especies){
      output_folder_tmp <- output_folder %>%
        path(sp)
      if (!file_exists(paste0(output_folder,'/',sp,'/pseudoabsence_',method,'_',sp,'_',n_pa,'.csv'))){
        dir_create(output_folder_tmp)
        
        n_pa2 <- df_pa[,sp] %>% sum()
        
        df_bg <- df_var %>%
          bind_cols(especie=df_pa[,sp]) %>%
          filter(especie==0) %>%
          sample_n(n_pa2) %>%
          select(-especie)
        
        backgrounds <- backgrounds %>% 
          append(list(df_bg) %>% set_names(sp))
        
        #backgrounds[[sp]] <- df_bg
        if('geometry' %in% colnames(df_bg)){
          df_bg <- df_bg %>% select(-'geometry')
        }
        
        write.csv(df_bg, paste0(output_folder,'/',sp,'/pseudoabsence_',method,'_',sp,'_',n_pa,'.csv'), row.names = F)
      } else {
        backgrounds[[sp]] <- read.csv(paste0(output_folder,'/',sp,'/pseudoabsence_',method,'_',sp,'_',n_pa,'.csv'), sep=',')
      }
    }
    return(backgrounds)
  }
}
