# Inclui o algoritmo de predição DRE_Mahal no pacote SDM do Naimi ----------------------------

# Basicamente, para inclui um novo método são necessário 6 passos:
# 1 - Declarar uma classe S4 que contém o modelo treinado, o método de ajuste do modelo e o método de predição
# 2 - Uma função para ajuste do modelo (.mahal_m_fit) que retorna um objeto da classe que representa onovo algoritmo de predição 
# 3 - Um novo método para predição deve ser declarado na classe que representa o novo algoritmo de predição
# 4 - Uma lista com as características do novo algoritmo de predição deve ser criado de acordo com as regras criadas pelo Naimi
# 5 - O novo método que faz a predição do novo algoritmo de predição deve ser alocado em um ".GlobalEnv" para ser acessível pelo pacote SDM do Naimi
# 6 - Finalmente, o novo método de predição do novo algoritmo é adicionado ao pacote SDM do Naimi

.prepare_data <- function(formula, data, ...){
  prep_data <- list()
  prep_data$varnames <- all.vars(formula)
  prep_data$nsp <- deparse(formula[[2]])
  prep_data$all_data <- data %>% 
    select_if(~ !is.factor(.))

  prep_data$presences <- prep_data$all_data %>% 
    filter(.[[prep_data$nsp]]==1) %>%
    select(-prep_data$nsp)
  
  prep_data$absences <- prep_data$all_data %>% 
    filter(.[[prep_data$nsp]]==0) %>%
    select(-prep_data$nsp)

  if (ncol(prep_data$presences) < 2) 
    stop("At least two continous variables are needed to fit the model!")
  
  return(prep_data)
}

# Passo 1
classe_predicao <- setClass("Mahal_m", 
                            representation(
                              features="character",
                              cov="matrix",
                              presence="data.frame"
                            )
)

# Passo 2 
.mahal_m_fit <- function(formula, data, ...) {
  prep_data <- .prepare_data(formula, data)
  
  mahal_m <- new("Mahal_m")
  mahal_m@features <- colnames(prep_data$presences)
  mahal_m@presence <- prep_data$presences
  mahal_m@cov <- solve(cov(mahal_m@presence), tol=1e-20)
  
  mahal_m
}

# Passo 3
# https://rpubs.com/melinatarituba/356726
# http://w3.ufsm.br/adriano/livro/Caderno%20dedatico%20multivariada%20-%20LIVRO%20FINAL%201.pdf
# https://docs.ufpr.br/~soniaisoldi/ce090/CE076AM_2010.pdf
# https://www.machinelearningplus.com/statistics/mahalanobis-distance/
# http://eric.univ-lyon2.fr/~ricco/tanagra/fichiers/en_Tanagra_Calcul_P_Value.pdf
metodo_predicao <- setMethod("predict", signature(object="Mahal_m"),
                             function(object, newdata,...) {
                               if (!all(object@features %in% colnames(newdata))) 
                                 stop("One or more variables in the model do not exist in the data!")
                               
                               newdata <- newdata[,object@features] %>% 
                                 as.data.frame()
                               
                               out <- newdata %>% 
                                 apply( 
                                   1, 
                                   FUN=function(z) {
                                     mahal_dist <- object@presence %>% 
                                       mahalanobis(z, object@cov, inverted = T) %>%
                                       min() %>%  
                                       pchisq(df=ncol(object@presence))
                                   }
                                 )
                               out <- 1 - out
                             }
)

# Passo 4
methodInfo <- list(name=c("Mahal_m"),
                   packages=NULL,
                   modelTypes = c("po"),
                   fitParams = list(formula="standard.formula",data="sdmDataFrame"),
                   fitSettings = NULL,
                   fitFunction = ".mahal_m_fit",
                   settingRules = NULL,
                   tuneParams = NULL,
                   predictParams=list(object="model",newdata="sdmDataFrame"),
                   predictSettings=NULL,
                   predictFunction="predict"
)

# Passo 5
environment(classe_predicao) <- environment(sdm)
environment(metodo_predicao) <- environment(sdm)

# Passo 6
if (!methodInfo$name %in% names(getmethodNames())) {
  add(methodInfo, w="sdm")
}
#---------------------------------------------------------------------------------------------
# Inclui o algoritmo de dist de mahalanobis para um centro ótimo -----------------------------
classe_predicao <- setClass("Mahal", 
                            representation(
                              features="character",
                              cov="matrix",
                              presence="data.frame"
                            )
)

.mahal_fit <- function(formula, data, ...) {
  prep_data <- .prepare_data(formula, data)
  
  mahal_m <- new("Mahal")
  mahal_m@features <- colnames(prep_data$presences)
  mahal_m@presence <- prep_data$presences
  mahal_m@cov <- solve(cov(mahal_m@presence), tol=1e-20)
  
  mahal_m
}

metodo_predicao <- setMethod("predict", signature(object="Mahal"),
                             function(object, newdata,...) {
                               if (!all(object@features %in% colnames(newdata))) 
                                 stop("One or more variables in the model do not exist in the data!")
                               
                               newdata <- newdata[,object@features]
                               data_center <- colMeans(object@presence, na.rm=TRUE)
                               
                               out <- newdata %>% 
                                 mahalanobis(data_center, object@cov, inverted = T) %>%
                                 pchisq(df=ncol(object@presence))
                               
                               1 - out
                             }
)

methodInfo <- list(name=c("Mahal"),
                   packages=NULL,
                   modelTypes = c("po"),
                   fitParams = list(formula="standard.formula",data="sdmDataFrame"),
                   fitSettings = NULL,
                   fitFunction = ".mahal_fit",
                   settingRules = NULL,
                   tuneParams = NULL,
                   predictParams=list(object="model",newdata="sdmDataFrame"),
                   predictSettings=NULL,
                   predictFunction="predict"
)


environment(classe_predicao) <- environment(sdm)
environment(metodo_predicao) <- environment(sdm)

if (!methodInfo$name %in% names(getmethodNames())) {
  add(methodInfo, w="sdm")
}
#---------------------------------------------------------------------------------------------
# Inclui o algoritmo de dist euclideana ------------------------------------------------------
# https://github.com/cran/OutlierDetection/blob/master/R/distance.R
classe_predicao <- setClass("DEuc_m", 
                            representation(
                              features="character",
                              presence="data.frame"
                            )
)

.deuc_m_fit <- function(formula, data, ...) {
  prep_data <- .prepare_data(formula, data)
  
  deuc_m <- new("DEuc_m")
  deuc_m@features <- colnames(prep_data$presences)
  deuc_m@presence <- prep_data$presences %>%
    as.data.frame()

  deuc_m
}

metodo_predicao <- setMethod("predict", signature(object="DEuc_m"),
                             function(object, newdata,...) {
                               if (!all(object@features %in% colnames(newdata))) 
                                 stop("One or more variables in the model do not exist in the data!")
                               
                               newdata <- newdata[,object@features] %>% 
                                 as.data.frame()
                               
                               out <-  newdata %>% 
                                 by_row( ~ min(cdist(., object@presence))) %>% 
                                 select(ncol(newdata) + 1) %>% 
                                 unlist() %>%
                                 unname() 
                               out <- out %>% 
                                 (out %>% ecdf())
                               
                               1 - out
                             }
)

methodInfo <- list(name=c("DEuc_m"),
                   packages=c("class"),
                   modelTypes = c("po"),
                   fitParams = list(formula="standard.formula",data="sdmDataFrame"),
                   fitSettings = NULL,
                   fitFunction = ".deuc_m_fit",
                   settingRules = NULL,
                   tuneParams = NULL,
                   predictParams=list(object="model",newdata="sdmDataFrame"),
                   predictSettings=NULL,
                   predictFunction="predict"
)

environment(classe_predicao) <- environment(sdm)
environment(metodo_predicao) <- environment(sdm)

if (!methodInfo$name %in% names(getmethodNames())) {
  add(methodInfo, w="sdm")
}

#---------------------------------------------------------------------------------------------
# Inclui o algoritmo de dist euclideana  de um centro ótimo ----------------------------------
classe_predicao <- setClass("DEuc", 
                            representation(
                              features="character",
                              presence="data.frame",
                              centro="matrix",
                              modelo="ANY"
                            )
)

.deuc_fit <- function(formula, data, ...) {
  prep_data <- .prepare_data(formula, data)
  
  deuc <- new("DEuc")
  deuc@features <- colnames(prep_data$presences)
  deuc@presence <- prep_data$presences %>% 
    as.data.frame()
  
  deuc@centro <- deuc@presence %>%
    colMeans(na.rm=TRUE) %>% 
    as.data.frame() %>%
    t()
  
  deuc@modelo <- cdist(deuc@presence, deuc@centro) %>% ecdf()
    
  deuc
}

metodo_predicao <- setMethod("predict", signature(object="DEuc"),
                             function(object, newdata,...) {
                               if (!all(object@features %in% colnames(newdata))) 
                                 stop("One or more variables in the model do not exist in the data!")
                               
                               newdata <- newdata[,object@features] %>%
                                 as.data.frame()
                               
                               out <- cdist(newdata, object@centro) %>%
                                 object@modelo()
                               
                               # out <-  newdata %>% 
                               #   by_row( ~ cdist(., object@centro)) %>% 
                               #   select(ncol(newdata) + 1) %>% 
                               #   unlist() %>%
                               #   unname()
                               
                               1 - out
                             }
)

methodInfo <- list(name=c("DEuc"),
                   packages=c("class"),
                   modelTypes = c("po"),
                   fitParams = list(formula="standard.formula",data="sdmDataFrame"),
                   fitSettings = NULL,
                   fitFunction = ".deuc_fit",
                   settingRules = NULL,
                   tuneParams = NULL,
                   predictParams=list(object="model",newdata="sdmDataFrame"),
                   predictSettings=NULL,
                   predictFunction="predict"
)

environment(classe_predicao) <- environment(sdm)
environment(metodo_predicao) <- environment(sdm)

if (!methodInfo$name %in% names(getmethodNames())) {
  add(methodInfo, w="sdm")
}

#---------------------------------------------------------------------------------------------
# Inclui algoritmo Svm one-class -------------------------------------------------------------
classe_predicao <- setClass("Svm_occ", 
                            representation(
                              features="character",
                              presence='data.frame',
                              filtro='ksvm'
                            )
)

.svm_occ_fit <- function(formula, data, ...) {
  prep_data <- .prepare_data(formula, data)
  
  svm_occ <- new("Svm_occ")
  svm_occ@features <- colnames(prep_data$presences)
  svm_occ@presence <- prep_data$presences
  svm_occ@filtro <- ksvm(
      svm_occ@presence %>% as.matrix(), 
      1 %>% rep(svm_occ@presence %>% nrow()), 
      type="one-svc", 
      kernel="rbfdot", 
      kpar="automatic"
    )
  
  svm_occ
}

metodo_predicao <- setMethod("predict", signature(object="Svm_occ"),
                             function(object, newdata,...) {
                               if (!all(object@features %in% colnames(newdata))) 
                                 stop("One or more variables in the model do not exist in the data!")
                               newdata <- newdata[,object@features]
                               
                               out <- predict(object@filtro, newdata) %>%
                                 as.integer()
                               
                               out
                             }
)

methodInfo <- list(name=c("Svm_occ"),
                   packages=NULL,
                   modelTypes = c('po'),
                   fitParams = list(formula="standard.formula",data="sdmDataFrame"),
                   fitSettings = NULL, 
                   fitFunction = ".svm_occ_fit",
                   settingRules = NULL,
                   tuneParams = NULL,
                   predictParams=list(object="model",newdata="sdmDataFrame"),
                   predictSettings=NULL,
                   predictFunction="predict"
)

environment(classe_predicao) <- environment(sdm)
environment(metodo_predicao) <- environment(sdm)

if (!methodInfo$name %in% names(getmethodNames())) {
  add(methodInfo, w="sdm")
}

#----------------------------------------------------------------------------------------------
# Inclui o algoritmo Enfa ----------------------------------
# https://www.researchgate.net/publication/251284999_Ecological-Niche_Factor_Analysis_How_to_Compute_Habitat-Suitability_Maps_without_Absence_Data/stats
classe_predicao <- setClass("Enfa", 
                            representation(
                              features="character",
                              bg="data.frame",
                              pr_aus="vector",
                              dudi_pca="ANY",
                              model="ANY"
                            )
)

.enfa <- function(dudi, pr, scannf = F, nf = 1)
{
  ## Verifications
  if (!inherits(dudi, "dudi"))
    stop("object of class dudi expected")
  call <- match.call()
  if (any(is.na(dudi$tab)))
    stop("na entries in table")
  if (!is.vector(pr))
    stop("pr should be a vector")
  
  ## Bases of the function
  prb <- pr
  pr <- pr/sum(pr)
  row.w <- dudi$lw/sum(dudi$lw)
  col.w <- dudi$cw
  Z <- as.matrix(dudi$tab)
  n <- nrow(Z)
  f1 <- function(v) sum(v * row.w)
  center <- apply(Z, 2, f1)
  Z <- sweep(Z, 2, center)
  
  
  ## multiply with the square root of the column weights
  Ze <- sweep(Z, 2, sqrt(col.w), "*")
  
  ## Inertia matrices S and G
  DpZ <- apply(Ze, 2, function(x) x*pr)
  
  ## Marginality computation
  mar <- apply(Z,2,function(x) sum(x*pr))
  me <- mar*sqrt(col.w)
  Se <- crossprod(Ze, DpZ)
  Ge <- crossprod(Ze, apply(Ze,2,function(x) x*row.w))
  
  ## Computation of S^(-1/2)
  eS <- eigen(Se)
  kee <- (eS$values > 1e-9)     ## keep only the positive values
  S12 <- eS$vectors[,kee] %*% diag(eS$values[kee]^(-0.5)) %*% t(eS$vectors[,kee])
  
  ## Passage to the third problem
  W <- S12 %*% Ge %*% S12
  x <- S12%*%me
  b <- x / sqrt(sum(x^2))
  
  ## Eigenstructure of H
  H <- (diag(ncol(Ze)) - b%*%t(b)) %*% W %*% (diag(ncol(Ze)) - b%*%t(b))
  s <- eigen(H)$values[-ncol(Z)]
  
  ## Number of eigenvalues
  if (scannf) {
    barplot(s)
    cat("Select the number of specialization axes: ")
    nf <- as.integer(readLines(n = 1))
  }
  if (nf <= 0 | nf > (ncol(Ze) - 1))
    nf <- 1
  
  ## coordinates of the columns on the specialization axes
  co <- matrix(nrow = ncol(Z), ncol = nf + 1)
  tt <- data.frame((S12 %*% eigen(H)$vectors)[, 1:nf])
  ww <- apply(tt, 2, function(x) x/sqrt(col.w))
  norw <- sqrt(diag(t(as.matrix(tt))%*%as.matrix(tt)))
  co[, 2:(nf + 1)] <- sweep(ww, 2, norw, "/")
  
  ## coordinates of the columns on the marginality axis
  m <- me/sqrt(col.w)
  co[, 1] <- m/sqrt(sum(m^2))
  
  ## marginality
  m <- sum(m^2)
  
  ## Coordinates of the rows on these axes
  li <- Z %*% apply(co, 2, function(x) x*col.w)
  
  ## Output
  co <- as.data.frame(co)
  li <- as.data.frame(li)
  names(co) <- c("Mar", paste("Spe", (1:nf), sep = ""))
  row.names(co) <- dimnames(dudi$tab)[[2]]
  names(li) <- c("Mar", paste("Spe", (1:nf), sep = ""))
  enfa <- list(call = call, tab = data.frame(Z), pr = prb, cw = col.w,
               nf = nf, m = m, s = s, lw = row.w, li = li,
               co = co, mar = mar)
  class(enfa) <- "enfa"
  return(invisible(enfa))
}


.enfa_fit <- function(formula, data, ...) {
  prep_data <- .prepare_data(formula, data)
  
  enfa <- new("Enfa")
  enfa@features <- colnames(prep_data$presences)
  enfa@bg <- prep_data$absences
  enfa@pr_aus <- prep_data$presences %>% 
    bind_rows(prep_data$absences)
  
  enfa@dudi_pca <- dudi.pca(prep_data$all_data %>% select(-prep_data$nsp), scannf=F)
  enfa@model <- .enfa(enfa@dudi_pca, prep_data$all_data %>% select(prep_data$nsp) %>% pull())
  
  enfa
}

metodo_predicao <- setMethod("predict", signature(object="Enfa"),
                             function(object, newdata,...) {
                               if (!all(object@features %in% colnames(newdata))) 
                                 stop("One or more variables in the model do not exist in the data!")
                               newdata <- newdata[,object@features]
                        
                               nf <- 1 #object@model$nf
                               Zli <- object@model$li[, 1:(1+nf)]
                               f1 <- function(x) {rep(x, object@model$pr)}
                               Sli <- apply(Zli, 2, f1) #extrai os valores de Zli para os pontos de ocorrencia da sp
                               
                               m.baseline <- apply(newdata, 2, mean)
                               sd.baseline <- apply(newdata, 2, sd)
                               newdata.scale <- sweep(newdata,2, m.baseline) #diminui cada coluna de 'new.climate' da sua respectiva media em 'm.baseline'
                               newdata.scale <- as.matrix(newdata.scale) %*% diag(1/sd.baseline)
                               new_Zli <- newdata.scale %*% as.matrix(object@model$co)
                               
                               #out <- mahalanobis(Zli, colMeans(Sli), cov(Sli)) 
                               out <- new_Zli %>% mahalanobis(colMeans(Sli), cov(Sli)) %>%
                                 pchisq(df=ncol(new_Zli))
                               
                               1 - out
                             }
)

methodInfo <- list(name=c('Enfa'),
                   packages=NULL,
                   modelTypes = c('pb'),
                   fitParams = list(formula="standard.formula",data="sdmDataFrame"),
                   fitSettings = NULL,
                   fitFunction = ".enfa_fit",
                   settingRules = NULL,
                   tuneParams = NULL,
                   predictParams=list(object="model",newdata="sdmDataFrame"),
                   predictSettings=NULL,
                   predictFunction="predict"
)


environment(classe_predicao) <- environment(sdm)
environment(metodo_predicao) <- environment(sdm)


if (!methodInfo$name %in% names(getmethodNames())) {
  add(methodInfo, w="sdm")
}
#----------------------------------------------------------------------------------------------
# Inclui o algoritmo Isolation Forest --------------------------------
classe_predicao <- setClass("IsoForest", 
                            representation(
                              features="character",
                              model="ANY"
                            )
)

.iso_forest_fit <- function(formula, data, ...) {
  prep_data <- .prepare_data(formula, data)
  
  iso_forest <- new("IsoForest")
  iso_forest@features <- colnames(prep_data$presences)
  
  iso_forest@model <- isolation.forest(prep_data$presences, ndim=2, prob_pick_avg_gain=1)
  
  iso_forest
}

# Passo 3
metodo_predicao <- setMethod("predict", signature(object="IsoForest"),
                             function(object, newdata,...) {
                               if (!all(object@features %in% colnames(newdata))) 
                                 stop("One or more variables in the model do not exist in the data!")
                               newdata <- newdata[,object@features]
                               
                               out <- predict(object@model, newdata)
                                
                               out
                             }
)

# Passo 4
methodInfo <- list(name=c("IsoForest"),
                   packages=c("isotree"),
                   modelTypes = c("po"),
                   fitParams = list(formula="standard.formula",data="sdmDataFrame"),
                   fitSettings = NULL,
                   fitFunction = ".iso_forest_fit",
                   settingRules = NULL,
                   tuneParams = NULL,
                   predictParams=list(object="model",newdata="sdmDataFrame"),
                   predictSettings=NULL,
                   predictFunction="predict"
)

# Passo 5
environment(classe_predicao) <- environment(sdm)
environment(metodo_predicao) <- environment(sdm)

# Passo 6
if (!methodInfo$name %in% names(getmethodNames())) {
  add(methodInfo, w="sdm")
}

#----------------------------------------------------------------------------------------------
# Inclui o algoritmo GMM do pacote ClusterR --------------------------------
classe_predicao <- setClass("gmm_clusterr", 
                            representation(
                              features="character",
                              model="ANY"
                            )
)

.gmm_clusterr_fit <- function(formula, data, ...) {
  prep_data <- .prepare_data(formula, data)
  
  gmm_clusterr <- new("gmm_clusterr")
  gmm_clusterr@features <- colnames(prep_data$presences)
  
  gmm_clusterr@model <- GMM(prep_data$presences, 2)
  
  gmm_clusterr
}

# Passo 3
metodo_predicao <- setMethod("predict", signature(object="gmm_clusterr"),
                             function(object, newdata,...) {
                               if (!all(object@features %in% colnames(newdata))) 
                                 stop("One or more variables in the model do not exist in the data!")
                               newdata <- newdata[,object@features]
                               
                               out <- newdata %>% predict_GMM(
                                   object@model$centroids, 
                                   object@model$covariance_matrices, 
                                   object@model$weights
                                 )
                                 
                               1 - out$cluster_proba[,1]
                             }
)

# Passo 4
methodInfo <- list(name=c("gmm_clusterr"),
                   packages=c("ClusterR"),
                   modelTypes = c("po"),
                   fitParams = list(formula="standard.formula",data="sdmDataFrame"),
                   fitSettings = NULL,
                   fitFunction = ".gmm_clusterr_fit",
                   settingRules = NULL,
                   tuneParams = NULL,
                   predictParams=list(object="model",newdata="sdmDataFrame"),
                   predictSettings=NULL,
                   predictFunction="predict"
)

# Passo 5
environment(classe_predicao) <- environment(sdm)
environment(metodo_predicao) <- environment(sdm)

# Passo 6
if (!methodInfo$name %in% names(getmethodNames())) {
  add(methodInfo, w="sdm")
}
#-----------------------------------------------------------------------------
# Inclui o algoritmo GMM do pacote mclust --------------------------------
classe_predicao <- setClass("gmm_mclust", 
                            representation(
                              features="character",
                              model="ANY"
                            )
)

.gmm_mclust_fit <- function(formula, data, ...) {
  prep_data <- .prepare_data(formula, data)
  
  gmm_mclust <- new("gmm_mclust")
  gmm_mclust@features <- colnames(prep_data$presences)
  
  gmm_mclust@model <- 
    Mclust(
      prep_data$presences %>% bind_rows(prep_data$absences), 
      initialization = list(
        noise =  c(rep(T, nrow(prep_data$presences)), rep(F, nrow(prep_data$absences)))
      ), 
      G=1
    )
  
  gmm_mclust
}


# Passo 3
metodo_predicao <- setMethod("predict", signature(object="gmm_mclust"),
                             function(object, newdata,...) {
                               if (!all(object@features %in% colnames(newdata))) 
                                 stop("One or more variables in the model do not exist in the data!")
                               newdata <- newdata[,object@features]
                               
                               
                               out <- predict(object@model, newdata)
                               
                               1 - out$z[,1]
                             }
)

# Passo 4
methodInfo <- list(name=c("gmm_mclust"),
                   packages=c("mclust"),
                   modelTypes = c("pb"),
                   fitParams = list(formula="standard.formula",data="sdmDataFrame"),
                   fitSettings = NULL,
                   fitFunction = ".gmm_mclust_fit",
                   settingRules = NULL,
                   tuneParams = NULL,
                   predictParams=list(object="model",newdata="sdmDataFrame"),
                   predictSettings=NULL,
                   predictFunction="predict"
)

# Passo 5
environment(classe_predicao) <- environment(sdm)
environment(metodo_predicao) <- environment(sdm)

# Passo 6
if (!methodInfo$name %in% names(getmethodNames())) {
  add(methodInfo, w="sdm")
}