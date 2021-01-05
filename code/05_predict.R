
k = 16
load(paste("resu/DimRed_FK_kNN", k, ".Rdata", sep = ""))
load(paste("resu/DimRed_FI_kNN", k, ".Rdata", sep = ""))
load(paste("resu/cca_results_k", k, ".RData", sep = ""))

# clustering is needed to find "homogeneous" regions for the spatial CV in the prediction step
if (file.exists(paste("ap_clust_k", k, ".RData", sep = "")) && !calc_from_scratch) {
    load(paste("ap_clust_k", k, ".RData", sep = ""))
  } else {
    # get the centroids of the MTBs
    lonlat = coordinates(gCentroid(MTB, byid=TRUE))
    lonlat = as.table(lonlat)

    # select coordinates that we also used throughout the study
    idx = which(MTB$NAME %in% idx_mtb, arr.ind = TRUE)
    lonlat_new = lonlat[idx, ]

    # geographical similarities among grid cells
    DGeo  = negDistMat(lonlat_new, r = 2)

    # similarities among leading canonical variates
    DX = negDistMat(Y_fk[ , 1:4], r = 2)

    # make an assymetric distance matrix containing dim red simililarities ...
    Dmix = DX

    # ... and geographical distances to induce spatial coherence
    Dmix[lower.tri(Dmix)] = DGeo[lower.tri(DGeo)]

    # cluster the geolocations by considering "floristic" and "real" vicinity with affinity propagation
    ap_result = apcluster(Dmix, q = 0)
    ap_result

    # convert the clustering results in vector
    clust_res = array(NA, length(idx_mtb))
    for (i in 1:length(ap_result)) {
      clust_res[ap_result[[i]]] = i
    }

    # save results
    save(clust_res, file = paste("ap_clust_k", k, ".RData", sep = ""))
}


# choose colormpa from here https://kwstat.github.io/pals/
# depending on length of cluster output
if (length(ap_result) <= 12) {
  palette = sample(tol(length(ap_result)), length(ap_result))
} else if (length(ap_result) > 12 && length(ap_result) < 6*4) {
  palette = sample(stepped(length(ap_result)), length(ap_result))
} else {
  palette = sample(alphabet(length(ap_result)), length(ap_result))
}

# set color range
scale_range = c(1, length(ap_result))

# generate vector to plot
dim2plot  = as.data.frame(cbind(idx_mtb, clust_res))
colnames(dim2plot) = c("NAME", "score")

# https://stackoverflow.com/questions/40276569/reverse-order-in-r-leaflet-continuous-legend
# check plot(1, 1, pch = 19, cex = 3, col = pal(1))
pal = colorNumeric(
  palette,
  domain = scale_range,
  na.color = "transparent",
  alpha = T,
  reverse = TRUE)

pal_rev = colorNumeric(
  palette,
  domain = scale_range,
  na.color = "transparent",
  alpha = T,
  reverse = FALSE)

the_Map = map_mtb(dim2plot, idx_mtb, MTB, pal, pal_rev, scale_range)
the_Map
saveWidget(the_Map, paste("Fig-Cluster4VCfolds-", idim, "-k", k, ".html", sep = ""))


## Now we come to the prediction step
clim_dat = read_csv("data/mtbschnitt_worldclim.csv")
soil_dat = read_csv("data/mtbschnitt_soilgrids.csv")

# clean data
clim_dat = clim_dat %>% mutate_all(~ replace(., . < -99, NA))
soil_dat = soil_dat %>% mutate_all(~ replace(., . < -99, NA))

# retain only values where we have MTBs
clim_dat = filter(clim_dat, clim_dat$mtbschnitt_id %in% idx_mtb)
soil_dat = filter(soil_dat, soil_dat$mtbschnitt_id %in% idx_mtb)

# sort
clim_dat = arrange(clim_dat, mtbschnitt_id)
soil_dat = arrange(soil_dat, mtbschnitt_id)

# just checking and if ok remove the name var
if (sum(clim_dat$mtbschnitt_id == idx_mtb) == length(idx_mtb)) {
  clim_dat = dplyr::select(clim_dat, -mtbschnitt_id)
  soil_dat = dplyr::select(soil_dat, -mtbschnitt_id)
}

clim_soil <- cbind(clim_dat, soil_dat)

## reducing the number of variables for the climate and soil system to 5 each
clim_pca  = PCA(clim_dat, scale.unit = TRUE, ncp = 4, graph = FALSE)
soil_pca  = PCA(soil_dat, scale.unit = TRUE, ncp = 4, graph = FALSE)
clim_soil_pca = PCA(cbind(clim_dat, soil_dat), scale.unit = TRUE, ncp = 10, graph = FALSE)

# Figure S3
pdf("Fig-pca-var-clim.pdf")
fviz_eig(clim_pca, addlabels = TRUE, ylim = c(0, 50), barfill = "darkcyan", barcolor = "darkcyan", cex = 1.5) +
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))
dev.off()

pdf("Fig-pca-var-soil.pdf")
fviz_eig(soil_pca, addlabels = TRUE, ylim = c(0, 50), barfill = "darkcyan", barcolor = "darkcyan", cex = 1.5) +
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))
dev.off()

# just checking
# fviz_eig(clim_soil, addlabels = TRUE, ylim = c(0, 50), barfill = "darkcyan", barcolor = "darkcyan", cex = 1.5) +
#   theme(text = element_text(size = 16),
#         axis.title = element_text(size = 16),
#         axis.text = element_text(size = 12))

fviz_pca_var(clim_pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)


fviz_pca_var(soil_pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)

# Figures S4-S7 (left and right column)
for (i in 1:4) {
  for (v in c("soil", "clim")) {
    figname = paste("Fig-pca-", v, "-contrib-dim-", i, ".pdf", sep = "")
    pdf(figname)
    p = fviz_contrib(eval(parse(text = paste(v, "_pca", sep = ""))), choice = "var", axes = i, top = 50, xtickslab.rt = 90, color = "darkcyan", fill = "darkcyan")+
      theme(text = element_text(size = 16),
            axis.text.x = element_text(angle = 0),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 12)) + 
      coord_flip() +
      labs(title = paste("Contrib to Dim ", i, sep = ""), x = "Variable")
    print(p)
    dev.off()
  }
}

# interpretation
#corrplot(t(get_pca_var(clim_pca)$cos2),    is.corr = FALSE)
#corrplot(t(get_pca_var(clim_pca)$contrib), is.corr = FALSE)
#corrplot(t(get_pca_var(soil_pca)$cos2),    is.corr = FALSE)
#corrplot(t(get_pca_var(soil_pca)$contrib), is.corr = FALSE)

PC_clim = as.data.frame(get_pca_ind(clim_pca)$coord)
PC_soil = as.data.frame(get_pca_ind(soil_pca)$coord)
# check
#plot(PC_clim[, 1], PC_clim[, 2])
#fviz_pca_ind(clim_pca, col.var = "black")


# Fig 3
truncate_P = function(P) {
  P_q = quantile(P, probs = c(0.01, 0.99))
  P[P > P_q[2]] = P_q[2]
  P[P < P_q[1]] = P_q[1]
  return(P)
}

# for Florkart and Flora Incognita
for (ploter in c("clim", "soil")) {

  # looks like k = 24 is best correlated, but lower k show better embedding
  # choose colormpa from here https://kwstat.github.io/pals/
  if (ploter == "clim") {
    palette = brewer.spectral(100)
  } else {
    palette = cubehelix(100) #cubicl(100)ocean.matter(100) #brewer.rdbu(100) # brewer.rdbu(100) #
  }


  # for the leading 5 dimensions
  for (idim in 1:4) {

    # we change the (anyway arbitrary) sign of the dimensions for visual consitency
    if (ploter == "clim") {
      P = PC_clim[ , idim]
    } else {
      P = PC_soil[ , idim]
    }

    # set the min/max colors to quantiles 0.01 and 0.99
    P = truncate_P(P)

    # set the color range
    scale_range = c(min(P), max(P))

    # megr by name
    dim2plot  = as.data.frame(cbind(idx_mtb, P))
    colnames(dim2plot) = c("NAME", "score")

    # https://stackoverflow.com/questions/40276569/reverse-order-in-r-leaflet-continuous-legend
    # check plot(1, 1, pch = 19, cex = 3, col = pal(1))
    pal = colorNumeric(
        palette,
        domain = scale_range,
        na.color = "transparent",
        alpha = T,
        reverse = TRUE)

    pal_rev = colorNumeric(
        palette,
        domain = scale_range,
        na.color = "transparent",
        alpha = T,
        reverse = FALSE)

    # generate the map
    the_Map = map_mtb(dim2plot, idx_mtb, MTB, pal, pal_rev, scale_range)
    #the_Map

    # save the map
    filename = paste("Fig-", ploter, "-Dim-", idim, sep = "")
    saveWidget(the_Map, paste(filename, ".html", sep = ""))
    webshot(url = paste("file://", getwd(), "/", filename, ".html", sep = ""),
            file = paste(filename, ".png", sep = ""),
            vwidth = 1100,
            vheight = 1300,
            debug = TRUE)
  }
}

# we also consider the population density as predictors
load("resu/jacdata.RData")
length(jacdata$Population)

model_params = list()
cases2pred = c("Y_fk", "Y_fi", "CCA_fk", "CCA_fi")

for (d in cases2pred) {
  for (i in 1:5){
    model_params[[length(model_params) + 1]]  <- list(d = d, i = i)
  }
}

f = function(d, i){
  print(paste("Case: ", d, "-dim-", i, sep = ""))

  # set the prediction target
  if (d == "Y_fk") {
    Y_target  = Y_fk[ , i]
  } else if (d == "Y_fi") {
    Y_target  = Y_fi[ , i]
  } else if (d == "CCA_fk") {
    Y_target  = cca.fit$canvarx[ , i]
  } else if (d == "CCA_fi") {
    Y_target  = cca.fit$canvary[ , i]
  }

  # make one data set
  learn_predict = as.data.frame(cbind(Y_target, PC_clim, PC_soil, jacdata$Population, clust_res))
  colnames(learn_predict) = c("Y", paste("PC_clim", 1:4, sep = ""), paste("PC_soil", 1:4, sep = ""), "Population", "SpatialClust")

  # make sure we have no na values
  idx_nan       = unique(which(!is.na(rowSums(learn_predict)), arr.ind = TRUE))
  idx_learn     = sample(idx_nan, 1000)
  train_data    = learn_predict[idx_learn, ]

  # indices for spatial folds
  idx_spcv      = CreateSpacetimeFolds(train_data, spacevar = "SpatialClust", k = length(unique(clust_res)))
  idx_pred      = setdiff(colnames(train_data), c("Y", "SpatialClust"))

  # variable selection the model
  set.seed(i)
  test = list()
  model_ffs <- ffs(train_data[, idx_pred],
                   train_data[, "Y"],
                   metric = "Rsquared",
                   method = "rf",
                   tuneLength = 1,
                   verbose = FALSE,
                   trControl=trainControl(method="cv",  index = idx_spcv$index))

  # then model with the selected vars
  set.seed(i)
  model_spcv <- train(train_data[, model_ffs$selectedvars],
                      train_data[, "Y"],
                      method = "rf",
                      tuneLength = 1,
                      importance = TRUE,
                      trControl=trainControl(method="cv",  index = idx_spcv$index))
  #model_spcv
  input_par = list(d = d, i = i)
  return_list = list(model= model_spcv, input_par = input_par)
  return(return_list)
}

mod2par <- function() {
  cl   = makeForkCluster(n_nodes)
  res  = parLapplyLB(cl, seq_along(model_params), function(x) do.call(f, model_params[[x]]))
  stopCluster(cl)
  return(res)
}

models_spcv_all <- mod2par()
save(models_spcv_all, file = "models2report.RData")

# the first FI dimension
i = 18
models_spcv_all[[i]]$input_par
models_spcv_all[[i]]$model$results

#Temperature
#Wind
#Precipitation
#Radiation
#HeatStress

#SoilWaterTexture
#OrganicC
#pH
#Nutrient
#Depth

require(lattice)

for (i in 1:length(models_spcv_all)) {
  val2plot = varImp(models_spcv_all[[i]]$model)
  df = as.data.frame(cbind(rownames(val2plot$importance), val2plot$importance))
  colnames(df) = c("Predictor", "Importance")
  df = arrange(df, -Importance)
  df

  figname = paste("figs/Fig-pred-", models_spcv_all[[i]]$input_par$d, "-dim-", models_spcv_all[[i]]$input_par$i, "-varImp.pdf", sep = "")
  pdf(figname)
  p  = ggplot(data = df, aes(reorder(Predictor, Importance), Importance)) +
    geom_bar(stat="identity",  fill="sienna3") +
    xlab("Selected Predictor") +
    coord_flip() +
    theme_minimal() +
    theme(text = element_text(size = 16),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 12))
  print(p)
  dev.off()
}




# or more elegantly

pdf("figs/Fig-varSelImp-spcv.pdf")
par(mfrow = c(2,2))

Y_fk_pred_tab = array(NA, c(4*2+1, 5))
rownames(Y_fk_pred_tab) = c("PC_clim1", "PC_clim2", "PC_clim3", "PC_clim4",
                            "PC_soil1", "PC_soil2", "PC_soil3", "PC_soil4", "Population")
colnames(Y_fk_pred_tab) = paste("Y_FK", 1:5)
for (i in 1:length(models_spcv_all)) {
  if (models_spcv_all[[i]]$input_par$d == "Y_fk") {
    val2plot = varImp(models_spcv_all[[i]]$model)
    df = as.data.frame(val2plot$importance)
    Y_fk_pred_tab[colnames(t(df)), models_spcv_all[[i]]$input_par$i] = t(df)
  }
}

c = corrplot(Y_fk_pred_tab[, 1:4], is.corr = FALSE,
             na.label = "square", na.label.col = "gray", tl.col = "black")
corrlegend(vertical = TRUE, ratio.colbar = 0.4, lim.segment = "auto",
align = c("c", "l", "r"), addlabels = TRUE)



Y_fi_pred_tab = array(NA, c(4*2+1, 5))
rownames(Y_fi_pred_tab) = c("PC_clim1", "PC_clim2", "PC_clim3", "PC_clim4",
                            "PC_soil1", "PC_soil2", "PC_soil3", "PC_soil4", "Population")
colnames(Y_fi_pred_tab) = paste("Y_FI", 1:5)
for (i in 1:length(models_spcv_all)) {
  if (models_spcv_all[[i]]$input_par$d == "Y_fi") {
    val2plot = varImp(models_spcv_all[[i]]$model)
    df = as.data.frame(val2plot$importance)
    Y_fi_pred_tab[colnames(t(df)), models_spcv_all[[i]]$input_par$i] = t(df)
  }
}

corrplot(Y_fi_pred_tab[, 1:4], is.corr = FALSE, na.label = "square", na.label.col = "gray", tl.col = "black")



CCA_fk_pred_tab = array(NA, c(4*2+1, 5))
rownames(CCA_fk_pred_tab) = c("PC_clim1", "PC_clim2", "PC_clim3", "PC_clim4",
                              "PC_soil1", "PC_soil2", "PC_soil3", "PC_soil4",
                               "Population")
colnames(CCA_fk_pred_tab) = paste("Z_FK", 1:5)
for (i in 1:length(models_spcv_all)) {
  if (models_spcv_all[[i]]$input_par$d == "CCA_fk") {
    val2plot = varImp(models_spcv_all[[i]]$model)
    df = as.data.frame(val2plot$importance)
    CCA_fk_pred_tab[colnames(t(df)), models_spcv_all[[i]]$input_par$i] = t(df)
  }
}

corrplot(CCA_fk_pred_tab[, 1:4], is.corr = FALSE, na.label = "square", na.label.col = "gray", tl.col = "black")


CCA_fi_pred_tab = array(NA, c(4*2+1, 5))
rownames(CCA_fi_pred_tab) = c("PC_clim1", "PC_clim2", "PC_clim3", "PC_clim4",
                              "PC_soil1", "PC_soil2", "PC_soil3", "PC_soil4", "Population")
colnames(CCA_fi_pred_tab) = paste("Z_FI", 1:5)
for (i in 1:length(models_spcv_all)) {
  if (models_spcv_all[[i]]$input_par$d == "CCA_fi") {
    val2plot = varImp(models_spcv_all[[i]]$model)
    df = as.data.frame(val2plot$importance)
    CCA_fi_pred_tab[colnames(t(df)), models_spcv_all[[i]]$input_par$i] = t(df)
  }
}

corrplot(CCA_fi_pred_tab[, 1:4], is.corr = FALSE, na.label = "square", na.label.col = "gray", tl.col = "black")

dev.off()



## joint prediction of the can axes

model_params = list()

vec = sort(rep(1:5, 10))
vec
for (i in 1:length(vec)){
    model_params[[length(model_params) + 1]] = list(i = vec[i])
  }



f = function(i){

  Y_target_1 = cca.fit$canvarx[ , i]
  Y_target_2 = cca.fit$canvary[ , i]

  # make one data set
  learn_predict_1 = as.data.frame(cbind(Y_target_1, PC_clim, PC_soil, jacdata$Population, clust_res))
  colnames(learn_predict_1) = c("Y", paste("PC_clim", 1:4, sep = ""), paste("PC_soil", 1:4, sep = ""), "Population", "SpatialClust")

  learn_predict_2 = as.data.frame(cbind(Y_target_1, PC_clim, PC_soil, jacdata$Population, clust_res))
  colnames(learn_predict_2) = c("Y", paste("PC_clim", 1:4, sep = ""), paste("PC_soil", 1:4, sep = ""), "Population", "SpatialClust")

  learn_predict = rbind(learn_predict_1, learn_predict_2)

  # make sure we have no na values
  idx_nan       = unique(which(!is.na(rowSums(learn_predict)), arr.ind = TRUE))
  idx_learn     = sample(idx_nan, 1000)
  train_data    = learn_predict[idx_learn, ]

  # indices for spatial folds
  idx_spcv      = CreateSpacetimeFolds(train_data, spacevar = "SpatialClust", k = length(unique(clust_res)))
  idx_pred      = setdiff(colnames(train_data), c("Y", "SpatialClust"))

  # variable selection the model
  #set.seed(i)
  test = list()
  model_ffs <- ffs(train_data[, idx_pred],
                   train_data[, "Y"],
                   metric = "Rsquared",
                   method = "rf",
                   tuneLength = 1,
                   verbose = FALSE,
                   trControl=trainControl(method="cv",  index = idx_spcv$index))

  # then model with the selected vars
  #set.seed(i)
  model_spcv <- train(train_data[, model_ffs$selectedvars],
                      train_data[, "Y"],
                      method = "rf",
                      tuneLength = 1,
                      importance = TRUE,
                      trControl=trainControl(method="cv",  index = idx_spcv$index))

  input_par = list(i = i)
  return_list = list(model= model_spcv, input_par = input_par)
  return(return_list)
}

mod2par <- function() {
  cl   = makeForkCluster(40)
  res  = parLapplyLB(cl, seq_along(model_params), function(x) do.call(f, model_params[[x]]))
  stopCluster(cl)
  return(res)
}

models_spcv_all_joint <- mod2par()
save(models_spcv_all_joint, file = "models2report_joint.RData")

#Temperature
#Wind
#Precipitation
#Radiation
#HeatStress

#SoilWaterTexture
#OrganicC
#pH
#Nutrient
#Depth

pred_R2 = array(0, c(5, 2))
for (i in 1:length(models_spcv_all_joint)) {
  dim2plot = models_spcv_all_joint[[i]]$input_par$i
  dim2plot
  pred_R2[dim2plot, 1] = pred_R2[dim2plot, 1] + t(models_spcv_all_joint[[i]]$model$results["Rsquared"])
  pred_R2[dim2plot, 2] = pred_R2[dim2plot, 2] + t(models_spcv_all_joint[[i]]$model$results["RMSE"])
}

round(pred_R2/10, 2)

#0.49, 0.51, 0.38, 0.41

res_mat = array(0, c(9, 5))
rownames(res_mat) = c("PC_clim1", "PC_clim2", "PC_clim3", "PC_clim4",
                      "PC_soil1", "PC_soil2", "PC_soil3", "PC_soil4",
                      "Population")
colnames(res_mat) = c("Z1", "Z2", "Z3", "Z4", "Z5")

for (i in 1:length(models_spcv_all_joint)) {
  models_spcv_all_joint[[i]]$model
  val2plot = as.data.frame(varImp(models_spcv_all_joint[[i]]$model)$importance)
  dim2plot = models_spcv_all_joint[[i]]$input_par$i
  res_mat[rownames(val2plot), dim2plot] = res_mat[rownames(val2plot), dim2plot] + t(val2plot)
  res_mat
}

rownames(res_mat) = c("PC1(clim)", "PC2(clim)", "PC3(Clim)", "PC4(Clim)",
                      "PC1(soil)", "PC2(soil)", "PC3(soil)", "PC4(soil)",
                      "Population")

res_mat = apply(res_mat, 2, function(x) x/max(x))
res_mat = apply(t(res_mat), 1, rev)
res_mat
max(res_mat)
res_melt = melt(res_mat)
names(res_melt) = c("Pred", "Dim", "Value")



for (i in 1:4) {
  figname = paste("Fig-pboth-Z-", i, "-varImp.png", sep = "")
  figname
  png(figname)
  var2plot = paste("Z", i, sep = "")
  res_tmp = filter(res_melt, Dim == var2plot)
  p  = ggplot(data = res_tmp, aes(Pred, Value)) +
    geom_bar(stat="identity",  fill="gray", width = .99) +
    xlab("Selected Predictor") +
    ylim(0, 1) +
    coord_flip() +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 20, colour = "black"))
  print(p)
  dev.off()
}
