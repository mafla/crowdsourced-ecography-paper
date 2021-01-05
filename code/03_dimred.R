# estimate standard distance matrix (do this once only as this is sloooow)
print("compute or load initial distance matrix")
if (file.exists("resu/D.Rdata") && !calc_from_scratch) {
  load("resu/D.Rdata")
} else {
  D_fk = vegdist(tibble_fk, "jaccard")
  D_fi = vegdist(tibble_fi, "jaccard")
  save(D_fk, D_fi, file = "resu/D.Rdata")
}

    # k-NNs
    kNN_vec = c(seq(6, 50, by = 2), seq(60, 100, by = 10), seq(200, 3000, by = 100))
    kNN_vec
    # number of components to keep
    n_Y = 15
    
    # function to compute Isomap embeddings for varying k values
    parallel_dimred = function(k) {
        file_name_fk <- paste("resu/DimRed_FK_kNN", k, ".Rdata", sep = "")
        file_name_fi <- paste("resu/DimRed_FI_kNN", k, ".Rdata", sep = "")
    
        if (file.exists(file_name_fk) && !calc_from_scratch){
          load(file_name_fk)
        } else {
          # florkart
          try({
            DG_fk = d2dgeo(D_fk, k)
            Y_fk  = cmdscale(DG_fk, k = n_Y)
            save(DG_fk, Y_fk, file = file_name_fk)
          })
        }
    
        if (file.exists(file_name_fi) && !calc_from_scratch) {
          load(file_name_fi)
        } else {
          # flora inkognita
          try({
            DG_fi = d2dgeo(D_fi, k)
            Y_fi  = cmdscale(DG_fi, k = n_Y)
            save(DG_fi, Y_fi, file = file_name_fi)
          })
        }
    }
    
    # initiate cluster
    cl = makeForkCluster(n_nodes)
    
    # apply function to compute Isomap embeddings for varying k values
    parLapplyLB(cl, kNN_vec, parallel_dimred)
    
    # stop the workers
    stopCluster(cl)
    
    ## Estimate explained variances
    # alocate matrix to save the explained variances
    print("estimate exlplained variances")
    res_var_fk = array(NA, dim = c(length(kNN_vec), n_Y))
    res_var_fi = array(NA, dim = c(length(kNN_vec), n_Y))
    rownames(res_var_fk) = kNN_vec
    rownames(res_var_fi) = kNN_vec
    colnames(res_var_fk) = 1:n_Y
    colnames(res_var_fi) = 1:n_Y
    
    # go through the embeddings and estimate explained variances
    if (file.exists("resu/ResidualVariances.Rdata") && !calc_from_scratch) {
      load("resu/ResidualVariances.Rdata")
    } else {
      for (k in kNN_vec) {
        # load the computed embedding
        load(paste("resu/DimRed_FK_kNN", k, ".Rdata", sep = ""))
        load(paste("resu/DimRed_FI_kNN", k, ".Rdata", sep = ""))
    
        # loop through the dimensions
        for (i in 1:n_Y) {
          # compute distance matrix of the embedding
          D_Y_fk = as.matrix(dist(Y_fk[, 1:i]))
          D_Y_fi = as.matrix(dist(Y_fi[, 1:i]))
          # estimate explained variances
          res_var_fk[toString(k), i] =  1-(cor(as.vector(DG_fk), as.vector(D_Y_fk))^2)
          res_var_fi[toString(k), i] =  1-(cor(as.vector(DG_fi), as.vector(D_Y_fi))^2)
        }
      }
      save(res_var_fk, res_var_fi, file = "resu/ResidualVariances.Rdata")
    }
    

### ----------------
### VISUALIZATION

# we have to mess around, as R has no useful plotting library
# trying to tick ggplot to do matplot
res_var_fk_df = as.data.frame(t(res_var_fk))
res_var_fi_df = as.data.frame(t(res_var_fi))
#variable for position in matrix
res_var_fk_df$dimension = 1:nrow(res_var_fk_df)
res_var_fi_df$dimension = 1:nrow(res_var_fi_df)
#reshape to long format
plot_resvar_fk = melt(res_var_fk_df, id.var = "dimension")
plot_resvar_fi = melt(res_var_fi_df, id.var = "dimension")


## Fig S1
p1 = ggplot(plot_resvar_fk, aes(x = dimension, y = value, group = variable, colour = variable)) +
  geom_line()   +
  ggtitle("Residual variance by k-NN") +
  #theme_ipsum_ps() +
  #theme_ipsum() +
  #theme_tufte() +#
  theme_minimal() +
  theme(legend.position="bottom") +
  ylim(0, 1) +
  xlim(1, 6) +
  ylab("Residual Variance") +
  xlab("Embedding Dimension")

p2 = ggplot(plot_resvar_fi, aes(x = dimension, y = value, group = variable, colour = variable)) +
  geom_line() +
  ggtitle("Residual variance by k-NN") +
  #theme_ipsum_ps() +
  #theme_ipsum() +
  #theme_tufte() +#
  theme_minimal() +
  theme(legend.position="bottom") +
  ylim(0, 1) +
  xlim(1, 6) +
  ylab("Residual Variance") +
  xlab("Embedding Dimension")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p1)

p3 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                               p2 + theme(legend.position="none"),
                               nrow=2, ncol=2),
                   mylegend, nrow=2, heights=c(30, 15))

p3 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                               p2 + theme(legend.position="none"),
                               mylegend,
                               nrow=2, ncol=2))
#grid.arrange(p1, p2 +theme(legend.position="none"),ncol=2, nrow = 2)
  
# we have learned we need some k<50 or so
# we want to know on which of the embeddings we find
# strong common patterns
kNN_vec_use = kNN_vec[kNN_vec <= 40]

# We want to understand which FK and FI dim red correlate best
res_cca = array(NA, c(length(kNN_vec_use), length(kNN_vec_use), 4))

kk = 0
for (k in kNN_vec_use) {
  kk = kk + 1
  ll = 0
  load(paste("resu/DimRed_FK_kNN", k, ".Rdata", sep = ""))
  for (l in kNN_vec_use) {
    ll = ll +1
    load(paste("resu/DimRed_FI_kNN", l, ".Rdata", sep = ""))
    res_cca[kk, ll, ] = yacca::cca(Y_fk[, 1:4], Y_fi[, 1:4])$corr[1:4]
  }
}

rownames(res_cca) = paste("FK k=", kNN_vec_use, sep = "")
colnames(res_cca) = paste("FI k=", kNN_vec_use, sep = "")

# very dirty hack to correct for this corrplot bug
# https://github.com/taiyun/corrplot/issues/122
palette = c(brewer.spectral(6), brewer.spectral(6))
corrplot(apply(res_cca, c(1, 2), function(x){sum(x[1:4])/4}),
         col =  palette,
         is.corr = FALSE)

dim(apply(res_cca, c(1, 2), function(x){sum(x[1:4])/4}))

# Fig 2

# plot for paper
dev.off()
par(cex=1.5, font=1)
y_fk = c(1, res_var_fk_df[1:10, "16"])
y_fi = c(1, res_var_fi_df[1:10, "16"])

plot(0:10, y_fk,
     type = "b",
     yaxs = "i",
     xaxs = "i",
     xlim = c(0, 10),
     ylim = c(0, 1),
     pch = 19,
     col = "coral",
     xaxt = "none", yaxt = "none",
     xlab = "", ylab = "", lwd = 2)
abline(h = seq(0, 1, 0.1), v = seq(0, 10, 1), lty = 3, col = "gray", lwd = 0.5)
axis(1, seq(0, 10, 1))
axis(2, seq(0, 1, 0.1), las=2)
lines(0:10, y_fi,
      type = "b",
      pch = 19,
      col="cornflowerblue",
      lty = 2, lwd = 2)
legend(6.5, 0.95,
       legend = c("Florkart", "", "Flora Incognita"),
       col = c("coral", "White", "cornflowerblue"),
       lty = 1:2,
       cex = 1.2,
       lwd = 2,
       box.lwd = 0,
       box.col = "white",
       bg = "white")
mtext(side = 1, line = 2, "Dimension", cex = 2)
mtext(side = 2, line = 3, "Residual variance",  cex=2)




# an alternative visualization ... more down to Earth (not used in the paper)
# some inspections
# this is why we go for rather higher k-NN values
# Lower k-NN lead + expl. variances, but
# Highe k-NN better correlations among the Isomap coordinates.
# Figures not used in the paper tough
#Corr_Y  = cor(Y_merge_all)
#corrplot(Corr_Y,
#         method = "square",
#         order = "AOE",
#         type = "upper")
#idx_x = grepl("FI", colnames(Corr_Y))
#idx_y = grepl("FK", colnames(Corr_Y))
#Corr_Y_sel_sort = abs(Corr_Y[idx_y, idx_x, drop=FALSE])
#Corr_Y_sel_sort = Corr_Y_sel_sort[order(rownames(Corr_Y_sel_sort)), order(colnames(Corr_Y_sel_sort))]
#col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
#                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
#                           "#4393C3", "#2166AC", "#053061"))
#corrplot(Corr_Y_sel_sort, method = "square", cl.lim = c(0, 0.65), col = col2(50), is.corr = FALSE)


# Fig 3

# once of the best embedding,
# once for maximal CCA
for (k in c(16)) {
  # looks like k = 24 is best correlated, but lower k show better embedding
  # choose colormpa from here https://kwstat.github.io/pals/
  palette = ocean.haline(100) #cubicl(100)ocean.matter(100) #brewer.rdbu(100) # brewer.rdbu(100) #

  load(paste("resu/DimRed_FK_kNN", k, ".Rdata", sep = ""))
  load(paste("resu/DimRed_FI_kNN", k, ".Rdata", sep = ""))

  # for Florkart and Flora Incognita
  for (ploter in c("fk", "fi")) {

    # for the leading 5 dimensions
    for (idim in 1:5) {

      # we change the (anyway arbitrary) sign of the dimensions for visual consitency
      if (ploter == "fk") {
        P = Y_fk[ , idim]
      } else if (ploter == "fi") {
        # we do this for dim-1 as we found out later that that it's these the dims that match
        if (idim == 1) {
          P = -Y_fi[ , idim]
        } else {
          if (cor(Y_fk[ , idim-1], Y_fi[ , idim]) < 0) {
            P = -Y_fi[ , idim]
          } else {
            P = Y_fi[ , idim]
          }
        }
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

      # save the map as html and png
      filename = paste(getwd(), "/figs/TEstNewFig-", ploter, "-Dim-", idim, "-k", k, sep = "")
      saveWidget(the_Map, paste(filename, ".html", sep = ""))
      webshot(url = paste("file://", filename, ".html", sep = ""),
              file = paste(filename, ".png", sep = ""),
              #vwidth = 1200,
              #vheight = 1300,
              vwidth = 500,
              vheight = 630,
              debug = TRUE)
    }
  }
}





