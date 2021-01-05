
library("yacca")

for (k in c(16)) {
  try(load(paste("resu/DimRed_FK_kNN", k, ".Rdata", sep = "")))
  try(load(paste("resu/DimRed_FI_kNN", k, ".Rdata", sep = "")))

  # compute the CCA
  cca.fit <- cca(Y_fk[, 1:5], Y_fi[ , 1:5])

  save(cca.fit, file = paste("resu/cca_results_k", k, ".RData", sep = ""))
}

load(paste("resu/cca_results_k", k, ".RData", sep = ""))

#View the results
round(cca.fit$corr, 2)
summary(cca.fit)
plot(cca.fit)


par(mfrow=c(1, 2), font=1, cex=.7,  mar = c(4, 5, 2, 2))
df = as.matrix(cca.fit$xstructcorrsq[, 1:4])
graphics::barplot(df,ylim = c(0, 1.17), col = viridis(5))
df = as.matrix(cca.fit$ystructcorrsq[, 1:4])
graphics::barplot(df, ylim = c(0, 1.17), col = viridis(5),legend = paste("Isomap Dim.", 1:5), 
                  args.legend = list(x = "topright", cex = 2))

# plot for paper
dev.off()
par(mfrow=c(4, 1), font=1, cex=.7,  mar = c(4, 5, 2, 2))

for (i in 1:4) {
  plot(1:5, cca.fit$xstructcorrsq[, i],
       type = "b",
       pch = 19,
       col = "coral",
       xaxt = "none", yaxt = "none",
       xlab = "", ylab = "", lwd = 2)

  abline(h = seq(0, 1, 0.1), v = seq(1, 5, 1), lty = 3, col = "gray", lwd = 0.5)

  lines(1:5,  cca.fit$ystructcorrsq[, i],
        type = "b",
        pch = 19,
        col="cornflowerblue",
        lty = 2, lwd = 2)

  axis(1, seq(1, 5, 1))
  axis(2, seq(0, 1, 0.1), las=2)

  if (i == 1) {
    legend(4, 0.95,
           legend = c("Florkart", "", "Flora Incognita"),
           col = c("coral", "White", "cornflowerblue"),
           lty = 1:2,
           cex = 1.3,
           lwd = 2,
           box.lwd = 0,
           box.col = "white",
           bg = "white")
    mtext(side = 2, line = 3, "Fractional variance",  cex=1.3)
  }
  if (i == 4) {
    mtext(side = 1, line = 2, "Canonical variate", cex = 1.3)
  }
}



# choose colormpa from here https://kwstat.github.io/pals/
palette = cubicl(100) #brewer.rdbu(100) # brewer.rdbu(100) #

# for the different k values
for (k in c(16)) {

 load(paste("resu/cca_results_k", k, ".RData", sep = ""))

  # for Florkart and Flora Incognita
  for (ploter in c("fk", "fi")) {

    # for the leading 5 dimensions
    for (idim in 1:5) {

      # we change the (anyway arbitrary) sign of the dimensions for visual consitency
      if (ploter == "fk") {
        P = -cca.fit$canvarx[ , idim]
      } else if (ploter == "fi") {
        P = -cca.fit$canvary[ , idim]
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
      #saveWidget(the_Map, paste("Fig-", ploter, "-CanonicalVar-", idim, "-k", k, ".html", sep = ""))
    
      filename = paste(getwd(), "/figs/NewFig-", ploter, "-CanonicalVar-", idim, "-k", k, sep = "")
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
