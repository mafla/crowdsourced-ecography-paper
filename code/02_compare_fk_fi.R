# number of grid cells
N = dim(tibble_fk)[1]

# more comment
# function to estimate at each of the N grid cells (MTBs):
# 1: Jaccard Distance between FD and FDI
# 2: Fraction of species in FK but not in the intersection of FK and FI
# 3: Fraction of species in FI but not in the intersection of FI and FK
jac_calc <- function() {
  cl <- makeForkCluster(n_nodes)
  res <- parLapplyLB(cl, 1:N,
                     function(i) {
                       # sum of species in FK
                       fk_sum = sum(tibble_fk[i, ])
                       # Intersection
                       ab_inter = sum(tibble_fk[i, ] * tibble_fi[i, ])
                       # Union
                       ab_union = sum(tibble_fk[i, ] | tibble_fi[i, ])
                       # Complements
                       a_compl  = sum(tibble_fk[i, tibble_fi[i, ] %in% 0])
                       b_compl  = sum(tibble_fi[i, tibble_fk[i, ] %in% 0])

                         list(jacdist = 1 - ab_inter / ab_union,
                            frac_fk = a_compl/ab_inter,
                            frac_fi = b_compl/ab_inter,
                            fk_sum  = fk_sum 
                            #, count_fk = sum(a_compl) - sum(ab_inter)
                            #, count_fi = sum(b_compl) - sum(ab_inter)
                            )
                     })
  stopCluster(cl)
  bind_rows(res)
}

# call function
jac_calc_res <- jac_calc()

### VISUALIZATION
# choose colormap e.g. from here https://kwstat.github.io/pals/
palette =  ocean.matter(100)

# loop through the three metrics:
for (ploter in c("jacdist", "frac_fk", "frac_fi")) {

  # we do an interactve plot that we save as html and then export from the browser
  # there is no elegant mapping tool in R ... sorry
  if (ploter == "jacdist") {
    P = jac_calc_res$jacdist
    scale_range = c(0.5, 1)
  } else if (ploter == "frac_fk") {
    P = jac_calc_res$frac_fk
    P[P > 10] = 10
    scale_range = c(0, max(P))
  } else if (ploter == "frac_fi") {
    P = jac_calc_res$frac_fi
    P[P > 10] = 10
    P[is.na(P)] = NA
    P[is.nan(P)] = NA
    P[is.infinite(P)] = NA
    scale_range = c(0, max(P, na.rm = T))
  }

  # prep data For plotting
  dim2plot  = as.data.frame(cbind(idx_mtb, P))
  colnames(dim2plot) = c("NAME", "score")

  # generate two colormap: one for plotting of for the legend
  # sorry for this hack, check here for why:
  # https://stackoverflow.com/questions/40276569/reverse-order-in-r-leaflet-continuous-legend
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

  # generate map
  the_Map = map_mtb(dim2plot, idx_mtb, MTB, pal, pal_rev, scale_range)
  the_Map

  # save the map as html and png
  filename = paste(getwd(), "/figs/Fig-", ploter, sep = "")
  saveWidget(the_Map, paste(filename, ".html", sep = ""))
  webshot(url = paste("file://", filename, ".html", sep = ""),
          file = paste(filename, ".png", sep = ""),
          vwidth = 1100,
          vheight = 1300,
          debug = TRUE)
}


# Function to plot color bar extra
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1, cex.axis = 3)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

# generate the legends extra - this is terribly stupid and we appologize for this,
# but we need to save the legends extra and join then ina graphics prog for the final plots
# reason is that the leaflet color legends cannot be scaled automatically
par(oma = c(0, 4, 1, 1))
color.bar(rev(ocean.matter(100)), 0.5, 1)
dev.off()

#the_Map
### Interpret the map
jacdata = data.frame(NAME = idx_mtb,
                      JaccardDistance = jac_calc_res$jacdist,
                      FlorkartExceedance = jac_calc_res$frac_fk,
                      FloraIncExceedance = jac_calc_res$frac_fi,
                      FlorkartSum = jac_calc_res$fk_sum)

# rename the zensus data to merge
names(zensus)[1] = "NAME"

# merge
jacdata = merge(zensus, jacdata)

# turn it into English :-)
jacdata = rename(jacdata, Population = Bevoelkerungszahl)

# we have selected a priori which geo locations to show  using
# plot(log10(x), log10(y), bty="l")
# identify(log10(x), log10(y),1)
# based on these we got the indices
idx2name = c(29, 57, 102, 107, 128, 206, 777, 886, 961, 1221, 1273, 1300, 1426, 1482, 1677, 1698, 1701, 1732, 1753, 1831, 2065, 2210, 2216, 2254, 2478, 2580, 2826, 2893, 2927, 2991, 2998)

# extract these names from the shape files
name2print = MTB@data[as.character(MTB@data$NAME) %in% sprintf("%04d", merge(zensus, jacdata)$NAME[idx2name]),"LONGNAME"]
name2print

# replace German umlauts
name2print = name2print %>% str_replace_all(c("Kniepsand (Insel Amrum)" = "Amrum_Island",
     "D\303\274sseldorf" = "Duesseldorf",
     "D\303\274ren" = "Dueren",
     "Ein\303\266dsbach" = "Einoedsbach"))

# do the plot (needs to be exported manually for nice margins I guess...)
pdf("figs/Fig1-def.pdf")
par(mfrow=c(1, 3), mar=c(5,6,4,1)+.1)
x = log10(jacdata$Population)
y = log10(jacdata$JaccardDistance)
plot(x, y, 'n',
     xlab = expression('log'[10]*'(Population count)'),
     ylab = expression('log'[10]*'(Jaccard Distance)'),
     cex.axis = 1.5,
     cex.lab = 2,
     xlim = c(-0.5, 6.5))
points(x[setdiff(1:N, idx2name)],
       y[setdiff(1:N, idx2name)],
       col="gray", pch = 19)
textplot(x[idx2name], y[idx2name], name2print, cex=1.1, new = FALSE)

y = log10(jacdata$FlorkartExceedance)
plot(x, y, 'n',
     xlab = expression('log'[10]*'(Population count)'),
     ylab = expression('log'[10]*'(Florkart exceedance)'),
     cex.axis = 1.5,
     cex.lab = 2,
     xlim = c(-0.5, 6.5))
points(x[setdiff(1:N, idx2name)],
       y[setdiff(1:N, idx2name)],
       col="gray", pch=19)

y = log10(jacdata$FloraIncExceedance)
plot(x, y, 'n',
     xlab = expression('log'[10]*'(Population count)'),
     ylab = expression('log'[10]*'(Flora Incognita exceedance)'),
     cex.axis = 1.5,
     cex.lab = 2,
     xlim = c(-0.5, 6.5))
points(x[setdiff(1:N, idx2name)],
       y[setdiff(1:N, idx2name)],
       col="gray", pch=19)

dev.off()
