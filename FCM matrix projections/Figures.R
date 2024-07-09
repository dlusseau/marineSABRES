
##############################################################
#### Figures                                      ############
##############################################################

# Packages
library(tiff)
library(patchwork)
library(ggplot2)

output.folder <- './FCM matrix projections'

# Graphs
file.name <- 'graphs'
tiff.dim <- c(2800, 1800)

tiff(file.path(output.folder, paste0(file.name, ".tiff")),
                   width = tiff.dim[1], height = tiff.dim[2],
                   units = "px", pointsize = 12,  res = 300,
                   compression=c("lzw"))
p.net.gr1 + p.net.gr2 + p.net.gr3a + p.net.gr4 + 
  plot_layout(ncol = 2, guides = 'collect')  +
  plot_annotation(tag_levels = 'A', tag_suffix = '.')

dev.off()

# Boolean graphs
library(magick)

file.name <- 'boolean_networks'
tiff.dim <- c(2800, 1800)

tiff(file.path(output.folder, paste0(file.name, ".tiff")),
     width = tiff.dim[1], height = tiff.dim[2],
     units = "px", pointsize = 12,  res = 300,
     compression=c("lzw"))

image_read("./FCM matrix projections/FCM_group1.svg") |> image_ggplot() + image_read("./FCM matrix projections/FCM_group2.svg") |> image_ggplot() + 
  image_read("./FCM matrix projections/FCM_group3a.svg") |> image_ggplot() + image_read("./FCM matrix projections/FCM_group4.svg") |> image_ggplot() +
  plot_layout(ncol = 2, guides = 'collect')  +
  plot_annotation(tag_levels = 'A', tag_suffix = '.')

dev.off()



# Network progress over time and PCA results

# First round all axis labels
install.packages("scales")                  # Install scales package
library("scales")       

p.time.gr1 <- p.pca.gr2 + scale_x_continuous(labels = number_format(accuracy = 0.001))

times <- (p.time.gr1 + p.time.gr2 + p.time.gr3a + p.time.gr4 + plot_layout(tag_level = 'new',ncol = 4, axis_titles = "collect", axes = "collect", guides = 'collect')) 
PCA1 <- (p.pca.gr1[[1]] + p.pca.gr2[[1]] + p.pca.gr3a[[1]] + p.pca.gr4[[1]] + plot_layout(tag_level = 'new',ncol = 4, axis_titles = "collect", axes = "collect", guides = 'collect')) 
PCA2 <-  (p.pca.gr1[[2]] + p.pca.gr2[[2]] + p.pca.gr3a[[2]] + p.pca.gr4[[2]] + plot_layout(tag_level = 'new',ncol = 4, axis_titles = "collect", axes = "collect", guides = 'collect'))
  
patchwork <- (times / PCA1 / PCA2) +
  plot_annotation(tag_levels = c('A', '1'),
                  tag_sep = '.', tag_suffix = '.')

file.name <- 'quant'
tiff.dim <- c(3000, 3800)

tiff(file.path(output.folder, paste0(file.name, ".tiff")),
     width = tiff.dim[1], height = tiff.dim[2],
     units = "px", pointsize = 12,  res = 300,
     compression=c("lzw"))

patchwork
dev.off()


