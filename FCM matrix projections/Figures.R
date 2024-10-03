##############################################################
#### Figures                                      ############
##############################################################

# Packages
library(tiff)
library(patchwork)
library(ggplot2)

output.folder <- "./FCM matrix projections"

#### Graphs
file.name <- "graphs"
tiff.dim <- c(2800, 1800)

tiff(file.path(output.folder, paste0(file.name, ".tiff")),
  width = tiff.dim[1], height = tiff.dim[2],
  units = "px", pointsize = 12, res = 300,
  compression = c("lzw")
)

p.net.gr1 + p.net.gr2 + p.net.gr3a + p.net.gr4 +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A", tag_suffix = ".")

dev.off()

#### Boolean graphs
library(magick)

gr1 <- image_read("./FCM matrix projections/FCM_group1.png")|> image_ggplot() + ggtitle('A. Group 1')
gr2 <- image_read("./FCM matrix projections/FCM_group2.png") |> image_ggplot() + ggtitle('B. Group 2')
gr3 <- image_read("./FCM matrix projections/FCM_group3a.png") |> image_ggplot() + ggtitle('C. Group 3')
gr4 <- image_read("./FCM matrix projections/FCM_group4.png") |> image_ggplot() + ggtitle('D. Group 4')


file.name <- "boolean_networks"
tiff.dim <- c(2800, 1800)

tiff(file.path(output.folder, paste0(file.name, ".tiff")),
  width = tiff.dim[1], height = tiff.dim[2],
  units = "px", pointsize = 12, res = 300,
  compression = c("lzw")
)

gr1 + gr2 + gr3 + gr4 +
  plot_layout(ncol = 2, guides = "collect") #+
  # plot_annotation(tag_levels = "A", tag_suffix = ".")

dev.off()

#### Network progress over time and PCA results

# First round all axis labels
install.packages("scales") # Install scales package
library("scales")

p.time.gr1 <- p.pca.gr2 + scale_x_continuous(labels = number_format(accuracy = 0.001))


times <- (p.time.gr1 + p.time.gr2 + p.time.gr3a + p.time.gr4 + plot_layout(tag_level = "new", ncol = 4, axis_titles = "collect", axes = "collect", guides = "collect"))

times.subs <- (p.time.gr1.subs + p.time.gr2.subs + p.time.gr3a.subs + p.time.gr4.subs + plot_layout(tag_level = "new", ncol = 4, axis_titles = "collect", axes = "collect", guides = "collect"))
PCA1 <- (p.pca.gr1[[1]] + p.pca.gr2[[1]] + p.pca.gr3a[[1]] + p.pca.gr4[[1]] + plot_layout(tag_level = "new", ncol = 4, axis_titles = "collect", axes = "collect", guides = "collect"))
PCA2 <- (p.pca.gr1[[2]] + p.pca.gr2[[2]] + p.pca.gr3a[[2]] + p.pca.gr4[[2]] + plot_layout(tag_level = "new", ncol = 4, axis_titles = "collect", axes = "collect", guides = "collect"))

patchwork <- (times / times.subs / PCA1 / PCA2) +
  plot_annotation(
    tag_levels = c("A", "1"),
    tag_sep = ".", tag_suffix = "."
  )

file.name <- "quant"
tiff.dim <- c(3000, 2500)

tiff(file.path(output.folder, paste0(file.name, ".tiff")),
  width = tiff.dim[1], height = tiff.dim[2],
  units = "px", pointsize = 12, res = 300,
  compression = c("lzw")
)

patchwork
dev.off()

#### Sensitivity initial conditions
file.name <- "ini.sensitivity"
tiff.dim <- c(1800, 2800)

tiff(file.path(output.folder, paste0(file.name, ".tiff")),
  width = tiff.dim[1], height = tiff.dim[2],
  units = "px", pointsize = 12, res = 300,
  compression = c("lzw")
)

initial.cond.sens.group1[[1]] + initial.cond.sens.group2[[1]] + initial.cond.sens.group3a[[1]] + initial.cond.sens.group4[[1]] +
  plot_layout(ncol = 2, guides = "collect", axis_titles = "collect") + # , axes = "collect"
  plot_annotation(tag_levels = "A", tag_suffix = ".")

dev.off()


#### Sensitivity matrix
file.name <- "system.sensitivity"
tiff.dim <- c(3000, 4500)

tiff(file.path(output.folder, paste0(file.name, ".tiff")),
  width = tiff.dim[1], height = tiff.dim[2],
  units = "px", pointsize = 12, res = 300,
  compression = c("lzw")
)

(resilience.sens.group1[[1]] + scale_y_continuous(breaks = scales::breaks_extended(n = 3))) + (resilience.sens.group2[[1]] + scale_y_continuous(breaks = scales::breaks_extended(n = 3))) +
  (resilience.sens.group3a[[1]] + scale_y_continuous(breaks = scales::breaks_extended(n = 3))) + (resilience.sens.group4[[1]] + scale_y_continuous(breaks = scales::breaks_extended(n = 3))) +
  plot_layout(ncol = 1, guides = "collect", axis_titles = "collect", axes = "collect") + #
  plot_annotation(tag_levels = "A", tag_suffix = ".")

dev.off()
