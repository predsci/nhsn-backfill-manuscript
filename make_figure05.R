rm(list=ls())
library(ggplot2)
library(gridExtra)
library(grid)

##############################################
# read the data for figure 5  and make plot 
##############################################

path = c('flu', 'cov', 'rsv')

npath = length(path)

path_title = c('(A) Influenza', '(B) COVID-19', '(C) RSV')

df_list = list()

for (ii in 1:npath) {
  filename = paste0('figure_05_',path[ii],'_data.rds')
  df_list[[path[[ii]]]] = readRDS(filename)
}

pl = list()

# Use only one plot to extract the legend (e.g., the first one)
get_legend <- function(my_plot) {
  g <- ggplotGrob(my_plot)
  legend_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  legend <- g$grobs[[legend_index]]
  return(legend)
}

glist <- list()

legend_plot_index <- 1  # just pick one for legend extraction

for (ii in 1:npath) {
  
  show_y_axis <- (ii == 1)  # only TRUE for the first plot
  
  show_legend <- (ii == legend_plot_index)
  
  # Create the plot
  p <- ggplot(df_list[[path[[ii]]]], aes(x = type, y = state, fill = frac)) +
    geom_tile() +
    geom_text(aes(label = round(frac, 2)), color = "black", size = 3) +
    scale_fill_gradient2(
      high = "#6baed6", low = "#d6936b", midpoint = 0.5, na.value = "grey50", #name = "Fraction\nOutpreform", 
      breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits=c(0,1),
      guide = guide_colourbar(barheight = unit(8, "cm"))
    ) +
    labs(
      title = path_title[ii],
      x = NULL,
      y = NULL,
      fill = "Fraction\nOutperform", 
    ) +
    # theme_minimal(base_size = 12) +
    theme_set(theme_minimal(base_size = 12) +
                theme(
                  plot.margin = margin(2, 2, 2, 2),
                  legend.box.margin = margin(0, 0, 0, 0),
                  legend.margin = margin(0, 0, 0, 0)
                )
    ) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = if (show_y_axis) element_text() else element_blank(),
      #axis.ticks.y = if (show_y_axis) element_line() else element_blank(),
      axis.title.y = if (show_y_axis) element_text() else element_blank(),
      panel.grid = element_blank()
    )
  
  # Save legend and blanked plot separately
  if (ii == legend_plot_index) {
    legend_shared <- get_legend(p)
    p <- p + theme(legend.position = "none")  # remove for layout
  } 
  glist[[ii]] <- ggplotGrob(p)
  
  
}

  # Normalize widths across all plots
  max_widths <- do.call(grid::unit.pmax, lapply(glist, function(g) g$widths))
  for (i in seq_along(glist)) {
    glist[[i]]$widths <- max_widths
  }
  
  # Combine plots and shared legend into final layout
  grid.arrange(
    arrangeGrob(grobs = glist, ncol = 3),
    legend_shared,
    ncol = 2,
    widths = c(6, 0.9)  # Adjust relative widths as needed
  )

  combined_grob <- arrangeGrob(
    arrangeGrob(grobs = glist, ncol = 3),
    legend_shared,
    ncol = 2,
    widths = c(6, 0.9)
  )
  ggsave("fig/figure_05.pdf", plot = combined_grob, width = 8.5, height = 11, units = "in")
