library(patchwork)
patch <- (p11 / legend_common / p21 + plot_layout(heights = c(5, 1, 5))) |
  (p12 / plot_spacer() / p22 + plot_layout(heights = c(5, 1, 5))) |
  (p13 / plot_spacer() / p23 + plot_layout(heights = c(5, 1, 5))) 

patch + plot_annotation(tag_levels = 'a')

  
