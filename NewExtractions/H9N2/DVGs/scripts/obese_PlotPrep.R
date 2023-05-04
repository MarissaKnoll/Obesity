PlotTheme1 = theme_bw() +
            theme(axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              strip.background = element_rect(colour="black", fill="white"),
              text = element_text(size = 14))

PlotTheme2 = theme_bw() +
                theme(legend.key = element_blank(),
                strip.background = element_rect(colour="black", fill="white"),
              text = element_text(size = 14))

# completely blank backgrounds
PlotTheme3 = theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                legend.key = element_blank(),
                strip.background = element_rect(colour="black", fill="white"),
              text = element_text(size = 14))

DAYS = c('stock','d02','d04','d06','d08','d10','d12')
weight_colors = c('black','#66CCff','#ff9933')
WEIGHTS = c('stock','lean','obese')
names(weight_colors) = WEIGHTS
weight_colFill <- scale_fill_manual(name = "weight", values = weight_colors)
weight_colScale <- scale_colour_manual(name = "weight", values = weight_colors)

diet = c("Obese","Lean","Control")
dietColors = c("#FF9933","#66CCFF","#606060")
names(dietColors) = diet
DietcolScale_fill <- scale_fill_manual(name = "grp",values = dietColors)
DietcolScale <- scale_colour_manual(name = "grp",values = dietColors)