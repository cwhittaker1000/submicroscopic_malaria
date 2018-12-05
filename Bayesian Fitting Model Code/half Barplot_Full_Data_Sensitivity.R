library(binom)

# Plotting barplot with equal sized groupings of prevalences
vector_PCR_prev_ordering <- order(full_data$PCR_Prev)
full_data_bar_plotting <- full_data[vector_PCR_prev_ordering, ]
full_data_bar_plotting$PCR_Prev[seq(1, 289, 32)]
breaks <- c(0, 0.032346869, 0.072202166, 0.123348018, 
            0.207373272, 0.344827586, 0.435483871, 0.558974359, 0.686431015, 0.967391304)

full_data_bar_plotting$Groups <- cut(x = full_data_bar_plotting$PCR_Prev, breaks = breaks)
full_data_bar_plotting$Sensitivity <- full_data_bar_plotting$Micro_Prev/full_data_bar_plotting$PCR_Prev

Sens_by_group <- tapply(full_data_bar_plotting$Sensitivity, full_data_bar_plotting$Groups, mean)
names_for_subsetting <- names(Sens_by_group)
names(Sens_by_group) = ""

lower_ci <- c()
upper_ci <- c()

for (i in 1:length(names_for_subsetting)) {
  
  sd <- sd(full_data_bar_plotting$Sensitivity[full_data_bar_plotting$Groups == names_for_subsetting[i]])
  stderr <- sd/sqrt(length(full_data_bar_plotting$Sensitivity[full_data_bar_plotting$Groups == names_for_subsetting[i]]))
  
  low_ci <- Sens_by_group[i] - 1.96 * stderr
  high_ci <- Sens_by_group[i] + 1.96 * stderr 
  
  lower_ci <- c(lower_ci, low_ci)
  upper_ci <- c(upper_ci, high_ci)
  
}


# set up empty plot with sensible x and y lims
plot(0, 0, ylim=c(0, 1), xlim = c(0, 1), xlab = "PCR Prevalence (%)", ylab = "Microscopy Sensitivity")

# draw data of data frame 1 and 2
for (i in 1:length(Sens_by_group)){
  
  rect(xleft = breaks[i], ybottom = 0, xright = breaks[i + 1], ytop = Sens_by_group[i])
  
  arrows(x0 = (breaks[i] + breaks[i + 1])/2, y0 = lower_ci[i], 
         x1 = (breaks[i] + breaks[i + 1])/2, y1 = upper_ci[i], 
         col=1, angle=90, code=3, length = 0.05)
  

}

