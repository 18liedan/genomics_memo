# This is a custom R script designed to plot SMC++ results.
# As an example, I am using two populations, with the assumption that the outputs from your SMC++ are named
# smcpp_pop1.csv and smcpp_pop2.csv, but please adjust accordingly to suit your needs :)

library(ggplot2)
library(dplyr)
library(readr)
library(scales)

# 1. Load Data
pop1_data <- read_csv("smcpp_pop1.csv")
pop2_data <- read_csv("smcpp_pop2.csv") # Updated label to match your plot
all_data <- bind_rows(pop1_data, pop2_data)

# 2. Plot
ggplot(all_data, aes(x = x, y = y, group = interaction(label, population), color = population)) +
  
  # Add the log ticks (little lines) at the bottom
  annotation_logticks(sides = "b") + 
  
  # Bootstraps
  geom_step(data = subset(all_data, type == "Bootstrap"), 
            alpha = 0.2, linewidth = 0.3) +
  
  # Main Estimates
  geom_step(data = subset(all_data, type == "Main"), 
            alpha = 1.0, linewidth = 0.8) +
  
  scale_color_manual(values = c("pop1" = "firebrick", "pop2" = "dodgerblue")) +
  
  # --- X-AXIS CONFIGURATION ---
  scale_x_log10(
    # 1. Manually force the breaks to be integers
    breaks = c(100, 1000, 10000, 100000, 1000000),
    
    # 2. Format labels as 10^x
    # We use trans_format combining log10 and math_format for the cleanest look
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  
  scale_y_log10(labels = scales::comma) +
  
  # Crop to view range
  coord_cartesian(xlim = c(100, 1000, 1000000)) + 
  
  theme_bw() +
  labs(x = "Years Ago (g=9.93)", # use your species' generation time
       y = "Effective Population Size (Ne)") +
  theme(
    panel.grid.minor = element_blank()
  )
