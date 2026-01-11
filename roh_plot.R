# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# --- 1. SETTINGS ---
FROH_THRESHOLD_MB <- 1.0  
hom_file <- "species_roh.final_roh.hom" 
sim_file <- "species_roh.all_homozygous_roh.hom"

# --- 2. LOAD & PREPARE MAIN DATA ---
df <- read.table(hom_file, header = TRUE, stringsAsFactors = FALSE)
df_sim <- read.table(sim_file, header = TRUE, stringsAsFactors = FALSE)

# Calculate Mb
df <- df %>% mutate(Mb = (POS2 - POS1) / 1e6)
df_sim <- df_sim %>% mutate(Mb = (POS2 - POS1) / 1e6)

# Create Length Classes
df$Class <- cut(df$Mb, 
                breaks = c(0.25, 0.5, 1.0, 1.5, 2.0, Inf), 
                labels = c("0.25-0.5", "0.5-1.0", "1.0-1.5", "1.5-2.0", ">2.0"),
                include.lowest = TRUE)

# --- 3. PLOT 1: N_ROH COUNTS (Boxplots per Population) ---
indiv_counts <- df %>%
  group_by(FID, IID, Class) %>%
  summarise(nROH = n(), .groups = "drop")

plot_counts <- ggplot(indiv_counts, aes(x = Class, y = nROH, fill = FID, color = FID)) +
  # Alpha 0.2 on fill makes the boxes lighter so the points are easy to see
  geom_boxplot(outlier.shape = NA, alpha = 0.2, position = position_dodge(width = 0.8)) +
  # jitterdodge ensures points stay with their specific population box
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1.2) +
  theme_bw() +
  labs(title = "Number of ROH Segments ($$N_{ROH}$$) by Population", 
       x = "Length Class (Mb)", y = "Count (n)")

ggsave("1_ROH_Counts_Boxplot.png", plot_counts, width = 11, height = 6)

# --- 4. PLOT 2: S_ROH vs N_ROH (Scatter Plot) ---
sn_sum <- df %>%
  group_by(FID, IID) %>%
  summarise(total_N = n(), total_S = sum(Mb), .groups = "drop")

plot_sn <- ggplot(sn_sum, aes(x = total_S, y = total_N, color = FID)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  labs(title = "$$N_{ROH}$$ vs $$S_{ROH}$$ per Individual", 
       x = "Total Length ($$S_{ROH}$$ in Mb)", y = "Total Number ($$N_{ROH}$$)")

ggsave("2_NROH_vs_SROH.png", plot_sn, width = 8, height = 6)

# --- 5. CALCULATE FROH & POPULATION PLOTS ---
L_cov <- df_sim %>% filter(Mb >= FROH_THRESHOLD_MB) %>% summarise(total = sum(Mb)) %>% pull(total)

pop_stats <- df %>%
  group_by(FID, IID) %>%
  summarise(
    S_ROH_1Mb = sum(Mb[Mb >= FROH_THRESHOLD_MB]),
    Avg_Len = mean(Mb),
    .groups = "drop"
  ) %>%
  mutate(Froh = replace_na(S_ROH_1Mb / L_cov, 0))

# Plot 3a: Froh by Population
plot_froh <- ggplot(pop_stats, aes(x = FID, y = Froh, fill = FID, color = FID)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste("Inbreeding Coefficient ($$F_{ROH}$$) >", FROH_THRESHOLD_MB, "Mb"), 
       x = "Population", y = "$$F_{ROH}$$")

ggsave("3a_Froh_by_Population.png", plot_froh, width = 8, height = 6)

# Plot 3b: Average ROH Length by Population
plot_avg_len <- ggplot(pop_stats, aes(x = FID, y = Avg_Len, fill = FID, color = FID)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Average ROH Segment Length by Population", 
       x = "Population", y = "Average Length (Mb)")

ggsave("3b_Avg_Length_by_Population.png", plot_avg_len, width = 8, height = 6)

# --- 6. VALIDATION PLOTS (Line Plots) ---

theme_val <- theme_bw() + theme(panel.grid.major = element_line(color = "grey80"), legend.position = "none")

# Density
val_dens <- read.table("species_roh.validation_density.txt", header = TRUE)
plot_dens <- ggplot(val_dens, aes(x = Value, y = Detected_BP, color = IID, group = IID)) +
  geom_line(alpha = 0.6) + geom_point() + theme_val +
  labs(title = "Sensitivity: Density (Flat lines = High SNP density)", x = "--homozyg-density", y = "Total BP")
ggsave("4_Val_Density.png", plot_dens, width = 8, height = 6)

# Gap (Coverage)
val_gap <- read.table("species_roh.validation_gap.txt", header = TRUE)
val_gap$Coverage <- (val_gap$Detected_BP / 1e6) / L_cov
plot_gap <- ggplot(val_gap, aes(x = Value, y = Coverage, color = IID, group = IID)) +
  geom_line(alpha = 0.6) + geom_point() + theme_val +
  labs(title = "Sensitivity: Gap (Coverage)", x = "--homozyg-gap (kb)", y = "Genome Coverage")
ggsave("5_Val_Gap.png", plot_gap, width = 8, height = 6)

# Window Het
val_het <- read.table("species_roh.validation_window_het.txt", header = TRUE)
plot_het <- ggplot(val_het, aes(x = Value, y = Detected_BP, color = IID, group = IID)) +
  geom_line(alpha = 0.6) + geom_point() + theme_val +
  labs(title = "Sensitivity: Window Heterozygosity", x = "--homozyg-window-het", y = "Total BP")
ggsave("6_Val_Het.png", plot_het, width = 8, height = 6)

# Window Missing
val_miss <- read.table("species_roh.validation_window_missing.txt", header = TRUE)
plot_miss <- ggplot(val_miss, aes(x = Value, y = Detected_BP, color = IID, group = IID)) +
  geom_line(alpha = 0.6) + geom_point() + theme_val +
  labs(title = "Sensitivity: Window Missing", x = "--homozyg-window-missing", y = "Total BP")
ggsave("7_Val_Miss.png", plot_miss, width = 8, height = 6)
