library(tidyverse)
library(viridis)
library(mgcv)
library(cowplot)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# ======================= #
# ---- Data processing ----
# ======================= #

# Specify simulation date + which sets have been simulated.
date <- "2023_07_19"
epsilon <- 1e-3

SPY <- seq(0.01, 0.05, 0.01)
SPA <- SPY


# Read first batch to generate empty data frame with proper variables.
dr <- read.table(str_c(date,"_Masculinizer_SPY_0.01_SPA_0.01/results_summarized.txt"), T) %>% 
  as_tibble() %>% filter(!MaleDeterminer %in% unique(MaleDeterminer))
dt <- dr
# Read data batches
for(i in SPY) for(j in SPA)
{
  filename <- str_c(date,"_Masculinizer_SPY_", i,"_SPA_",j,"/results_summarized.txt")
  dt <- read.table(filename, T) %>% 
    as_tibble()
  
  dr <- dr %>% full_join(dt, by = names(dt))
}
write.table(dr, file = "results_masculinizer_combined.txt", col.names = T, row.names = F)


# Repeat procedure above for A-as-feminizer data.
dr <- read.table(str_c(date,"_Feminizer_SPY_0.01_SPA_0.01/results_summarized.txt"), T) %>% 
  as_tibble() %>% filter(!MaleDeterminer %in% unique(MaleDeterminer))
dt <- dr

for(i in SPY) for(j in SPA)
{
  filename <- str_c(date,"_Feminizer_SPY_", i,"_SPA_",j,"/results_summarized.txt")
  dt <- read.table(filename, T) %>% 
    as_tibble()
  
  dr <- dr %>% full_join(dt, by = names(dt))
}
write.table(dr, file = "results_feminizer_combined.txt", col.names = T, row.names = F)

# ==================================== #
# ---- Data analysis - Masculinizer ----
# ==================================== #

# Read masculinizer data
dr <- read.table("results_masculinizer_combined.txt", T)

# # Plot raw data on Y frequencies
dr %>% ggplot(aes(AY, AA, color = MXYPat)) + geom_point() + facet_grid(CY ~ CA)

# Determine which SD gene is most-prevalent
dr3 <- dr %>%
  mutate(Y = ifelse(MXYPat > 0.1, T, F),
         A = ifelse(MAPat > 0.1, T, F)) %>% 
  as_tibble()

# Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
mMY <- gam(Y ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)),
                    k = 3, m = 2, bs = "ts"),
          data = dr3,
          full = T, family = "binomial")

# Check GAM for oddities
summary(mMY)
gam.check(mMY) # no oddities

# Save GAM
saveRDS(object = mMY, file = "2023_07_19_GAM_Masculinizer_SPY_SPA_RY_RA_MY.RDS")
mMY <- readRDS("2023_07_19_GAM_Masculinizer_SPY_SPA_RY_RA_RA_MY.RDS")

# Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
mMA <- gam(A ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)), 
                  k = 3, m = 2, bs = "ts"), 
           data = dr3, 
           full = T, family = "binomial")

# Check GAM for oddities
summary(mMA)
gam.check(mMA) # no oddities

# Save GAM
saveRDS(object = mMA, file = "2023_07_19_GAM_Masculinizer_SPY_SPA_RY_RA_MA.RDS")
mMA <- readRDS("2023_07_19_GAM_Masculinizer_SPY_SPA_RY_RA_MA.RDS")

# Generate data set to obtain predicted values across parameter space.
dp <- expand.grid(SPY = factor(SPY), SPA = factor(SPA),
                  RY = seq(0, 0.5, length.out = 300), RA = seq(0, 0.5, length.out = 300))

# Add predicted values from model
dpp <- dp %>% transform(Y = predict(mMY, newdata = dp, type = "response"),
                        A = predict(mMA, newdata = dp, type = "response"))

# Add nice facet strip labels (NB: forgot to rename to fSPY/fSPA instead of fCY/fCA)
dpp <- dpp %>% dplyr::mutate(fCY = str_c("italic(s)[XY]==", SPY),
                     fCA = str_c("italic(s)[A]==", SPA))

write.table(x = dpp, col.names = T, row.names = F,
            file = "2023_07_19_Masculinizer_SPY_SPA_RY_RA_predictions.txt")
dpp <- read.table("2023_07_19_Masculinizer_SPY_SPA_RY_RA_predictions.txt", T)

# ========================================= #
# ---- Data visualization - Masculinizer ----
# ========================================= #

# Make graph of raw data to aid in visual check of GAM fit - sometimes formal checks of GAM fits
# produce errors despite a nearly-flawless fit to the raw data
pm_raw_data <- dr %>% dplyr::mutate(fSPY = str_c("italic(s)[XY]==", SPY),
                                    fSPA = str_c("italic(s)[A]==", SPA)) %>% 
  ggplot(aes(RY, RA, color = MAPat)) + geom_point(size = 0.5) + 
  facet_grid(fSPA~fSPY, labeller = label_parsed) + 
  scale_x_continuous(limits = c(0, 0.5), expand = c(0,0),
                     breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
  scale_y_continuous(limits = c(0, 0.5), expand = c(0,0),
                     breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
  labs(x = bquote(italic(r)[XY]), y = bquote(italic(r)[A]), color = "Frequency") + 
  scale_color_viridis(option = "magma", begin = 0.1, end = 0.9, limits = c(0,1)) + 
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"))

# Write to output.
png(file = "2023_07_19_Masculinizer_SPY_SPA_RY_RA_raw_data.png", width = 9, height = 8, units = "in", res = 1000)
pm_raw_data
dev.off()

pdf(file = "2023_07_19_Masculinizer_SPY_SPA_RY_RA_raw_data.pdf", width = 9, height = 8)
pm_raw_data
dev.off()


# Check for iffy simulations (turns out to be a minute fraction, so not an isssue)
dpp %>% group_by(Y > 0.5, A > 0.5) %>% summarize(n()) # 165 / 2250000 where neither Y nor A > 0.5.

# Write predictions to output (only done for masculinizer data, not feminizer data).
pm <- dpp %>% filter(Y > 0.5 | A > 0.5) %>% ggplot(aes(RY, RA, fill = interaction(Y > 0.5, A > 0.5))) +
  geom_tile() + 
  facet_grid(fCA~fCY, labeller = label_parsed) + 
  scale_x_continuous(limits = c(min(dpp$RY), max(dpp$RY)), expand = c(0,0),
                     breaks = c(0,0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) + 
  scale_y_continuous(limits = c(min(dpp$RA), max(dpp$RA)), expand = c(0,0),
                     breaks = c(0,0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) + 
  scale_fill_viridis(discrete = T, option = "magma", begin = 0.9, end = 0.1,
                     labels = c(bquote(italic(Y)),
                                                bquote(italic(A)),
                                                bquote(italic(Y)*" & "*italic(A)))) + 
  labs(x = bquote(italic(r)[XY]), y = bquote(italic(r)[A]), fill = "SD gene") + 
  
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"))

png(file = "2023_07_19_Masculinizer_SPY_SPA_RY_RA.png", width = 9, height = 8, units = "in", res = 1000)
pm
dev.off()


pdf(file = "2023_07_19_Masculinizer_SPY_SPA_RY_RA.pdf", width = 9, height = 8)
pm
dev.off()

# ================================= #
# ---- Data analysis - Feminizer ----
# ==================================== #

# Read feminizer data
dr <- read.table("results_feminizer_combined.txt", T)

# Plot raw data on A frequencies
dr %>% ggplot(aes(RY, RA, color = FAMat)) + geom_point() + facet_grid(SPY ~ SPA)

dr3 <- dr %>%
  mutate(Y = ifelse(FAMat < 0.1, T, F),
         A = ifelse(FAMat > 0.9, T, F)) %>% 
  as_tibble()

# Generate GAM with specific beta estimates for EACH combination of s_XY and s_A
mMY <- gam(Y ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)),
                    k = 3, m = 2, bs = "ts"),
          data = dr3,
          full = T, family = "binomial")

# Check GAM for oddities
summary(mMY)
gam.check(mMY) # no oddities

# Save GAM
saveRDS(object = mMY, file = "2023_07_19_GAM_Feminizer_SPY_SPA_RY_RA_Y.RDS")
mMY <- readRDS("2023_07_19_GAM_Feminizer_SPY_SPA_RY_RA_Y.RDS")

# Generate GAM with specific beta estimates for EACH combination of s_XY and s_A
mMA <- gam(A ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)), 
                  k = 3, m = 2, bs = "ts"), 
           data = dr3, 
           full = T, family = "binomial")

# Check GAM for oddities
summary(mMA)
gam.check(mMA) # no oddities

# Save GAM
saveRDS(object = mMA, file = "2023_07_19_GAM_Feminizer_SPY_SPA_RY_RA_A.RDS")
mMA <- readRDS("2023_07_19_GAM_Feminizer_SPY_SPA_RY_RA_A.RDS")

# Generate data set to obtain predicted values across parameter space.
dp <- expand.grid(SPY = factor(SPY), SPA = factor(SPA),
                  RY = seq(0, 0.5, length.out = 300), RA = seq(0, 0.5, length.out = 300))

# Add predicted values from model
dpp <- dp %>% transform(Y = predict(mMY, newdata = dp, type = "response"),
                        A = predict(mMA, newdata = dp, type = "response"))

# # Add predicted values from model
# dpp <- dp %>% transform(Turnover = predict(mF, newdata = dp, type = "response"))

# Add nice facet strip labels
dpp <- dpp %>% dplyr::mutate(fSPY = str_c("italic(s)[XY]==", SPY),
                             fSPA = str_c("italic(s)[A]==", SPA))

write.table(x = dpp, col.names = T, row.names = F,
            file = "2023_07_19_Feminizer_SPY_SPA_RY_RA_predictions.txt")
dpp <- read.table("2023_07_19_Feminizer_SPY_SPA_RY_RA_predictions.txt", T)


# ========================================= #
# ---- Data visualization - Feminizer ----
# ========================================= #

dr <- read.table("results_feminizer_combined.txt", T)

pf_raw_data <- dr %>% dplyr::mutate(fSPY = str_c("italic(s)[XY]==", SPY),
                                    fSPA = str_c("italic(s)[A]==", SPA)) %>% 
  ggplot(aes(RY, RA, color = FAMat)) + geom_point(size = 0.5) + 
  facet_grid(fSPA~fSPY, labeller = label_parsed) + 
  scale_x_continuous(limits = c(0, 0.5), expand = c(0,0),
                     breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
  scale_y_continuous(limits = c(0, 0.5), expand = c(0,0),
                     breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
  labs(x = bquote(italic(r)[XY]), y = bquote(italic(r)[A]), color = "Frequency") + 
  scale_color_viridis(option = "magma", begin = 0.1, end = 0.9, limits = c(0,1)) + 
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"))


png(file = "2023_07_19_Feminizer_SPY_SPA_RY_RA_raw_data.png", width = 9, height = 8, units = "in", res = 1000)
pf_raw_data
dev.off()


pdf(file = "2023_07_19_Feminizer_SPY_SPA_RY_RA_raw_data.pdf", width = 9, height = 8)
pf_raw_data
dev.off()

# ================ #
# ---- Figure 2 ----
# ================ #

dppf <- read.table("2023_07_19_Feminizer_SPY_SPA_RY_RA_predictions.txt", T)
dppm <- read.table("2023_07_19_Masculinizer_SPY_SPA_RY_RA_predictions.txt", T)

dpf2 <- dppf %>% mutate(MaleDeterminer = F)
dpm2 <- dppm %>% mutate(MaleDeterminer = T)

dp <- dpf2 %>% full_join(dpm2)

dp2 <- dp %>% filter(SPY %in% c(0.01, 0.03, 0.05), SPA %in% c(0.01, 0.03, 0.05))


# Set up part A
pa <- dp2 %>% filter(MaleDeterminer == T) %>% 
  filter(Y > 0.5 | A > 0.5) %>% ggplot(aes(RY, RA, fill = interaction(Y > 0.5, A > 0.5))) +
  geom_tile() + 
  facet_grid(fCA~fCY, labeller = label_parsed) + 
  scale_x_continuous(limits = c(min(dpp$RY), max(dpp$RY)), expand = c(0,0),
                     breaks = c(0,0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) + 
  scale_y_continuous(limits = c(min(dpp$RA), max(dpp$RA)), expand = c(0,0),
                     breaks = c(0,0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) + 
  scale_fill_viridis(discrete = T, option = "magma", begin = 0.9, end = 0.1,
                     labels = c(bquote(italic(Y)),
                                bquote(italic(A)),
                                bquote(italic(Y)*" & "*italic(A)))) + 
  labs(x = bquote(italic(r)[XY]), y = bquote(italic(r)[A]), fill = "SD gene") + 
  
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"))


# Set up part B (note: logic rules here are inverted due to different GAM setup for feminizing A)
pb <- dp2 %>% filter(MaleDeterminer == F) %>%  filter(Y < 0.5 | A < 0.5) %>% 
  ggplot(aes(RY, RA, fill = interaction(A < 0.5, Y < 0.5))) +
  geom_tile() + 
  facet_grid(fSPA~fSPY, labeller = label_parsed) + 
  scale_x_continuous(limits = c(min(dpp$RY), max(dpp$RY)), expand = c(0,0),
                     breaks = c(0,0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) + 
  scale_y_continuous(limits = c(min(dpp$RA), max(dpp$RA)), expand = c(0,0),
                     breaks = c(0,0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) + 
  scale_fill_viridis(discrete = T, option = "magma", begin = 0.9, end = 0.1,
                     labels = c(bquote(italic(Y)),
                                bquote(italic(A)),
                                bquote(italic(Y)*" & "*italic(A)))) + 
  labs(x = bquote(italic(r)[XY]), y = bquote(italic(r)[A]), fill = "SD gene") + 
  
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"))


legend <- get_legend(pa)

p <- plot_grid(pa + theme(legend.position = "none"),
          pb + theme(legend.position = "none"),
          legend,
          labels = c("A", "B", ""),
          rel_widths = c(1,1,0.2), ncol = 3)


png(file = "2023_07_19_Figure_2.png", width = 9, height = 4.5, units = "in", res = 1000)
p
dev.off()


pdf(file = "2023_07_19_Figure_2.pdf", width = 9, height = 4.5)
p
dev.off()
