library(tidyverse)
library(viridis)
library(mgcv)
library(cowplot)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# Specify simulation date + which sets have been simulated.
date <- "2023_08_18"
epsilon <- 1e-3

SPY <- seq(0.01, 0.05, 0.02)
SPA <- SPY

# ======================= #
# ---- Data processing ----
# ======================= #


# PS Data
# Read first batch
dr <- read.table(str_c(date,"_Mas_PS_SPY_0.01_SPA_0.01/results_summarized.txt"), T) %>% 
  as_tibble() %>% filter(!MaleDeterminer %in% unique(MaleDeterminer))
dt <- dr
# Read remaining batches
for(i in SPY) for(j in SPA)
{
  filename <- str_c(date,"_Mas_PS_SPY_", i,"_SPA_",j,"/results_summarized.txt")
  dt <- read.table(filename, T) %>% 
    as_tibble()
  
  dr <- dr %>% full_join(dt, by = names(dt))
}
write.table(dr, file = "results_Mas_PS_combined.txt", col.names = T, row.names = F)


# Read first batch
dr <- read.table(str_c(date,"_Fem_PS_SPY_0.01_SPA_0.01/results_summarized.txt"), T) %>% 
  as_tibble() %>% filter(!MaleDeterminer %in% unique(MaleDeterminer))
dt <- dr
# Read remaining batches
for(i in SPY) for(j in SPA)
{
  filename <- str_c(date,"_Fem_PS_SPY_", i,"_SPA_",j,"/results_summarized.txt")
  dt <- read.table(filename, T) %>% 
    as_tibble()
  
  dr <- dr %>% full_join(dt, by = names(dt))
}
write.table(dr, file = "results_Fem_PS_combined.txt", col.names = T, row.names = F)



# SP data
# Read first batch
dr <- read.table(str_c(date,"_Mas_SP_SPY_0.01_SPA_0.01/results_summarized.txt"), T) %>% 
  as_tibble() %>% filter(!MaleDeterminer %in% unique(MaleDeterminer))
dt <- dr
# Read remaining batches
for(i in SPY) for(j in SPA)
{
  filename <- str_c(date,"_Mas_SP_SPY_", i,"_SPA_",j,"/results_summarized.txt")
  dt <- read.table(filename, T) %>% 
    as_tibble()
  
  dr <- dr %>% full_join(dt, by = names(dt))
}
write.table(dr, file = "results_Mas_SP_combined.txt", col.names = T, row.names = F)


# Read first batch
dr <- read.table(str_c(date,"_Fem_SP_SPY_0.01_SPA_0.01/results_summarized.txt"), T) %>% 
  as_tibble() %>% filter(!MaleDeterminer %in% unique(MaleDeterminer))
dt <- dr
# Read remaining batches
for(i in SPY) for(j in SPA)
{
  filename <- str_c(date,"_Fem_SP_SPY_", i,"_SPA_",j,"/results_summarized.txt")
  dt <- read.table(filename, T) %>% 
    as_tibble()
  
  dr <- dr %>% full_join(dt, by = names(dt))
}
write.table(dr, file = "results_Fem_SP_combined.txt", col.names = T, row.names = F)





# ========================================= #
# ---- Data analysis - Masculinizer - PS ----
# ========================================= #

# Read masculinizer data
dr <- read.table("results_Mas_PS_combined.txt", T)

# # Plot raw data on Y frequencies
# dr %>% ggplot(aes(AY, AA, color = MXYPat)) + geom_point() + facet_grid(CY ~ CA)

# Filter out runs with mixed systems (defined as both A and Y have frequency over >0.001)
dr3 <- dr %>%
  dplyr::mutate(Y = ifelse(MXYPat > 0.1, T, F),
                A = ifelse(MAPat > 0.1, T, F)) %>% 
  as_tibble()


dr3 %>% ggplot(aes(RY, RA, color = interaction(Y, A))) + geom_point() + facet_grid(SPA ~ SPY)
# Plot raw data on turnover rates
# dr3 %>% ggplot(aes(RY, RA, color = Turnover)) + geom_point() + facet_grid(SPY ~ SPA)

# Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
mMY <- gam(Y ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)),
                    k = 3, m = 2, bs = "ts"),
          data = dr3,
          full = T, family = "binomial")

# Check GAM for oddities
summary(mMY)
gam.check(mMY) # no oddities

# Save GAM
saveRDS(object = mMY, file = "2023_08_18_GAM_Mas_PS_SPY_SPA_RY_RA_RA_MY.RDS")
mMY <- readRDS("2023_08_18_GAM_Mas_PS_SPY_SPA_RY_RA_RA_MY.RDS")

# Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
mMA <- gam(A ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)),
                  k = 3, m = 2, bs = "ts"),
           data = dr3,
           full = T, family = "binomial")

# Check GAM for oddities
summary(mMA)
gam.check(mMA) # no oddities

# Save GAM
saveRDS(object = mMA, file = "2023_08_18_GAM_Mas_PS_SPY_SPA_RY_RA_MA.RDS")
mMA <- readRDS("2023_08_18_GAM_Mas_PS_SPY_SPA_RY_RA_MA.RDS")



# Generate data set to obtain predicted values across parameter space.
dp <- expand.grid(SPY = factor(SPY), SPA = factor(SPA),
                  RY = seq(0, 0.5, length.out = 300), RA = seq(0, 0.5, length.out = 300))

# Add predicted values from model
dpp <- dp %>% transform(Y = predict(mMY, newdata = dp, type = "response"),
                        A = predict(mMA, newdata = dp, type = "response"))

# Add nice facet strip labels
dpp <- dpp %>% dplyr::mutate(fCY = str_c("italic(u)[1]==", SPY),
                             fCA = str_c("italic(u)[2]==", SPA))

write.table(x = dpp, col.names = T, row.names = F,
            file = "2023_11_14_Mas_PS_SPY_SPA_RY_RA_predictions.txt")

# ====================================== #
# ---- Data analysis - Feminizer - PS ----
# ====================================== #

# Read feminizer data
dr <- read.table("results_Fem_PS_combined.txt", T)

# Plot raw data on Y frequencies
dr %>% ggplot(aes(RY, RA, color = FAMat)) + geom_point() + facet_grid(SPY ~ SPA)

dr3 <- dr %>%
  dplyr::mutate(Y = ifelse(FAMat < 0.1, T, F),
                A = ifelse(FAMat > 0.9, T, F)) %>% 
  as_tibble()

dr3 %>% ggplot(aes(RY, RA, color = interaction(Y,A))) + geom_point() + facet_grid(SPY ~ SPA)

# Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
mMY <- gam(Y ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)),
                  k = 3, m = 2, bs = "ts"),
           data = dr3,
           full = T, family = "binomial")

# Check GAM for oddities
summary(mMY)
gam.check(mMY) # no oddities

# Save GAM
saveRDS(object = mMY, file = "2023_11_14_GAM_Fem_PS_SPY_SPA_RY_RA_Y.RDS")
mMY <- readRDS("2023_11_14_GAM_Fem_PS_SPY_SPA_RY_RA_Y.RDS")

# Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
mMA <- gam(A ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)), 
                  k = 3, m = 2, bs = "ts"), 
           data = dr3, 
           full = T, family = "binomial")

# Check GAM for oddities
summary(mMA)
gam.check(mMA) # no oddities

# Save GAM
saveRDS(object = mMA, file = "2023_11_14_GAM_Fem_PS_SPY_SPA_RY_RA_A.RDS")
mMA <- readRDS("2023_11_14_GAM_Fem_PS_SPY_SPA_RY_RA_A.RDS")

# Generate data set to obtain predicted values across parameter space.
dp <- expand.grid(SPY = factor(SPY), SPA = factor(SPA),
                  RY = seq(0, 0.5, length.out = 300), RA = seq(0, 0.5, length.out = 300))

# Add predicted values from model
dpp <- dp %>% transform(Y = predict(mMY, newdata = dp, type = "response"),
                        A = predict(mMA, newdata = dp, type = "response"))

# Add nice facet strip labels
dpp <- dpp %>% dplyr::mutate(fSPY = str_c("italic(u)[1]==", SPY),
                             fSPA = str_c("italic(u)[2]==", SPA))

write.table(x = dpp, col.names = T, row.names = F,
            file = "2023_11_14_Fem_PS_SPY_SPA_RY_RA_predictions.txt")

# ========================================= #
# ---- Data analysis - Masculinizer - SP ----
# ========================================= #

# Read masculinizer data
dr <- read.table("results_Mas_SP_combined.txt", T)

# # Plot raw data on Y frequencies
# dr %>% ggplot(aes(AY, AA, color = MXYPat)) + geom_point() + facet_grid(CY ~ CA)

# Filter out runs with mixed systems (defined as both A and Y have frequency over >0.001)
dr3 <- dr %>%
  dplyr::mutate(Y = ifelse(MXYPat > 0.1, T, F),
         A = ifelse(MAPat > 0.1, T, F)) %>% 
  as_tibble()


dr3 %>% ggplot(aes(RY, RA, color = interaction(Y, A))) + geom_point() + facet_grid(SPA ~ SPY)
# Plot raw data on turnover rates
# dr3 %>% ggplot(aes(RY, RA, color = Turnover)) + geom_point() + facet_grid(SPY ~ SPA)

# # Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
# mMY <- gam(Y ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)),
#                     k = 3, m = 2, bs = "ts"),
#           data = dr3,
#           full = T, family = "binomial")
# 
# # Check GAM for oddities
# summary(mMY)
# gam.check(mMY) # no oddities
# 
# # Save GAM
# saveRDS(object = mMY, file = "2023_08_18_GAM_Mas_SP_SPY_SPA_RY_RA_RA_MY.RDS")
mMY <- readRDS("2023_08_18_GAM_Mas_SP_SPY_SPA_RY_RA_RA_MY.RDS")



# # Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
# mMYk4 <- gam(Y ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)),
#                   k = 4, m = 2, bs = "ts"), 
#            data = dr3, 
#            full = T, family = "binomial")
# 
# # Check GAM for oddities
# summary(mMYk4)
# gam.check(mMYk4) # no oddities
# 
# # Save GAM
# saveRDS(object = mMY, file = "2023_08_18_GAM_Mas_SP_SPY_SPA_RY_RA_MYk4.RDS")
# mMY <- readRDS("2023_08_18_GAM_Mas_SP_SPY_SPA_RY_RA_MYk4.RDS")
# 
# # Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
# mMA <- gam(A ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)), 
#                   k = 3, m = 2, bs = "ts"), 
#            data = dr3, 
#            full = T, family = "binomial")
# 
# # Check GAM for oddities
# summary(mMA)
# gam.check(mMA) # no oddities
# 
# # Save GAM
# saveRDS(object = mMA, file = "2023_08_18_GAM_Mas_SP_SPY_SPA_RY_RA_MA.RDS")
mMA <- readRDS("2023_08_18_GAM_Mas_SP_SPY_SPA_RY_RA_MA.RDS")



# Generate data set to obtain predicted values across parameter space.
dp <- expand.grid(SPY = factor(SPY), SPA = factor(SPA),
                  RY = seq(0, 0.5, length.out = 300), RA = seq(0, 0.5, length.out = 300))

# Add predicted values from model
dpp <- dp %>% transform(Y = predict(mMY, newdata = dp, type = "response"),
                        A = predict(mMA, newdata = dp, type = "response"))

# Add nice facet strip labels
dpp <- dpp %>% dplyr::mutate(fCY = str_c("italic(u)[1]==", SPY),
                     fCA = str_c("italic(u)[2]==", SPA))

write.table(x = dpp, col.names = T, row.names = F,
            file = "2023_11_14_Mas_SP_SPY_SPA_RY_RA_predictions.txt")
dpp <- read.table("2023_11_14_Mas_SP_SPY_SPA_RY_RA_predictions.txt", T)
summary(dpp)

# ====================================== #
# ---- Data analysis - Feminizer - SP ----
# ====================================== #

# Read feminizer data
dr <- read.table("results_Fem_SP_combined.txt", T)

# Plot raw data on Y frequencies
dr %>% ggplot(aes(RY, RA, color = FAMat)) + geom_point() + facet_grid(SPY ~ SPA)

dr3 <- dr %>%
  dplyr::mutate(Y = ifelse(FAMat < 0.1, T, F),
         A = ifelse(FAMat > 0.9, T, F)) %>% 
  as_tibble()

dr3 %>% ggplot(aes(RY, RA, color = interaction(Y,A))) + geom_point() + facet_grid(SPY ~ SPA)

# Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
mMY <- gam(Y ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)),
                    k = 3, m = 2, bs = "ts"),
          data = dr3,
          full = T, family = "binomial")

# Check GAM for oddities
summary(mMY)
gam.check(mMY) # no oddities

# Save GAM
saveRDS(object = mMY, file = "2023_11_14_GAM_Fem_SP_SPY_SPA_RY_RA_Y.RDS")
mMY <- readRDS("2023_11_14_GAM_Fem_SP_SPY_SPA_RY_RA_Y.RDS")

# Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
mMA <- gam(A ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)), 
                  k = 3, m = 2, bs = "ts"), 
           data = dr3, 
           full = T, family = "binomial")

# Check GAM for oddities
summary(mMA)
gam.check(mMA) # no oddities

# Save GAM
saveRDS(object = mMA, file = "2023_11_14_GAM_Fem_SP_SPY_SPA_RY_RA_A.RDS")
mMA <- readRDS("2023_11_14_GAM_Fem_SP_SPY_SPA_RY_RA_A.RDS")

# Generate data set to obtain predicted values across parameter space.
dp <- expand.grid(SPY = factor(SPY), SPA = factor(SPA),
                  RY = seq(0, 0.5, length.out = 300), RA = seq(0, 0.5, length.out = 300))

# Add predicted values from model
dpp <- dp %>% transform(Y = predict(mMY, newdata = dp, type = "response"),
                        A = predict(mMA, newdata = dp, type = "response"))

# Add nice facet strip labels
dpp <- dpp %>% dplyr::mutate(fSPY = str_c("italic(u)[1]==", SPY),
                             fSPA = str_c("italic(u)[2]==", SPA))

write.table(x = dpp, col.names = T, row.names = F,
            file = "2023_11_14_Fem_SP_SPY_SPA_RY_RA_predictions.txt")
dpp <- read.table("2023_11_14_Fem_SP_SPY_SPA_RY_RA_predictions.txt", T)
summary(dpp)


# ================ #
# ---- Figure 3 ----
# ================ #

dppf_PS <- read.table("2023_11_14_Fem_PS_SPY_SPA_RY_RA_predictions.txt", T)
dppm_PS <- read.table("2023_11_14_Mas_PS_SPY_SPA_RY_RA_predictions.txt", T)
dppf_SP <- read.table("2023_11_14_Fem_SP_SPY_SPA_RY_RA_predictions.txt", T)
dppm_SP <- read.table("2023_11_14_Mas_SP_SPY_SPA_RY_RA_predictions.txt", T)

dpf_PS2 <- dppf_PS %>% dplyr::mutate(MaleDeterminer = F, Type = "PS")
dpm_PS2 <- dppm_PS %>% dplyr::mutate(MaleDeterminer = T, Type = "PS")
dpf_SP2 <- dppf_SP %>% dplyr::mutate(MaleDeterminer = F, Type = "SP")
dpm_SP2 <- dppm_SP %>% dplyr::mutate(MaleDeterminer = T, Type = "SP")

dp <- dpf_PS2 %>% full_join(dpm_PS2) %>% full_join(dpf_SP2) %>% full_join(dpm_SP2)

dp2 <- dp %>% filter(SPY %in% c(0.01, 0.03, 0.05), SPA %in% c(0.01, 0.03, 0.05))

# Set up part A
pa <- dp2 %>% filter(MaleDeterminer == T, Type == "PS") %>% 
  filter(Y > 0.5 | A > 0.5) %>% ggplot(aes(RY, RA, fill = interaction(Y > 0.5, A > 0.5))) +
  geom_tile() + 
  facet_grid(fCA~fCY, labeller = label_parsed) + 
  scale_x_continuous(limits = c(min(dp2$RY), max(dp2$RY)), expand = c(0,0),
                     breaks = c(0,0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) + 
  scale_y_continuous(limits = c(min(dp2$RA), max(dp2$RA)), expand = c(0,0),
                     breaks = c(0,0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) + 
  scale_fill_viridis(discrete = T, option = "magma", begin = 0.9, end = 0.1,
                     labels = c(bquote(italic(s)[12]),
                                bquote(italic(s)[22]),
                                bquote(italic(s)[12]*" & "*italic(s)[22]))) + 
  labs(x = bquote(italic(r)[1]), y = bquote(italic(r)[2]), fill = "SD gene") + 
  
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 10))


pb <- dp2 %>% filter(MaleDeterminer == T, Type == "SP") %>% 
  filter(Y > 0.5 | A > 0.5) %>% ggplot(aes(RY, RA, fill = interaction(Y > 0.5, A > 0.5))) +
  geom_tile() + 
  facet_grid(fCA~fCY, labeller = label_parsed) + 
  scale_x_continuous(limits = c(min(dp2$RY), max(dp2$RY)), expand = c(0,0),
                     breaks = c(0,0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) + 
  scale_y_continuous(limits = c(min(dp2$RA), max(dp2$RA)), expand = c(0,0),
                     breaks = c(0,0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) + 
  scale_fill_viridis(discrete = T, option = "magma", begin = 0.9, end = 0.1,
                     labels = c(bquote(italic(s)[12]),
                                bquote(italic(s)[22]),
                                bquote(italic(s)[12]*" & "*italic(s)[22]))) + 
  labs(x = bquote(italic(r)[1]), y = bquote(italic(r)[2]), fill = "SD gene") + 
  
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 10))

# Set up part B (note: logic rules here are inverted due to different GAM setup for feminizing A)
pc <- dp2 %>% filter(MaleDeterminer == F, Type == "PS") %>%  filter(Y < 0.5 | A < 0.5) %>% 
  ggplot(aes(RY, RA, fill = interaction(A < 0.5, Y < 0.5))) +
  geom_tile() + 
  facet_grid(fSPA~fSPY, labeller = label_parsed) + 
  scale_x_continuous(limits = c(min(dp2$RY), max(dp2$RY)), expand = c(0,0),
                     breaks = c(0,0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) + 
  scale_y_continuous(limits = c(min(dp2$RA), max(dp2$RA)), expand = c(0,0),
                     breaks = c(0,0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) + 
  scale_fill_viridis(discrete = T, option = "magma", begin = 0.9, end = 0.1,
                     labels = c(bquote(italic(s)[12]),
                                bquote(italic(s)[22]),
                                bquote(italic(s)[12]*" & "*italic(s)[22]))) + 
  labs(x = bquote(italic(r)[1]), y = bquote(italic(r)[2]), fill = "SD gene") + 
  
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 10))

pd <- dp2 %>% filter(MaleDeterminer == F, Type == "SP") %>%  filter(Y < 0.5 | A < 0.5) %>% 
  ggplot(aes(RY, RA, fill = interaction(A < 0.5, Y < 0.5))) +
  geom_tile() + 
  facet_grid(fSPA~fSPY, labeller = label_parsed) + 
  scale_x_continuous(limits = c(min(dp2$RY), max(dp2$RY)), expand = c(0,0),
                     breaks = c(0,0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) + 
  scale_y_continuous(limits = c(min(dp2$RA), max(dp2$RA)), expand = c(0,0),
                     breaks = c(0,0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) + 
  scale_fill_viridis(discrete = T, option = "magma", begin = 0.9, end = 0.1,
                     labels = c(bquote(italic(s)[12]),
                                bquote(italic(s)[22]),
                                bquote(italic(s)[12]*" & "*italic(s)[22]))) + 
  labs(x = bquote(italic(r)[1]), y = bquote(italic(r)[2]), fill = "SD gene") + 
  
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 10))


text_left <- ggplot(data = NULL, aes(1, 1)) + 
  geom_text(label = expression("  "~italic(P)[1]~"parentally-antagonistic,"~italic(P)[2]~"sexually-antagonistic")) + 
  theme_nothing()

text_right <- ggplot(data = NULL, aes(1, 1)) + 
  geom_text(label = expression("  "~italic(P)[1]~"sexually-antagonistic,"~italic(P)[2]~"parentally-antagonistic")) + 
  theme_nothing()

text_top <- ggplot(data = NULL, aes(1, 1)) + 
  geom_text(label = expression("  Male-determining "~italic(s)[22]), angle = 90) + 
  theme_nothing()

text_bottom <- ggplot(data = NULL, aes(1, 1)) + 
  geom_text(label = expression("  Female-determining "~italic(s)[22]), angle = 90) + 
  theme_nothing()

topleft_filler <- ggplot(data = NULL, aes(1, 1)) + 
  geom_text(label = "") + 
  theme_nothing()

legend <- get_legend(pa)

# p <- plot_grid(pa + theme(legend.position = "none"),
#               
#                pc + theme(legend.position = "none"),
#                pd + theme(legend.position = "none"),
#                labels = c("A", "B", "C", "D"),
#                ncol = 2, nrow = 2)

p <- plot_grid(topleft_filler, text_left, text_right, 
          text_top, pa + theme(legend.position = "none"),  pb + theme(legend.position = "none"),
          text_bottom,  pc + theme(legend.position = "none"),  pd + theme(legend.position = "none"),
          
          ncol = 3, nrow = 3, rel_widths = c(0.1, 2, 2), rel_heights = c(0.1, 2, 2))

p2 <- plot_grid(p, legend, rel_widths = c(2, 0.2))


png(file = "2024_07_12_Figure_4.png", width = 9, height = 8.5, units = "in", res = 1000)
p2
dev.off()


pdf(file = "2024_07_12_Figure_4.pdf", width = 9, height = 8.5)
p2
dev.off()

