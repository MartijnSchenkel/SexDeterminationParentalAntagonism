library(tidyverse)
library(viridis)
library(mgcv)
library(cowplot)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# Specify simulation date + which sets have been simulated.
date <- "2024_07_02"
epsilon <- 1e-3

SPY <- seq(0.01, 0.05, 0.01)
SPA <- SPY

# =========================================== #
# ---- Data processing - Sexual antagonism ----
# =========================================== #

# Read first batch
dr <- read.table(str_c(date,"_Mas_SS_SPY_0.01_SPA_0.01/results_summarized.txt"), T) %>% 
  as_tibble() %>% filter(!MaleDeterminer %in% unique(MaleDeterminer))
dt <- dr
# Read remaining batches
for(i in SPY) for(j in SPA)
{
  filename <- str_c(date,"_Mas_SS_SPY_", i,"_SPA_",j,"/results_summarized.txt")
  dt <- read.table(filename, T) %>% 
    as_tibble()
  
  dr <- dr %>% full_join(dt, by = names(dt))
}
write.table(dr, file = "results_masculinizer_SS_combined.txt", col.names = T, row.names = F)


# Read first batch
dr <- read.table(str_c(date,"_Fem_SS_SPY_0.01_SPA_0.01/results_summarized.txt"), T) %>% 
  as_tibble() %>% filter(!MaleDeterminer %in% unique(MaleDeterminer))
dt <- dr
# Read remaining batches
for(i in SPY) for(j in SPA)
{
  filename <- str_c(date,"_Fem_SS_SPY_", i,"_SPA_",j,"/results_summarized.txt")
  dt <- read.table(filename, T) %>% 
    as_tibble()
  
  dr <- dr %>% full_join(dt, by = names(dt))
}
write.table(dr, file = "results_feminizer_SS_combined.txt", col.names = T, row.names = F)

# ========================================= #
# ---- Data analysis - SA - Masculinizer ----
# ========================================= #

# Read masculinizer data
dr <- read.table("results_masculinizer_SS_combined.txt", T)

# # Plot raw data on Y frequencies
# dr %>% ggplot(aes(AY, AA, color = MXYPat)) + geom_point() + facet_grid(CY ~ CA)

# Filter out runs with mixed systems (defined as both A and Y have frequency over >0.001)
dr3 <- dr %>%
  mutate(Y = ifelse(MXYPat > 0.1, T, F),
         A = ifelse(MAPat > 0.1, T, F)) %>% 
  as_tibble()


dr3 %>% ggplot(aes(RY, RA, color = interaction(Y, A))) + geom_point() + facet_grid(SPY ~ SPA)
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
saveRDS(object = mMY, file = str_c(date, "_GAM_Masculinizer_SS_SPY_SPA_RY_RA_MY.RDS"))
mMY <- readRDS(str_c(date, "_GAM_Masculinizer_SS_SPY_SPA_RY_RA_MY.RDS"))

# Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
mMA <- gam(A ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)), 
                  k = 3, m = 2, bs = "ts"), 
           data = dr3, 
           full = T, family = "binomial")

# Check GAM for oddities
summary(mMA)
gam.check(mMA) # no oddities

# Save GAM
saveRDS(object = mMA, file = str_c(date, "_GAM_Masculinizer_SS_SPY_SPA_RY_RA_MA.RDS"))
mMA <- readRDS(str_c(date, "_GAM_Masculinizer_SS_SPY_SPA_RY_RA_MA.RDS"))



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
            file = str_c(date, "_Masculinizer_SS_SPY_SPA_RY_RA_predictions.txt"))
dpp <- read.table(str_c(date, "_Masculinizer_SS_SPY_SPA_RY_RA_predictions.txt"), T)
summary(dpp)

# ====================================== #
# ---- Data analysis - SA - Feminizer ----
# ====================================== #

# Read feminizer data
dr <- read.table("results_feminizer_SS_combined.txt", T)

# Plot raw data on Y frequencies
dr %>% ggplot(aes(RY, RA, color = FAMat)) + geom_point() + facet_grid(SPY ~ SPA)


dr3 <- dr %>%
  mutate(Y = ifelse(FAMat < 0.1, T, F),
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
saveRDS(object = mMY, file = str_c(date, "_GAM_Feminizer_SS_SPY_SPA_RY_RA_Y.RDS"))
mMY <- readRDS(str_c(date, "_GAM_Feminizer_SS_SPY_SPA_RY_RA_Y.RDS"))

# Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
mMA <- gam(A ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)), 
                  k = 3, m = 2, bs = "ts"), 
           data = dr3, 
           full = T, family = "binomial")

# Check GAM for oddities
summary(mMA)
gam.check(mMA) # no oddities

# Save GAM
saveRDS(object = mMA, file = str_c(date, "__GAM_Feminizer_SPY_SPA_RY_RA_A.RDS"))
mMA <- readRDS(str_c(date, "_GAM_Feminizer_SPY_SPA_RY_RA_A.RDS"))

# Generate data set to obtain predicted values across parameter space.
dp <- expand.grid(SPY = factor(SPY), SPA = factor(SPA),
                  RY = seq(0, 0.5, length.out = 300), RA = seq(0, 0.5, length.out = 300))

# Add predicted values from model
dpp <- dp %>% transform(Y = predict(mMY, newdata = dp, type = "response"),
                        A = predict(mMA, newdata = dp, type = "response"))

# # Add predicted values from model
# dpp <- dp %>% transform(Turnover = predict(mF, newdata = dp, type = "response"))

# Add nice facet strip labels
dpp <- dpp %>% dplyr::mutate(fSPY = str_c("italic(u)[1]==", SPY),
                             fSPA = str_c("italic(u)[2]==", SPA))

write.table(x = dpp, col.names = T, row.names = F,
            file = str_c(date, "_Feminizer_SS_SPY_SPA_RY_RA_predictions.txt"))
dpp <- read.table(str_c(date, "_Feminizer_SS_SPY_SPA_RY_RA_predictions.txt"), T)
summary(dpp)

# ============================================ #
# ---- REPEAT FOR PARENTAL ANTAGONISM PLOTS ----
# ============================================ #

# Specify simulation date + which sets have been simulated.
date <- "2023_10_20"
epsilon <- 1e-3

SPY <- seq(0.01, 0.05, 0.01)
SPA <- SPY

# ======================= #
# ---- Data processing ----
# ======================= #

# Read first batch
dr <- read.table(str_c(date,"_Mas_SPY_0.01_SPA_0.01/results_summarized.txt"), T) %>% 
  as_tibble() %>% filter(!MaleDeterminer %in% unique(MaleDeterminer))
dt <- dr
# Read remaining batches
for(i in SPY) for(j in SPA)
{
  filename <- str_c(date,"_Mas_SPY_", i,"_SPA_",j,"/results_summarized.txt")
  dt <- read.table(filename, T) %>% 
    as_tibble()
  
  dr <- dr %>% full_join(dt, by = names(dt))
}
write.table(dr, file = "results_masculinizer_combined.txt", col.names = T, row.names = F)


# Read first batch
dr <- read.table(str_c(date,"_Fem_SPY_0.01_SPA_0.01/results_summarized.txt"), T) %>% 
  as_tibble() %>% filter(!MaleDeterminer %in% unique(MaleDeterminer))
dt <- dr
# Read remaining batches
for(i in SPY) for(j in SPA)
{
  filename <- str_c(date,"_Fem_SPY_", i,"_SPA_",j,"/results_summarized.txt")
  dt <- read.table(filename, T) %>% 
    as_tibble()
  
  dr <- dr %>% full_join(dt, by = names(dt))
}
write.table(dr, file = "results_feminizer_combined.txt", col.names = T, row.names = F)

# ========================================= #
# ---- Data analysis - SA - Masculinizer ----
# ========================================= #

# Read masculinizer data
dr <- read.table("results_masculinizer_combined.txt", T)

# Filter out runs with mixed systems (defined as both A and Y have frequency over >0.001)
dr3 <- dr %>%
  mutate(Y = ifelse(MXYPat > 0.1, T, F),
         A = ifelse(MAPat > 0.1, T, F)) %>% 
  as_tibble()


dr3 %>% ggplot(aes(RY, RA, color = interaction(Y, A))) + geom_point() + facet_grid(SPY ~ SPA)

# Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
mMY <- gam(Y ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)),
                    k = 3, m = 2, bs = "ts"),
          data = dr3,
          full = T, family = "binomial")

# Check GAM for oddities
summary(mMY)
gam.check(mMY) # no oddities

# Save GAM
saveRDS(object = mMY, file = str_c(date, "_GAM_Masculinizer_SPY_SPA_RY_RA_MY.RDS"))
mMY <- readRDS(str_c(date, "_GAM_Masculinizer_SPY_SPA_RY_RA_MY.RDS"))



# Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
mMA <- gam(A ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)), 
                  k = 3, m = 2, bs = "ts"), 
           data = dr3, 
           full = T, family = "binomial")

# Check GAM for oddities
summary(mMA)
gam.check(mMA) # no oddities

# Save GAM
saveRDS(object = mMA, file = str_c(date, "_GAM_Masculinizer_SPY_SPA_RY_RA_MA.RDS"))
mMA <- readRDS(str_c(date, "_GAM_Masculinizer_SPY_SPA_RY_RA_MA.RDS"))

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
            file = str_c(date, "_Masculinizer_SPY_SPA_RY_RA_predictions.txt"))
dpp <- read.table(str_c(date, "_Masculinizer_SPY_SPA_RY_RA_predictions.txt"), T)
summary(dpp)

# ================================= #
# ---- Data analysis - Feminizer ----
# ==================================== #

# Read feminizer data
dr <- read.table("results_feminizer_combined.txt", T)

# Plot raw data on Y frequencies
dr %>% ggplot(aes(RY, RA, color = FAMat)) + geom_point() + facet_grid(SPY ~ SPA)

# Filter out runs with mixed systems (defined as both A and Y have frequency over >0.001)
# dr3 <- dr %>% filter(FAMat < epsilon | FAMat > (1 - epsilon)) %>% 
#   mutate(Turnover = (FAMat == 1)) %>% 
#   as_tibble()

dr3 <- dr %>%
  mutate(Y = ifelse(FAMat < 0.1, T, F),
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
saveRDS(object = mMY, file = str_c(date, "_GAM_Feminizer_SS_SPY_SPA_RY_RA_Y.RDS"))
mMY <- readRDS(str_c(date, "_GAM_Feminizer_SPY_SPA_RY_RA_Y.RDS"))

# Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
mMA <- gam(A ~ t2(RY, RA, by = interaction(factor(SPY), factor(SPA)), 
                  k = 3, m = 2, bs = "ts"), 
           data = dr3, 
           full = T, family = "binomial")

# Check GAM for oddities
summary(mMA)
gam.check(mMA) # no oddities

# Save GAM
saveRDS(object = mMA, file = str_c(date, "_GAM_Feminizer_SS_SPY_SPA_RY_RA_A.RDS"))
mMA <- readRDS(str_c(date, "_GAM_Feminizer_SPY_SS_SPA_RY_RA_A.RDS"))

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
            file = str_c(date, "_Feminizer_SPY_SPA_RY_RA_predictions.txt"))
dpp <- read.table(str_c(date, "_Feminizer_SPY_SPA_RY_RA_predictions.txt"), T)
summary(dpp)


# ================ #
# ---- Figure 2 ----
# ================ #

dppfp <- read.table("2023_10_20_Feminizer_SPY_SPA_RY_RA_predictions.txt", T)
dppmp <- read.table("2023_10_20_Masculinizer_SPY_SPA_RY_RA_predictions.txt", T)

dppfs <- read.table("2024_07_02_Feminizer_SS_SPY_SPA_RY_RA_predictions.txt", T)
dppms <- read.table("2024_07_02_Masculinizer_SS_SPY_SPA_RY_RA_predictions.txt", T)


dpf2p <- dppfp %>% dplyr::mutate(MaleDeterminer = F, Type = "Parental antagonism")
dpm2p <- dppmp %>% dplyr::mutate(MaleDeterminer = T, Type = "Parental antagonism")


dpf2s <- dppfs %>% dplyr::mutate(MaleDeterminer = F, Type = "Sexual antagonism")
dpm2s <- dppms %>% dplyr::mutate(MaleDeterminer = T, Type = "Sexual antagonism")

dp <- dpf2p %>% full_join(dpm2p) %>% full_join(dpf2s) %>% full_join(dpm2s)

dp2 <- dp %>% filter(SPY %in% c(0.01, 0.03, 0.05), SPA %in% c(0.01, 0.03, 0.05))

# Set up part A
pa <- dp2 %>% filter(MaleDeterminer == T, Type == "Parental antagonism") %>% 
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
        strip.text = element_text(color = "white"))


# Set up part B (note: logic rules here are inverted due to different GAM setup for feminizing A)
pb <- dp2 %>% filter(MaleDeterminer == F, Type == "Parental antagonism") %>%  filter(Y < 0.5 | A < 0.5) %>% 
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
        strip.text = element_text(color = "white"))

# Set up part C
pc <- dp2 %>% filter(MaleDeterminer == T, Type == "Sexual antagonism") %>% 
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
        strip.text = element_text(color = "white"))


# Set up part D (note: logic rules here are inverted due to different GAM setup for feminizing A)
pd <- dp2 %>% filter(MaleDeterminer == F, Type == "Sexual antagonism") %>%  filter(Y < 0.5 | A < 0.5) %>% 
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
        strip.text = element_text(color = "white"))


legend <- get_legend(pa)

p1 <- plot_grid(pa + theme(legend.position = "none"),
          pb + theme(legend.position = "none"),
          pc + theme(legend.position = "none"),
          pd + theme(legend.position = "none"),
          ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))

p2 <- plot_grid(p1, legend,
          rel_widths = c(1,0.1), ncol = 2)



text_left <- ggplot(data = NULL, aes(1, 1)) + 
  geom_text(label = expression("  "~italic(P)[1]~"and"~italic(P)[2]~"parentally-antagonistic")) + 
  theme_nothing()

text_right <- ggplot(data = NULL, aes(1, 1)) + 
  geom_text(label = expression("  "~italic(P)[1]~"and"~italic(P)[2]~"sexually-antagonistic")) + 
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
               text_top, pa + theme(legend.position = "none"),  pc + theme(legend.position = "none"),
               text_bottom,  pb + theme(legend.position = "none"),  pd + theme(legend.position = "none"),
               
               ncol = 3, nrow = 3, rel_widths = c(0.1, 2, 2), rel_heights = c(0.1, 2, 2))

p2 <- plot_grid(p, legend, rel_widths = c(2, 0.2))






png(file = str_c(date, "_Figure_2.png"), width = 9, height = 8.5, units = "in", res = 1000)
p2
dev.off()


pdf(file = str_c(date, "_Figure_2.pdf"), width = 9, height = 8.5)
p2
dev.off()

