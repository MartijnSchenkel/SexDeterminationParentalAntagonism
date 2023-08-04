library(tidyverse)
library(viridis)
library(mgcv)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# ======================= #
# ---- Data processing ----
# ======================= #

# Specify simulation date + which sets have been simulated.
date <- "2023_07_14"
epsilon <- 1e-3

CY <- seq(-1.1, -1.9, -0.2)
CA <- CY


# Read first batch
dr <- read.table(str_c(date,"_CY_-1.1_CA_-1.1/results_summarized.txt"), T) %>% 
  as_tibble() %>% filter(!MaleDeterminer %in% unique(MaleDeterminer))
dt <- dr
# Read remaining batches
for(i in CY) for(j in CA)
{
  filename <- str_c(date,"_CY_", i,"_CA_",j,"/results_summarized.txt")
  dt <- read.table(filename, T) %>% 
    as_tibble()
  
  dr <- dr %>% full_join(dt, by = names(dt))
}
write.table(dr, file = "results_masculinizer_combined.txt", col.names = T, row.names = F)


# Read first batch
dr <- read.table(str_c(date,"_Feminizer_CY_-1.1_CA_-1.1/results_summarized.txt"), T) %>% 
  as_tibble() %>% filter(!MaleDeterminer %in% unique(MaleDeterminer))
dt <- dr
# Read remaining batches
for(i in CY) for(j in CA)
{
  filename <- str_c(date,"_Feminizer_CY_", i,"_CA_",j,"/results_summarized.txt")
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

# Plot raw data on Y frequencies
dr %>% ggplot(aes(AY, AA, color = MXYPat)) + geom_point() + facet_grid(CY ~ CA)

# Filter out runs with mixed systems (defined as both A and Y have frequency over >0.001)
dr3 <- dr %>% filter(MXYPat < epsilon | MXYPat > (1 - epsilon)) %>% 
  mutate(Turnover = 1 - (round(MXYPat) == 1)) %>% 
  as_tibble()

# Plot raw data on turnover rates
dr3 %>% ggplot(aes(AY, AA, color = Turnover)) + geom_point() + facet_grid(CY ~ CA)

# Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
mM <- gam(Turnover ~ t2(AY, AA, by = interaction(factor(CY), factor(CA)), 
                    k = 3, m = 2, bs = "ts"), 
          data = dr3, 
          full = T, family = "binomial")

# Check GAM for oddities
summary(mM)
gam.check(mM) # no oddities

# Save GAM
saveRDS(object = mM, file = "2023_07_14_GAM_Masculinizer_CY_CA_AY_AA.RDS")
mM <- readRDS("2023_07_14_GAM_Masculinizer_CY_CA_AY_AA.RDS")

# Generate data set to obtain predicted values across parameter space.
dp <- expand.grid(CY = factor(CY), CA = factor(CA),
                  AY = seq(-1, 0, length.out = 300), AA = seq(-1, 0, length.out = 300))

# Add predicted values from model
dpp <- dp %>% transform(Turnover = predict(mM, newdata = dp, type = "response"))

# Add nice facet strip labels
dpp <- dpp %>% dplyr::mutate(fCY = str_c("italic(c)[XY]==", CY),
                     fCA = str_c("italic(c)[A]==", CA))

write.table(x = dpp, col.names = T, row.names = F,
            file = "2023_07_14_Masculinizer_CY_CA_AY_AA_predictions.txt")
dpp <- read.table("2023_07_14_Masculinizer_CY_CA_AY_AA_predictions.txt", T)
summary(dpp)

# ========================================= #
# ---- Data visualization - Masculinizer ----
# ========================================= #

dpp %>% as_tibble() %>% filter(Turnover > 0.5) %>% 
  group_by(AY, AA, CY) %>% summarize(mCA = max(CA))



pm <- dpp %>% filter(CA < -1.1) %>% ggplot(aes(AY, AA, z = Turnover)) + 
  geom_contour_filled(aes(fill = factor(CA)), breaks = c(0,0.5)) + 
  facet_wrap(~fCY, labeller = label_parsed, ncol = 5) + 
  scale_x_continuous(limits = c(min(dpp$AY), max(dpp$AY)), expand = c(0,0),
                     breaks = c(-1,-0.75, -0.5, -0.25, 0),
                     labels = c("-1", "-0.75", "-0.5", "-0.25", "0")) + 
  scale_y_continuous(limits = c(min(dpp$AA), max(dpp$AA)), expand = c(0,0)) + 
  labs(x = bquote(italic(a)[XY]), y = bquote(italic(a)[A]), fill = bquote(italic(c)[A])) + 
  scale_fill_viridis(discrete = T, option = "magma", begin = 0.9, end = 0.1) + 
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"))

png(file = "2023_07_14_Masculinizer_CY_CA_AY_AA.png", width = 8, height = 2, units = "in", res = 1000)
pm
dev.off()


pdf(file = "2023_07_14_Masculinizer_CY_CA_AY_AA.pdf", width = 8, height = 2)
pm
dev.off()

# ================================= #
# ---- Data analysis - Feminizer ----
# ==================================== #

# Read masculinizer data
dr <- read.table("results_feminizer_combined.txt", T)

# Plot raw data on Y frequencies
dr %>% ggplot(aes(AY, AA, color = MXYPat)) + geom_point() + facet_grid(CY ~ CA)

# Filter out runs with mixed systems (defined as both A and Y have frequency over >0.001)
dr3 <- dr %>% filter(MAMat < epsilon | MAMat > (1 - epsilon)) %>% 
  mutate(Turnover = (MAMat == 1)) %>% 
  as_tibble()

# Plot raw data on turnover rates
dr3 %>% ggplot(aes(AY, AA, color = Turnover)) + geom_point() + facet_grid(CY ~ CA)

# Generate GAM with specific beta estimates for EACH combination of c_XY and c_A
mF <- gam(Turnover ~ t2(AY, AA, by = interaction(factor(CY), factor(CA)), 
                        k = 3, m = 2, bs = "ts"), 
          data = dr3, 
          full = T, family = "binomial")

# Check GAM for oddities
summary(mF)
gam.check(mF) # no oddities

# Save GAM
saveRDS(object = mF, file = "2023_07_14_GAM_Feminizer_CY_CA_AY_AA.RDS")
mF <- readRDS("2023_07_14_GAM_Feminizer_CY_CA_AY_AA.RDS")

# Generate data set to obtain predicted values across parameter space.
dp <- expand.grid(CY = factor(CY), CA = factor(CA),
                  AY = seq(-1, 0, length.out = 300), AA = seq(-1, 0, length.out = 300))

# Add predicted values from model
dpp <- dp %>% transform(Turnover = predict(mF, newdata = dp, type = "response"))

# Add nice facet strip labels
dpp <- dpp %>% dplyr::mutate(fCY = str_c("italic(c)[XY]==", CY),
                             fCA = str_c("italic(c)[A]==", CA))

write.table(x = dpp, col.names = T, row.names = F,
            file = "2023_07_14_Feminizer_CY_CA_AY_AA_predictions.txt")
dpp <- read.table("2023_07_14_Feminizer_CY_CA_AY_AA_predictions.txt", T)
summary(dpp)

# ========================================= #
# ---- Data visualization - Feminizer ----
# ========================================= #

dpp %>% as_tibble() %>% filter(Turnover > 0.5) %>% 
  group_by(AY, AA, CY) %>% summarize(mCA = max(CA))



pm <- dpp %>% ggplot(aes(AY, AA)) + 
  geom_tile(aes(fill = factor(Turnover))) + 
  facet_grid(fCA~fCY, labeller = label_parsed) + 
  scale_x_continuous(limits = c(min(dpp$AY), max(dpp$AY)), expand = c(0,0),
                     breaks = c(-1,-0.75, -0.5, -0.25, 0),
                     labels = c("-1", "-0.75", "-0.5", "-0.25", "0")) + 
  scale_y_continuous(limits = c(min(dpp$AA), max(dpp$AA)), expand = c(0,0)) + 
  labs(x = bquote(italic(a)[XY]), y = bquote(italic(a)[A]), fill = "SD gene" ) + 
  scale_fill_viridis(discrete = T, option = "magma", begin = 0.9, end = 0.1, 
                     labels = c(bquote(italic(A)),bquote(italic(Y)),bquote(italic(Y)*"and"*(italic(A))))) + 
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white")) + 
  guides(fill = guide_legend(labels = c(letters[1:3])))
pm

png(file = "2023_07_14_Feminizer_CY_CA_AY_AA.png", width = 9, height = 8, units = "in", res = 1000)
pm
dev.off()


pdf(file = "2023_07_14_Feminizer_CY_CA_AY_AA.pdf", width = 9, height = 8)
pm
dev.off()