# ----------------------- #
# Administrative setup ----
# ----------------------- #
library(tidyverse)

# NOTE: This script performs 5,000 simulations of the possible invasion of a 
# novel male-determining allele. The parameter values are set below under
# "Parameter values". Note that the values for r_XY (RY) and r_A (RA) are
# set at a random value between 0 and 0.5; this is because this script is used
# to determine the effect of S_XY-P_XY and S_A-P_A recombination on the spread
# of a male-determining A.

# setwd(dirname(rstudioapi::getActiveDocumentContext()))

simdate <- "2023_07_19_Masculinizer"
simvars <- "SPY_SPA"

simset <- str_c(simdate, "_", simvars)

dir.create("all_sims")


ds <- tibble(FAMat = numeric(),
             FAPat = numeric(),
             FXYMat = numeric(),
             FXYPat = numeric(),
             MAMat = numeric(),
             MAPat = numeric(),
             MXYMat = numeric(),
             MXYPat = numeric(),
             SPY = numeric(),
             SPA = numeric(),
             AY = numeric(),
             AA = numeric(),
             BY = numeric(),
             BA = numeric(),
             RY = numeric(),
             RA = numeric(),
             CY = numeric(),
             CA = numeric(),

             MaleDeterminer = numeric())

for(nsim in 1:5000)
{
# ------------------- #
# Parameter values ---- 
# ------------------- #
  
# Randomized parameters
RY <- runif(min = 0, max = 0.5, n = 1)
RA <- runif(min = 0, max = 0.5, n = 1)
  

# Global parameters
max_t <- 50000
gen_skip <- 100

# De novo SD mutation parameters
mu_gen <- 10000
muSD <- 0.001

# Selection effect parameters
s <- c(0.01, 0.01)
# a <- c(-0.1, -0.3)
# b <- c(3.11, 3.6)

CY <- 1.5
CA <- 1.5

a <- c(-0.1, -0.9)
b <- c(-CY * a[1] + 2, -CA * a[2] + 2)


# Recombination rates
r <- c(RY, RA)

# Novel SD gene is male-determiner (T) or female-determiner (F)
NovelSexDeterminer <- T

# Frequency of focal pOne at all loci for [1] maternally and [2] paternally inherited copies

# Allele order
# 1/5: Y mat/pat
# 2/6: A mat/pat
# 3/7: PY mat/pat
# 4/8: PA mat/pat

initial_frequencies <- c(0, 0, 0.5, 0.5, 
                         0.5, 0.0, 0.5, 0.5)

pOne <- array(data = initial_frequencies, dim = c(length(r),2,2))



# -------------------- #
# Model preparation ----
# -------------------- #

# Fitness arrays
# Configurations:
# w[1,1] = mat 1, pat 1
# w[2,1] = mat 2, pat 1
# w[1,2] = mat 1, pat 2
# w[2,2] = mat 2, pat 2
wI <- array(data = c(1, 1+a[1]*s[1], 1+b[1]*s[1], 1+s[1]), dim = c(2,2))
wII <- array(data = c(1, 1+a[2]*s[2], 1+b[2]*s[2], 1+s[2]), dim = c(2,2))

# Calculate haplotype 
pHaplo <- array(dim=c(length(r),4,4,4))
for(i in 1:length(r))
{
  pHaplo[i,1,1,] <- c(1, 0, 0, 0)
  pHaplo[i,1,2,] <- c(1/2, 1/2, 0, 0) 
  pHaplo[i,1,3,] <- c(1/2, 0, 1/2, 0)
  pHaplo[i,1,4,] <- c((1-r[i])/2, r[i]/2, r[i]/2, (1-r[i])/2)
  pHaplo[i,2,1,] <- c(1/2, 1/2, 0, 0)
  pHaplo[i,2,2,] <- c(0, 1, 0, 0)
  pHaplo[i,2,3,] <- c(r[i]/2, (1-r[i])/2, (1-r[i])/2, r[i]/2)
  pHaplo[i,2,4,] <- c(0, 1/2, 0, 1/2)
  pHaplo[i,3,1,] <- c(1/2, 0, 1/2, 0)
  pHaplo[i,3,2,] <- c(r[i]/2, (1-r[i])/2, (1-r[i])/2, r[i]/2)
  pHaplo[i,3,3,] <- c(0, 0, 1, 0)
  pHaplo[i,3,4,] <- c(0, 0, 1/2, 1/2)
  pHaplo[i,4,1,] <- c((1-r[i])/2, r[i]/2, r[i]/2, (1-r[i])/2)
  pHaplo[i,4,2,] <- c(0, 1/2, 0, 1/2)
  pHaplo[i,4,3,] <- c(0, 0, 1/2, 1/2)
  pHaplo[i,4,4,] <- c(0, 0, 0, 1)
}


# Create gametogenesis array
G2G <- array(dim=c(4,4,4,4,4,4))
for(i in 1:4)
{
  for(j in 1:4)
  {
    for(k in 1:4)
    {
      for(l in 1:4)
      {
        for(m in 1:4)
        {
          for(n in 1:4)
          {
            G2G[i,j,k,l,m,n] <- 
              pHaplo[1,i,k,m] * 
              pHaplo[2,j,l,n]
          }
        }
      }
    }
  }
}


# Function to extract gamete frequencies from G2G array.
MakeGametes <- function(G){
  A[G[1], G[2], G[3], G[4]] * G2G[G[1], G[2], G[3], G[4],,]
}

# Mutagenesis function (introduce novel SD gene on II)
mutate <- function(G)
{
  G[,3] <- G[,1] * muSD
  G[,1] <- G[,1] * (1 - muSD)
  G[,4] <- G[,2] * muSD
  G[,2] <- G[,2] * (1 - muSD)
  G
}


# ----------------------- #
# Model initialization ----
# ----------------------- #

# Calculate initial haplotype frequencies
GS <- array(dim = c(length(r),2,4))
for(i in 1:length(r))
{
  GS[i,1,] <- c((1 - pOne[i,1,1]) * (1 - pOne[i,2,1]), 
                (1 - pOne[i,1,1]) * (pOne[i,2,1]), 
                (pOne[i,1,1]) * (1 - pOne[i,2,1]), 
                (pOne[i,1,1]) * (pOne[i,2,1]))
  GS[i,2,] <- c((1 - pOne[i,1,2]) * (1 - pOne[i,2,2]),
                (1 - pOne[i,1,2]) * (pOne[i,2,2]), 
                (pOne[i,1,2]) * (1 - pOne[i,2,2]), 
                (pOne[i,1,2]) * (pOne[i,2,2]))
}

# Set up initial population and count number of focal alleles per genotype.
A <- array(dim = c(4, 4, 4, 4))
nOneM <- array(dim=c(4, 4, 4, 4, 4)) # Number of focal alleles in maternal haplotype
nOneP <- nOneM # number of focal alleles in paternal haplotype

nSD <- c(0,0,1,1)
nSA <- c(0,1,0,1)

# Calculates fitness based on genotypes at PY and PA separately, to be combined later..
wIP <- array(dim = c(4,4,4,4))
wIIP <- array(dim = c(4,4,4,4))

for(i in 1:4)
{
  for(j in 1:4) 
  {
    for(k in 1:4)
    {
      for(l in 1:4)
      {
        # Initialize population.
        A[i,j,k,l] <- GS[1,1,i] * GS[2,1,j] * GS[1,2,k] * 
          GS[2,2,l]
      
        # Count number of focal alleles per genotype.
        nOneM[1,i,j,k,l] <- nSD[i] # freq of Y on maternal haplotype
        nOneM[2,i,j,k,l] <- nSA[i] # freq of PY on maternal haplotype
        nOneM[3,i,j,k,l] <- nSD[j] # freq of A on maternal haplotype
        nOneM[4,i,j,k,l] <- nSA[j] # freq of PA on maternal haplotype
        nOneP[1,i,j,k,l] <- nSD[k]
        nOneP[2,i,j,k,l] <- nSA[k]
        nOneP[3,i,j,k,l] <- nSD[l]
        nOneP[4,i,j,k,l] <- nSA[l]
        wIP[i,j,k,l] <- wI[nOneM[2,i,j,k,l] + 1, nOneP[2,i,j,k,l] + 1]
        wIIP[i,j,k,l] <- wII[nOneM[4,i,j,k,l] + 1, nOneP[4,i,j,k,l] + 1]
        
      }
    }
  }
}
nOne <- nOneM + nOneP 
nOneM[4,4,,4,]




# Create arrays with sex determination info
if(NovelSexDeterminer)
{
  IsMale <- (nOne[1,,,,] > 0 | nOne[3,,,,] > 0)
} else {
  IsMale <- (nOne[1,,,,] > 0 & nOne[3,,,,] == 0)
}
IsFemale <- 1 - IsMale


# Calculate fitness values of all genotypes.
W <- wIP * wIIP

# min(wIP);max(wIP)
# min(wIIP);max(wIIP)
# min(W);max(W)

# List allele type info
Alleles <- c("Y", "PY", "A", "PA")
LinkageGroup <- c("I", "I", "II", "II")
# Type <- c("SD", "PA", "SD", "PA")

# # Generate data frame for writing results.
# d <- tibble(Generation = numeric(),
#             SexRatio =  numeric(),
#             FreqMatF = numeric(),
#             FreqPatF = numeric(),
#             FreqMatM = numeric(),
#             FreqPatM = numeric(),
#             Chromosome = character(),
#             Allele = character())
# 
# # Write initial generation
# for(l in 1:4)
# {
#   dt <- tibble(Generation = 0,
#                SexRatio = A %*% IsMale,
#                FreqMatF = (A * IsFemale) %*% nOneM[l,,,,] / (1 - SexRatio),
#                FreqPatF = (A * IsFemale) %*% nOneP[l,,,,] / (1 - SexRatio),
#                FreqMatM = (A * IsMale) %*% nOneM[l,,,,] / SexRatio,
#                FreqPatM = (A * IsMale) %*% nOneP[l,,,,] / SexRatio,
#                Chromosome = LinkageGroup[l],
#                Allele = Alleles[l])
#   d <- d %>% full_join(dt, by = c("Generation", "SexRatio", 
#                                   "FreqMatF", "FreqPatF", 
#                                   "FreqMatM", "FreqPatM", 
#                                   "Chromosome", "Allele"))
# }

# ---------------- #
# ---- THE LOOP ----
# ---------------- #
for(t in 1:max_t)
{
  # Determine which genotypes have non-zero frequency and are female c.q. male.
  Females <- which(A > 0 & IsFemale == T,arr.ind=T,useNames=F)
  Males <- which(A > 0 & IsMale == T,arr.ind=T,useNames=F)

  # Define gamete vectors
  Gf <- rep(0, 4 * 4)  
  Gm <- rep(0, 4 * 4)
  
  
  # Make gametes for both sexes separately.
  for(i in 1:nrow(Females))
  {
    Gf <- Gf + MakeGametes(Females[i,])
  }
  
  for(i in 1:nrow(Males))
  {
    Gm <- Gm + MakeGametes(Males[i,])
  }
  
  # At generation mu_gen, mutate male and female gametes to introduce A
  if(t == mu_gen)
  {
    Gf <- mutate(Gf)
    Gm <- mutate(Gm)
  }
  
  # Reproduce & normalize
  J <- Gf %o% Gm # Random fusion
  J <- J * W     # Selection
  J <- J/sum(J)  # Normalization
  A <- J         # Maturation (non-overlapping generations)
  
  
  # Write output every gen_skip generations.
  if(t %% gen_skip == 0)
  {

  }
  
}

# # Write output set 1
# for(l in 1:4)
# {
#   dt <- tibble(Generation = t,
#                SexRatio = A %*% IsMale,
#                FreqMatF = (A * IsFemale) %*% nOneM[l,,,,] / (1 - SexRatio),
#                FreqPatF = (A * IsFemale) %*% nOneP[l,,,,] / (1 - SexRatio),
#                FreqMatM = (A * IsMale) %*% nOneM[l,,,,] / SexRatio,
#                FreqPatM = (A * IsMale) %*% nOneP[l,,,,] / SexRatio,
#                Chromosome = LinkageGroup[l],
#                Allele = Alleles[l])
#   d <- d %>% full_join(dt, by = c("Generation", "SexRatio", 
#                                   "FreqMatF", "FreqPatF", 
#                                   "FreqMatM", "FreqPatM", 
#                                   "Chromosome", "Allele"))
# }



nHaps <- array(dim = c(length(r), 2, 2 * 2, 4,4,4,4))
# dimensions = LG, origin, haplotype variant, 
# haplotype maternal XY, haplotype maternal A, haplotype paternal XY, haplotype paternal A 


Haplotype <- c("00", "01", "10", "11")
LG <- c("XY", "A")
Origin <- c("Maternal", "Paternal")
nHaps[1,1,1,,,,] <- nOneM[1,,,,] == F & nOneM[2,,,,] == F
nHaps[1,1,2,,,,] <- nOneM[1,,,,] == F & nOneM[2,,,,] == T
nHaps[1,1,3,,,,] <- nOneM[1,,,,] == T & nOneM[2,,,,] == F
nHaps[1,1,4,,,,] <- nOneM[1,,,,] == T & nOneM[2,,,,] == T

nHaps[1,2,1,,,,] <- nOneP[1,,,,] == F & nOneP[2,,,,] == F
nHaps[1,2,2,,,,] <- nOneP[1,,,,] == F & nOneP[2,,,,] == T
nHaps[1,2,3,,,,] <- nOneP[1,,,,] == T & nOneP[2,,,,] == F
nHaps[1,2,4,,,,] <- nOneP[1,,,,] == T & nOneP[2,,,,] == T

nHaps[2,1,1,,,,] <- nOneM[3,,,,] == F & nOneM[4,,,,] == F
nHaps[2,1,2,,,,] <- nOneM[3,,,,] == F & nOneM[4,,,,] == T
nHaps[2,1,3,,,,] <- nOneM[3,,,,] == T & nOneM[4,,,,] == F
nHaps[2,1,4,,,,] <- nOneM[3,,,,] == T & nOneM[4,,,,] == T

nHaps[2,2,1,,,,] <- nOneP[3,,,,] == F & nOneP[4,,,,] == F
nHaps[2,2,2,,,,] <- nOneP[3,,,,] == F & nOneP[4,,,,] == T
nHaps[2,2,3,,,,] <- nOneP[3,,,,] == T & nOneP[4,,,,] == F
nHaps[2,2,4,,,,] <- nOneP[3,,,,] == T & nOneP[4,,,,] == T


d <- tibble(Generation = numeric(),
            SexRatio = numeric(),
            Females = numeric(),
            Males = numeric(),
            LinkageGroup = character(),
            Origin = character(),
            Haplotype = character(),
            SPY = numeric(),
            SPA = numeric(),
            AY = numeric(),
            AA = numeric(),
            BY = numeric(),
            BA = numeric(),
            RY = numeric(),
            RA = numeric(),
            CY = numeric(),
            CA = numeric(),
            MaleDeterminer = numeric())
for(l in 1:length(LG)) for(o in 1:length(Origin)) for(h in 1:length(Haplotype))
{
  dt <- tibble(Generation = t,
               SexRatio = A %*% IsMale,
               Females = (A * IsFemale) %*% nHaps[l,o,h,,,,],
               Males = (A * IsMale) %*% nHaps[l,o,h,,,,],
               LinkageGroup = LG[l],
               Origin = Origin[o],
               Haplotype = Haplotype[h],
               SPY = s[1],
               SPA = s[2],
               AY = a[1],
               AA = a[2],
               BY = b[1],
               BA = b[2],
               RY = r[1],
               RA = r[2],
               CY = CY,
               CA = CA,

               MaleDeterminer = as.numeric(NovelSexDeterminer)
  )
  d <- d %>% full_join(dt, by = names(dt))
}
write.table(x = d, file = str_c("all_sims/results_sim_", nsim, ".txt"), row.names = F, col.names = T)


d %>% filter(Origin == "Maternal", Haplotype %in% c("10", "11")) %>% 
  select(-Generation, -Males) %>% group_by(LinkageGroup) %>% summarize(sum(Females) / (1 - SexRatio))

dp <- d %>% dplyr::mutate(Females = Females / (1-SexRatio),
                    Males = Males / (SexRatio)) %>% 
  gather(Females, Males, key = "Sex", value = "Frequency") %>% 
  filter(Haplotype %in% c("10", "11")) %>%
  group_by(Sex, Origin, LinkageGroup) %>% 
  summarize(Tot = sum(Frequency)) %>% 
  dplyr::mutate(SLO = str_c(str_sub(Sex, 1,1),
                     LinkageGroup,
                     str_sub(Origin, 1,3))) %>% 
  dplyr::ungroup(Sex, Origin, LinkageGroup) %>% 
  select(SLO, Tot) %>% 
  spread(key = SLO, value = Tot)


dp <- dp %>% dplyr::mutate(SPY = s[1],
                           SPA = s[2],
                           AY = a[1],
                           AA = a[2],
                           BY = b[1],
                           BA = b[2],
                           RY = r[1],
                           RA = r[2],
                           CY = CY,
                           CA = CA,
                           MaleDeterminer = as.numeric(NovelSexDeterminer))
ds <- ds %>% full_join(dp, by = names(dp))
}
write.table(x = ds, file = "results_summarized.txt", col.names = T, row.names = F)


