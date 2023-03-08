#####
## DFO and NOAA NMFS species mapping
#####
library(tidyverse)

# DFO species information
load(here::here("data/dfo/raw/RV.GSSPECIES.Rdata"))

names(GSSPECIES)

sp_dat <- GSSPECIES %>%
    dplyr::select(., SPEC, COMM, CODE, NMFS) %>%
    distinct() %>%
    drop_na()

sp_dat %>%
    filter(., grepl("REDFISH", COMM))


# NOAA species information
noaa_spp <- read_csv(here::here("data/supporting/species_table.csv")) %>%
    dplyr::select(., NMFS_SVSPP, Group) %>%
    left_join(., sp_dat, by = c("NMFS_SVSPP" = "NMFS")) %>%
    dplyr::select(., NMFS_SVSPP, CODE, COMM, SPEC, Group)
names(noaa_spp)[2]<- "DFO_SPEC"

write.csv(noaa_spp, here::here("data/supporting/species_table.csv"))
