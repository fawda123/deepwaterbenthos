library(tidyverse)
library(lubridate)

# species level data
biospp <- read.csv('ignore/clean benthos v3.csv', stringsAsFactors = FALSE) %>% 
  filter(Month %in% c(7, 8, 9)) %>% 
  unite('SampleID', StationID, Year, sep = '_', remove = F) %>% 
  group_by(SampleID) %>% 
  mutate(n = length(Taxon)) %>% 
  ungroup %>% 
  filter(n > 5) %>% 
  rename(nm = Taxon) %>% 
  select(SampleID, StationID, Year, nm, Abun)

# genus level data
biogen <- read.csv('ignore/Clean benthos v3 - Genus2.csv', stringsAsFactors = FALSE) %>% 
  filter(Month %in% c(7, 8, 9)) %>% 
  unite('SampleID', StationID, Year, sep = '_', remove = F) %>% 
  group_by(SampleID) %>% 
  mutate(n = length(Genus)) %>% 
  ungroup %>% 
  filter(n > 5) %>% 
  rename(
    nm = Genus
    ) %>% 
  select(SampleID, StationID, Year, nm, G_Abun) %>% 
  group_by(SampleID, StationID, Year, nm) %>%
  summarise(Abun = round(mean(G_Abun), 0)) %>% 
  ungroup

# family level data
biofam <- read.csv('ignore/Clean benthos v3 - family2.csv', stringsAsFactors = FALSE) %>% 
  filter(Month %in% c(7, 8, 9)) %>% 
  unite('SampleID', StationID, Year, sep = '_', remove = F) %>% 
  group_by(SampleID) %>% 
  mutate(n = length(Family)) %>% 
  ungroup %>% 
  filter(n > 5) %>% 
  rename(
    nm = Family,
    Abun = F_Abun
    ) %>% 
  select(SampleID, StationID, Year, nm, Abun)

# environmental data
env <- read.csv('ignore/clean stations v2 gradients.csv', stringsAsFactors = FALSE) %>% 
  select(-Source, -NFail1, -NFail2, -ERM_N, -PlatformDist_km) %>% 
  mutate(Longitude = -1 * abs(Longitude)) %>% 
  unite('SampleID', StationID, Year, sep = '_', remove = F)

# save all
save(biospp, file = 'data/biospp.RData')
save(biogen, file = 'data/biogen.RData')
save(biofam, file = 'data/biofam.RData')
save(env, file = 'data/env.RData')

