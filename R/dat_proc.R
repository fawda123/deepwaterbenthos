library(tidyverse)
library(lubridate)

bio <- read.csv('ignore/clean benthos v3.csv', stringsAsFactors = FALSE) %>% 
  filter(Month %in% c(7, 8, 9)) %>% 
  unite('SampleID', StationID, Year, sep = '_', remove = F) %>% 
  group_by(SampleID) %>% 
  mutate(n = length(Taxon)) %>% 
  ungroup %>% 
  filter(n > 5) %>% 
  select(SampleID, StationID, Year, Taxon, Abun)

env <- read.csv('ignore/clean stations v2 gradients.csv', stringsAsFactors = FALSE) %>% 
  select(-Source, -NFail1, -NFail2, -ERM_N, -PlatformDist_km) %>% 
  mutate(Longitude = -1 * abs(Longitude)) %>% 
  unite('SampleID', StationID, Year, sep = '_', remove = F)

save(bio, file = 'data/bio.RData')
save(env, file = 'data/env.RData')

