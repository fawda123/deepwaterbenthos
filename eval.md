
```r
knitr::opts_chunk$set(message = F, warning = F, dev = 'tiff', dev.args = list(tiff = list(compression = 'lzw', family = 'serif')), dpi = 400, out.width = '60%')

library(tidyverse)
library(vegan)
library(ggmap)
library(GGally)
library(ggdendro)
library(dendextend)
library(tibble)
library(d3heatmap)
library(NbClust)
library(ggord)
library(RColorBrewer)
library(gridExtra)
source('R/funcs.R')

data(biogen)
data(env)

# globals
rrs <- 5
evvrs <- c('StationWaterDepth', 'Latitude', 'Clay', 'Sand', 'Silt', 'TN', 'TOC')
ngrps <- 5
```

# Genus-level results

## Exploratory plots


```r
toplo <- table(biogen$nm) %>% 
  sort %>% 
  data.frame %>%
  mutate(Var1 = factor(Var1, levels = as.character(.$Var1)))
ln <- max(toplo$Freq) * rrs / 100
ggplot(toplo, aes(x = Var1, y = Freq)) +
  geom_bar(stat = 'identity') + 
  coord_flip() + 
  theme_bw() +
  theme(
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.y = element_text(size = 5), 
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.y = element_blank()
  ) + 
  geom_hline(yintercept = ln)
```

<img src="eval_files/figure-html/explore-1.tiff" width="60%" />


```r
toplo <- biogen %>% 
  group_by(SampleID) %>% 
  summarise(rich = length(nm)) %>% 
  ungroup %>% 
  left_join(env, by = c('SampleID')) %>% 
  select(-SampleID, -StationID, -Year)

ggpairs(toplo)
```

<img src="eval_files/figure-html/pairs-1.tiff" width="60%" />

## Clustering and ordination


```r
# abundance data in a matrix
dat <- biogen %>% 
  select(SampleID, nm, Abun) %>% 
  group_by(nm) %>% 
  mutate(n = length(nm)) %>% 
  filter(n > ln) %>% 
  ungroup %>% 
  select(-n) %>% 
  spread(nm, Abun, fill = 0) %>% 
  data.frame %>% 
  remove_rownames %>% 
  column_to_rownames('SampleID')

# cluster with euclidean dissim, ward grouping
dend <- dat %>% 
  decostand(method = 'log') %>% 
  vegdist(method = 'bray') %>% 
  hclust(method = 'average')

grps <- cutree(dend, k = ngrps)

p1 <- dend %>% 
  as.dendrogram %>% 
  set("branches_k_color", k = ngrps) %>%
  set("labels_colors", k = ngrps) %>%
  set("labels_cex", 0.4) 
circlize_dendrogram(p1) 
```

<img src="eval_files/figure-html/dendro-1.tiff" width="60%" />


```r
# combine grps 5 and 4 into 3
cols <- gg_color_hue(5)[c(2, 3, 1, 4, 5)]
grps[grps %in% c(4, 5)] <- 3
cols <- cols[1:3]

# ordination
ord <- dat %>% 
  decostand(method = 'log') %>% 
  metaMDS(distance = 'bray', trace = 0, autotransform = F, k = 3, trymax = 200)
ord
```

```
## 
## Call:
## metaMDS(comm = ., distance = "bray", k = 3, trymax = 200, autotransform = F,      trace = 0) 
## 
## global Multidimensional Scaling using monoMDS
## 
## Data:     . 
## Distance: bray 
## 
## Dimensions: 3 
## Stress:     0.1935309 
## Stress type 1, weak ties
## Two convergent solutions found after 20 tries
## Scaling: centring, PC rotation, halfchange scaling 
## Species: expanded scores based on '.'
```

```r
ggord(ord, grp_in = as.character(grps), axes = c("1", "2"), arrow = NULL, txt = NULL, size = 4, 
      obslab = F, alpha = 0.8, cols = cols)
```

<img src="eval_files/figure-html/nms-1.tiff" width="50%" />

```r
# ggord(ord, grp_in = as.character(grps), axes = c("2", "3"), arrow = NULL, txt = NULL, size = 4, 
#       obslab = F, alpha = 0.8, cols = cols)
# 
# ggord(ord, grp_in = as.character(grps), axes = c("1", "3"), arrow = NULL, txt = NULL, size = 4, 
#       obslab = F, alpha = 0.8, cols = cols)
```


```r
bbx <- make_bbox(lon = env$Longitude, lat = env$Latitude, f = 0.1)
bsmap <- get_stamenmap(bbx, maptype = "toner-lite", zoom = 8)

envloc <- env %>% 
  select(StationID, Latitude, Longitude) %>% 
  unique
toplo <- grps %>% 
  data.frame(ngrps = .) %>% 
  rownames_to_column('SampleID') %>% 
  left_join(biogen[, c('SampleID', 'StationID')], ., by = 'SampleID') %>% 
  left_join(envloc, by = 'StationID') %>% 
  mutate(
    Groups = factor(ngrps)
  )

ggmap(bsmap) +   
  geom_point(data = toplo, aes(x = Longitude, y = Latitude, colour = Groups), 
             alpha = 0.6, size = 2) + 
  scale_colour_manual(values = cols) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank()
  )
```

<img src="eval_files/figure-html/spatial-1.tiff" width="60%" />


```r
# abundance data in a matrix
datbio <- biogen %>% 
  select(SampleID, nm, Abun) %>% 
  group_by(nm) %>% 
  mutate(n = length(nm)) %>% 
  filter(n > ln) %>% 
  mutate(
    valt = log(1 + Abun)
  ) %>% 
  ungroup %>% 
  select(-Abun, -n) %>% 
  spread(nm, valt, fill = 0) %>% 
  data.frame %>% 
  arrange(SampleID) %>% 
  remove_rownames %>% 
  column_to_rownames('SampleID')

datenv <- env[, c('SampleID', evvrs)] %>%
  gather('var', 'val', -SampleID) %>% 
  group_by(var) %>% 
  mutate(val = ifelse(is.na(val), mean(val, na.rm = T), val)) %>% 
  ungroup %>% 
  spread(var, val) %>% 
  filter(SampleID %in% row.names(datbio)) %>% 
  arrange(SampleID) %>% 
  data.frame(stringsAsFactors = F) %>% 
  remove_rownames %>% 
  column_to_rownames('SampleID')

ord <- cca(datbio, datenv)

ggord(ord, grp_in = as.character(grps), axes = c('1', '2'), vec_ext = 3, ptslab = T, 
      parse = T, cols = cols)
```

<img src="eval_files/figure-html/cca-1.tiff" width="60%" />

```r
# ggord(ord, grp_in = as.character(grps), axes = c('2', '3'), vec_ext = 3, ptslab = T, 
#       parse = T, cols = cols)
# ggord(ord, grp_in = as.character(grps), axes = c('1', '3'), vec_ext = 3, ptslab = T, 
#       parse = T, cols = cols)
```

## Genera abundance by group


```r
# group data as data frame
grpdat <- grps %>% 
  data.frame(grps = .) %>% 
  rownames_to_column('SampleID')

# colors
colfun <- RColorBrewer::brewer.pal(9, 'Blues') %>% 
  colorRampPalette

# most abundant species, all sites
toplo1 <- biogen %>% 
  left_join(grpdat, by = 'SampleID') %>% 
  dplyr::select(grps, nm, Abun) %>% 
  group_by(grps, nm) %>% 
  summarise(Abun = sum(Abun)) %>% 
  group_by(grps) %>% 
  nest %>% 
  mutate(
    subdat = map(data, function(x){
      
    arrange(x, -Abun) %>% 
        .[1:10, ]
      
    })
  ) %>% 
  dplyr::select(-data) %>% 
  unnest %>% 
  ungroup %>% 
  mutate(
    cols = colfun(length(Abun))[rank(Abun)]
  ) %>% 
  split(., .[, c('grps')])

for(i in seq_along(names(toplo1))){
  
  tmp <- toplo1[[i]] %>% 
    arrange(Abun) %>% 
    mutate(nm = factor(nm, levels = nm))
  
  p <- ggplot(tmp, aes(x = nm, y = Abun, fill = Abun)) + 
    geom_bar(stat = 'identity', colour = 'black', fill = tmp$cols) + 
    ylab('Abundance') + 
    facet_wrap( ~ grps, ncol = 2, scales = 'free') + 
    theme_bw() + 
    scale_fill_gradientn(colours = brewer.pal(9, 'Blues')) +
    theme(legend.position = "none", 
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 8), 
      axis.title.x = element_blank()
      ) + 
    # scale_y_continuous(expand = c(0, 0)) +
    coord_flip()
  
  assign(paste0('p', i), p) 
  
}

# Get the widths
pA <- ggplot_gtable(ggplot_build(p1))
pB <- ggplot_gtable(ggplot_build(p2))
pC <- ggplot_gtable(ggplot_build(p3))
maxWidth = grid::unit.pmax(pA$widths[2:3], pB$widths[2:3], pC$widths[2:3])

# Set the widths
pA$widths[2:3] <- maxWidth
pB$widths[2:3] <- maxWidth
pC$widths[2:3] <- maxWidth

grid.arrange(pA, pB, pC, ncol = 1, bottom = 'Abundance')
```

<img src="eval_files/figure-html/abund-1.tiff" width="60%" />

### Environmental variables by group


```r
# group data as data frame
grpdat <- grps %>% 
  data.frame(grps = .) %>% 
  rownames_to_column('SampleID')
```
