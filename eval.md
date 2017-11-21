
```r
knitr::opts_chunk$set(message = F, warning = F)

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

data(biospp)
data(biogen)
data(biofam)
data(env)

# globals
rrs <- 5
evvrs <- c('StationWaterDepth', 'Latitude', 'Clay', 'Sand', 'Silt', 'TN', 'TOC')
ngrps <- 5
```

# {.tabset}

## Species 

### Exploratory plots


```r
toplo <- table(biospp$nm) %>% 
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

![](eval_files/figure-html/unnamed-chunk-2-1.png)<!-- -->


```r
bbx <- make_bbox(lon = env$Longitude, lat = env$Latitude, f = 0.2)
bsmap <- get_stamenmap(bbx, maptype = "toner-lite", zoom = 9)

envloc <- env %>% 
  select(StationID, Latitude, Longitude, TOC) %>% 
  group_by(StationID, Latitude, Longitude) %>% 
  summarise(TOC = mean(TOC, na.rm = T)) %>% 
  unique
toplo <- biospp %>% 
  group_by(SampleID, StationID) %>% 
  summarise(rich = length(nm)) %>% 
  group_by(StationID) %>% 
  summarise(rich = mean(rich)) %>% 
  ungroup %>% 
  mutate(rich = round(rich, 0)) %>% 
  left_join(envloc, by = 'StationID')

ggmap(bsmap) +   
  geom_point(data = toplo, aes(x = Longitude, y = Latitude, size = rich, colour = TOC), alpha = 0.6) + 
  theme_bw()
```

![](eval_files/figure-html/unnamed-chunk-3-1.png)<!-- -->


```r
toplo <- biospp %>% 
  group_by(SampleID) %>% 
  summarise(rich = length(nm)) %>% 
  ungroup %>% 
  left_join(env, by = c('SampleID')) %>% 
  select(-SampleID, -StationID, -Year)

ggpairs(toplo)
```

![](eval_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

### Clustering and ordination


```r
# abundance data in a matrix
dat <- biospp %>% 
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

![](eval_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


```r
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
## Stress:     0.189162 
## Stress type 1, weak ties
## Two convergent solutions found after 53 tries
## Scaling: centring, PC rotation, halfchange scaling 
## Species: expanded scores based on '.'
```

```r
ggord(ord, grp_in = as.character(grps), axes = c("1", "2"), arrow = NULL, txt = NULL, size = 4, 
      obslab = F, ellipse = F, hull = T, alpha = 0.8)
```

![](eval_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
ggord(ord, grp_in = as.character(grps), axes = c("2", "3"), arrow = NULL, txt = NULL, size = 4, 
      obslab = F, ellipse = F, hull = T, alpha = 0.8)
```

![](eval_files/figure-html/unnamed-chunk-6-2.png)<!-- -->

```r
ggord(ord, grp_in = as.character(grps), axes = c("1", "3"), arrow = NULL, txt = NULL, size = 4, 
      obslab = F, ellipse = F, hull = T, alpha = 0.8)
```

![](eval_files/figure-html/unnamed-chunk-6-3.png)<!-- -->


```r
bbx <- make_bbox(lon = env$Longitude, lat = env$Latitude, f = 0.2)
bsmap <- get_stamenmap(bbx, maptype = "toner-lite", zoom = 8)

envloc <- env %>% 
  select(StationID, Latitude, Longitude) %>% 
  unique
toplo <- grps %>% 
  data.frame(ngrps = .) %>% 
  rownames_to_column('SampleID') %>% 
  left_join(biospp[, c('SampleID', 'StationID')], ., by = 'SampleID') %>% 
  left_join(envloc, by = 'StationID')

ggmap(bsmap) +   
  geom_point(data = toplo, aes(x = Longitude, y = Latitude, colour = factor(ngrps)), 
             alpha = 0.6, size = 2) + 
  theme_bw()
```

![](eval_files/figure-html/unnamed-chunk-7-1.png)<!-- -->


```r
# abundance data in a matrix
datbio <- biospp %>% 
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
      parse = T, hull = T, ellipse = F)
```

![](eval_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
ggord(ord, grp_in = as.character(grps), axes = c('2', '3'), vec_ext = 3, ptslab = T, 
      parse = T, hull = T, ellipse = F)
```

![](eval_files/figure-html/unnamed-chunk-8-2.png)<!-- -->

```r
ggord(ord, grp_in = as.character(grps), axes = c('1', '3'), vec_ext = 3, ptslab = T, 
      parse = T, hull = T, ellipse = F)
```

![](eval_files/figure-html/unnamed-chunk-8-3.png)<!-- -->

## Genus 

### Exploratory plots


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

![](eval_files/figure-html/unnamed-chunk-9-1.png)<!-- -->


```r
bbx <- make_bbox(lon = env$Longitude, lat = env$Latitude, f = 0.2)
bsmap <- get_stamenmap(bbx, maptype = "toner-lite", zoom = 9)

envloc <- env %>% 
  select(StationID, Latitude, Longitude, TOC) %>% 
  group_by(StationID, Latitude, Longitude) %>% 
  summarise(TOC = mean(TOC, na.rm = T)) %>% 
  unique
toplo <- biogen %>% 
  group_by(SampleID, StationID) %>% 
  summarise(rich = length(nm)) %>% 
  group_by(StationID) %>% 
  summarise(rich = mean(rich)) %>% 
  ungroup %>% 
  mutate(rich = round(rich, 0)) %>% 
  left_join(envloc, by = 'StationID')

ggmap(bsmap) +   
  geom_point(data = toplo, aes(x = Longitude, y = Latitude, size = rich, colour = TOC), alpha = 0.6) + 
  theme_bw()
```

![](eval_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


```r
toplo <- biogen %>% 
  group_by(SampleID) %>% 
  summarise(rich = length(nm)) %>% 
  ungroup %>% 
  left_join(env, by = c('SampleID')) %>% 
  select(-SampleID, -StationID, -Year)

ggpairs(toplo)
```

![](eval_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

### Clustering and ordination


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

![](eval_files/figure-html/unnamed-chunk-12-1.png)<!-- -->


```r
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
## Stress:     0.193531 
## Stress type 1, weak ties
## Two convergent solutions found after 20 tries
## Scaling: centring, PC rotation, halfchange scaling 
## Species: expanded scores based on '.'
```

```r
ggord(ord, grp_in = as.character(grps), axes = c("1", "2"), arrow = NULL, txt = NULL, size = 4, 
      obslab = F, ellipse = F, hull = T, alpha = 0.8)
```

![](eval_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

```r
ggord(ord, grp_in = as.character(grps), axes = c("2", "3"), arrow = NULL, txt = NULL, size = 4, 
      obslab = F, ellipse = F, hull = T, alpha = 0.8)
```

![](eval_files/figure-html/unnamed-chunk-13-2.png)<!-- -->

```r
ggord(ord, grp_in = as.character(grps), axes = c("1", "3"), arrow = NULL, txt = NULL, size = 4, 
      obslab = F, ellipse = F, hull = T, alpha = 0.8)
```

![](eval_files/figure-html/unnamed-chunk-13-3.png)<!-- -->


```r
bbx <- make_bbox(lon = env$Longitude, lat = env$Latitude, f = 0.2)
bsmap <- get_stamenmap(bbx, maptype = "toner-lite", zoom = 8)

envloc <- env %>% 
  select(StationID, Latitude, Longitude) %>% 
  unique
toplo <- grps %>% 
  data.frame(ngrps = .) %>% 
  rownames_to_column('SampleID') %>% 
  left_join(biogen[, c('SampleID', 'StationID')], ., by = 'SampleID') %>% 
  left_join(envloc, by = 'StationID')

ggmap(bsmap) +   
  geom_point(data = toplo, aes(x = Longitude, y = Latitude, colour = factor(ngrps)), 
             alpha = 0.6, size = 2) + 
  theme_bw()
```

![](eval_files/figure-html/unnamed-chunk-14-1.png)<!-- -->


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
      parse = T, hull = T, ellipse = F)
```

![](eval_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

```r
ggord(ord, grp_in = as.character(grps), axes = c('2', '3'), vec_ext = 3, ptslab = T, 
      parse = T, hull = T, ellipse = F)
```

![](eval_files/figure-html/unnamed-chunk-15-2.png)<!-- -->

```r
ggord(ord, grp_in = as.character(grps), axes = c('1', '3'), vec_ext = 3, ptslab = T, 
      parse = T, hull = T, ellipse = F)
```

![](eval_files/figure-html/unnamed-chunk-15-3.png)<!-- -->

## Family

### Exploratory plots


```r
toplo <- table(biofam$nm) %>% 
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

![](eval_files/figure-html/unnamed-chunk-16-1.png)<!-- -->


```r
bbx <- make_bbox(lon = env$Longitude, lat = env$Latitude, f = 0.2)
bsmap <- get_stamenmap(bbx, maptype = "toner-lite", zoom = 9)

envloc <- env %>% 
  select(StationID, Latitude, Longitude, TOC) %>% 
  group_by(StationID, Latitude, Longitude) %>% 
  summarise(TOC = mean(TOC, na.rm = T)) %>% 
  unique
toplo <- biofam %>% 
  group_by(SampleID, StationID) %>% 
  summarise(rich = length(nm)) %>% 
  group_by(StationID) %>% 
  summarise(rich = mean(rich)) %>% 
  ungroup %>% 
  mutate(rich = round(rich, 0)) %>% 
  left_join(envloc, by = 'StationID')

ggmap(bsmap) +   
  geom_point(data = toplo, aes(x = Longitude, y = Latitude, size = rich, colour = TOC), alpha = 0.6) + 
  theme_bw()
```

![](eval_files/figure-html/unnamed-chunk-17-1.png)<!-- -->


```r
toplo <- biofam %>% 
  group_by(SampleID) %>% 
  summarise(rich = length(nm)) %>% 
  ungroup %>% 
  left_join(env, by = c('SampleID')) %>% 
  select(-SampleID, -StationID, -Year)

ggpairs(toplo)
```

![](eval_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

### Clustering and ordination


```r
# abundance data in a matrix
dat <- biofam %>% 
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

![](eval_files/figure-html/unnamed-chunk-19-1.png)<!-- -->


```r
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
## Stress:     0.1991924 
## Stress type 1, weak ties
## Two convergent solutions found after 75 tries
## Scaling: centring, PC rotation, halfchange scaling 
## Species: expanded scores based on '.'
```

```r
ggord(ord, grp_in = as.character(grps), axes = c("1", "2"), arrow = NULL, txt = NULL, size = 4, 
      obslab = F, ellipse = F, hull = T, alpha = 0.8)
```

![](eval_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

```r
ggord(ord, grp_in = as.character(grps), axes = c("2", "3"), arrow = NULL, txt = NULL, size = 4, 
      obslab = F, ellipse = F, hull = T, alpha = 0.8)
```

![](eval_files/figure-html/unnamed-chunk-20-2.png)<!-- -->

```r
ggord(ord, grp_in = as.character(grps), axes = c("1", "3"), arrow = NULL, txt = NULL, size = 4, 
      obslab = F, ellipse = F, hull = T, alpha = 0.8)
```

![](eval_files/figure-html/unnamed-chunk-20-3.png)<!-- -->


```r
bbx <- make_bbox(lon = env$Longitude, lat = env$Latitude, f = 0.2)
bsmap <- get_stamenmap(bbx, maptype = "toner-lite", zoom = 8)

envloc <- env %>% 
  select(StationID, Latitude, Longitude) %>% 
  unique
toplo <- grps %>% 
  data.frame(ngrps = .) %>% 
  rownames_to_column('SampleID') %>% 
  left_join(biofam[, c('SampleID', 'StationID')], ., by = 'SampleID') %>% 
  left_join(envloc, by = 'StationID')

ggmap(bsmap) +   
  geom_point(data = toplo, aes(x = Longitude, y = Latitude, colour = factor(ngrps)), 
             alpha = 0.6, size = 2) + 
  theme_bw()
```

![](eval_files/figure-html/unnamed-chunk-21-1.png)<!-- -->


```r
# abundance data in a matrix
datbio <- biofam %>% 
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
      parse = T, hull = T, ellipse = F)
```

![](eval_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

```r
ggord(ord, grp_in = as.character(grps), axes = c('2', '3'), vec_ext = 3, ptslab = T, 
      parse = T, hull = T, ellipse = F)
```

![](eval_files/figure-html/unnamed-chunk-22-2.png)<!-- -->

```r
ggord(ord, grp_in = as.character(grps), axes = c('1', '3'), vec_ext = 3, ptslab = T, 
      parse = T, hull = T, ellipse = F)
```

![](eval_files/figure-html/unnamed-chunk-22-3.png)<!-- -->
