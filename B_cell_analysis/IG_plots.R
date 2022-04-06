#!/usr/bin/env Rscript

### title: Print out heatmaps of IG combination usage in B cells 
### authors: Yiping Wang, Jana Biermann
library(Seurat)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)

seu <- readRDS("data/MBPM/data_MBPM_scn.rds")
seu <- subset(seu, cell_type_int %in%c('Plasma cells','B cells'))
uniqueIDs <- unique(seu$patient)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7")

# read in data into IGtable variable, filtering out rows with NAs for both heavy
# and light chains
IGtable <- read.csv("data/MBPM/IG_analysis/table_IG_MBPM_bcells_chains.csv")
IGtable <- IGtable[!is.na(IGtable$light) & !is.na(IGtable$heavy), ]

# read in order of heavy and light chains from IG combination frequency heatmap
uniquelight = rev(read.table("data/MBPM/IG_analysis/lightorder.txt", header = F, quote = "")$V1)
uniqueheavy = read.table("data/MBPM/IG_analysis/heavyorder.txt", header = F, quote = "")$V1

# store indices of where IGtable light and heavy chains match uniquelight and
# uniqueheavy
IGtable$lightidxs = match(IGtable$light, uniquelight)
IGtable$heavyidxs = match(IGtable$heavy, uniqueheavy)

# make list of all light and heavy chain combinations that occur in IGtable,
# store in combs make list of combinations that are occur more than once, store
# in dupcombs
combs = c()
dupcombs = c()
for (i in 1:length(IGtable$light)) {
  if (sum(combs == paste0(IGtable$light[i], " ", IGtable$heavy[i])) != 0 & 
      sum(dupcombs == paste0(IGtable$light[i], " ", IGtable$heavy[i])) == 0) {
    dupcombs = append(dupcombs, paste0(IGtable$light[i], " ", IGtable$heavy[i]))
  }
  combs = append(combs, paste0(IGtable$light[i], " ", IGtable$heavy[i]))
}


#### constants table
# for duplicated combinations, it is possible that each entry may have different
# values in the constant field find all unique values in constant field for each
# duplicated combination in the cell for each duplicated combination, subdivide
# the cell along the light chain axis, so that each subcell will represent a
# unique constant value store the endpoints of each subcell in lightidxs and
# lightidxsceil fields
duptable = data.frame(light = character(), heavy = character(), constant = character(), 
                      lightidxs = integer(), lightidxsceil = integer(), heavyidxs = integer())
for (i in 1:length(dupcombs)) {
  words = strsplit(dupcombs[i], " ")[[1]]
  smalltable = na.omit(IGtable[IGtable$light == words[1] & IGtable$heavy == words[2], ])
  smalltableconstant = unique(smalltable$constant)
  print(length(smalltableconstant))
  if (length(smalltableconstant) > 1) {
    for (j in 1:length(smalltableconstant)) {
      duptable = rbind(duptable, 
                       data.frame(light = words[1], 
                                  heavy = words[2], 
                                  constant = smalltableconstant[j], 
                                  lightidxs = match(words[1], uniquelight) + (j - 1)/length(smalltableconstant), 
                                  lightidxsceil = match(words[1], uniquelight) + j/length(smalltableconstant), 
                                  heavyidxs = match(words[2], uniqueheavy)))
    }
  }
}

# create heatmap of all light ahd heavy chain combinations, with subcells for
# multiple constant values appearing in one combination
p = ggplot(IGtable, aes(x = heavyidxs + 0.5, y = lightidxs + 0.5, fill = constant, width = 1, height = 1)) + 
  geom_tile()

if (dim(duptable)[1] != 0) {
  p = p + geom_rect(data = duptable, aes(xmin = heavyidxs, xmax = heavyidxs + 1, ymin = lightidxs, 
                                         ymax = lightidxsceil, fill = constant), 
                    colour = "black", size = 0.3)
}

p = p + scale_x_continuous("Heavy Chain", breaks = 1:(length(uniqueheavy) + 1), labels = append(uniqueheavy, "")) + 
  scale_y_continuous("Light Chain", breaks = 1:(length(uniquelight) + 1), labels = append(uniquelight, "")) + 
  scale_fill_manual(breaks = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHGP", "IGHM"), 
                    values = cbPalette[1:8]) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, margin = margin(r = 0)), 
        axis.text.y = element_text(vjust = 0, margin = margin(r = 0)), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.grid = element_line(size = 0.1, colour = "darkgrey"), 
        panel.ontop = T)

ggsave("data/MBPM/IG_analysis/plots_IG_combination_with_constant.pdf", height = 12)


#### occurrence in patients
# count number of times each IG chain combination occurs in IGtable, store in
# sharedIGtable, along with lightidxs and heavyidxs
sharedIGtable = data.frame(heavy = character(), light = character(), occurrences = integer(), 
                           heavyidxs = integer(), lightidxs = integer())
for (i in 1:length(uniqueIDs)) {
  IGtablesmall = IGtable[IGtable$barcode %in% 
                           names(seu$patient[seu$patient == uniqueIDs[i]]), ]
  if (dim(IGtablesmall)[1] != 0) {
    checktable = data.frame(heavy = character(), light = character())
    for (j in 1:dim(IGtablesmall)[1]) {
      if (!is.na(IGtablesmall$heavy[j]) & !is.na(IGtablesmall$light[j])) {
        if (sum(checktable$heavy == IGtablesmall$heavy[j] & checktable$light == 
                IGtablesmall$light[j]) == 0) {
          smalltable = data.frame(heavy = IGtablesmall$heavy[j], light = IGtablesmall$light[j])
          checktable = rbind(checktable, smalltable)
          
          matchingarr = (sharedIGtable$heavy == IGtablesmall$heavy[j] & sharedIGtable$light == IGtablesmall$light[j])
          if (sum(matchingarr) == 0) {
            smalltable = data.frame(heavy = IGtablesmall$heavy[j], 
                                    light = IGtablesmall$light[j], 
                                    occurrences = 1, 
                                    heavyidxs = match(IGtablesmall$heavy[j], uniqueheavy), 
                                    lightidxs = match(IGtablesmall$light[j], uniquelight))
            sharedIGtable = rbind(sharedIGtable, smalltable)
          } else {
            sharedIGtable$occurrences[matchingarr] = sharedIGtable$occurrences[matchingarr] + 1
          }
        }
      }
    }
  }
}

# print heatmap of number of times each IG chain combination occurs in IGtable
sharedIGtable$occurrences = as.character(sharedIGtable$occurrences)
ggplot(sharedIGtable, aes(x = heavyidxs + 0.5, y = lightidxs + 0.5, fill = occurrences, width = 1, height = 1)) + 
  geom_tile() + 
  scale_x_continuous("Heavy Chain", breaks = 1:(length(uniqueheavy) + 1), labels = append(uniqueheavy, "")) + 
  scale_y_continuous("Light Chain", breaks = 1:(length(uniquelight) + 1), labels = append(uniquelight, "")) + 
  scale_fill_manual(breaks = c("1", "2", "3", "4", "5"), values = cbPalette[1:5]) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, margin = margin(r = 0)), 
        axis.text.y = element_text(vjust = 0, margin = margin(r = 0)), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.grid = element_line(size = 0.1, colour = "darkgrey"), 
        panel.ontop = T)
ggsave("data/MBPM/IG_analysis/plots_IG_combination_occurence_in_patient.pdf")


#### individual IG combis
# create heatmaps of IG combinations that occur in each individual sample
for(i in 1:length(uniqueIDs)) {
  IGtablesmall = IGtable[IGtable$barcode %in% names(seu$patient[seu$patient == uniqueIDs[i]]), ]
  
  if(dim(IGtablesmall)[1] != 0) {
    combssmall = c()
    dupcombssmall = c()
    duptablesmall = data.frame(light = character(), 
                               heavy = character(), 
                               constant = character(), 
                               lightidxs = integer(), 
                               lightidxsceil = integer(), 
                               heavyidxs = integer())
    
    for(j in 1:length(IGtablesmall$light)) {
      if(sum(combssmall == paste0(IGtablesmall$light[j], " ", IGtablesmall$heavy[j])) != 0 & 
         sum(dupcombssmall == paste0(IGtablesmall$light[j], " ", IGtablesmall$heavy[j])) == 0) {
        dupcombssmall = append(dupcombssmall, 
                               paste0(IGtablesmall$light[j], " ", IGtablesmall$heavy[j]))
      }
      
      combssmall = append(combssmall, paste0(IGtablesmall$light[j], " ", IGtablesmall$heavy[j]))
    }
    
    if(length(dupcombssmall) != 0) {
      for (k in 1:length(dupcombssmall)) {
        words = strsplit(dupcombssmall[k], " ")[[1]]
        smalltable = na.omit(IGtablesmall[IGtablesmall$light == words[1] & 
                                            IGtablesmall$heavy == words[2], ])
        smalltableconstant = unique(smalltable$constant)
        if (length(smalltableconstant) > 1) {
          for (l in 1:length(smalltableconstant)) {
            duptablesmall = rbind(duptablesmall, 
                                  data.frame(light = words[1], 
                                             heavy = words[2], 
                                             constant = smalltableconstant[l], 
                                             lightidxs = match(words[1], uniquelight) + (l - 1)/length(smalltableconstant), 
                                             lightidxsceil = match(words[1], uniquelight) + l/length(smalltableconstant), 
                                             heavyidxs = match(words[2], uniqueheavy)))
          }
        }
      }
    }
    
    p = ggplot(IGtablesmall, aes(x = heavyidxs + 0.5, y = lightidxs + 0.5, fill = constant, width = 1, height = 1)) + 
      geom_tile()
    
    if(dim(duptablesmall)[1] > 0) {
      p = p + geom_rect(data = duptablesmall, 
                        aes(xmin = heavyidxs, xmax = heavyidxs + 1, ymin = lightidxs, ymax = lightidxsceil, 
                            fill = constant), colour = "black", size = 0.3)
    }
    
    p = p + scale_x_continuous("Heavy Chain", breaks = 1:(length(uniqueheavy) + 1), 
                               labels = append(uniqueheavy, ""), limits = c(1, length(uniqueheavy))) + 
      scale_y_continuous("Light Chain", breaks = 1:(length(uniquelight) + 1), labels = append(uniquelight, ""), 
                         limits = c(1, length(uniquelight))) + 
      scale_fill_manual(breaks = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHGP", "IGHM"), 
                        values = cbPalette[1:8]) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, margin = margin(r = 0)), 
            axis.text.y = element_text(vjust = 0, margin = margin(r = 0)), 
            panel.grid.minor = element_blank(), panel.background = element_blank(), 
            panel.grid = element_line(size = 0.1, colour = "darkgrey"), panel.ontop = T)  # + xlim(0,100) + ylim(0,100)
    ggsave(paste0("data/MBPM/IG_analysis/IG_combination_", uniqueIDs[i], ".pdf"), height = 10)
    print(uniqueIDs[i])
  }
}


### organ tile plot
duptable2 = data.frame(light = character(), heavy = character(), organ = character(), 
                      lightidxs = integer(), lightidxsceil = integer(), heavyidxs = integer())
for (i in 1:length(dupcombs)) {
  words = strsplit(dupcombs[i], " ")[[1]]
  smalltable = na.omit(IGtable[IGtable$light == words[1] & IGtable$heavy == words[2], ])
  smalltableorgan = unique(smalltable$organ)
  print(length(smalltableorgan))
  if (length(smalltableorgan) > 1) {
    for (j in 1:length(smalltableorgan)) {
      duptable2 = rbind(duptable2, 
                       data.frame(light = words[1], 
                                  heavy = words[2], 
                                  organ = smalltableorgan[j], 
                                  lightidxs = match(words[1], uniquelight) + (j - 1)/length(smalltableorgan), 
                                  lightidxsceil = match(words[1], uniquelight) + j/length(smalltableorgan), 
                                  heavyidxs = match(words[2], uniqueheavy)))
    }
  }
}

# create heatmap of all light ahd heavy chain combinations, with subcells for
# multiple organ values appearing in one combination
p = ggplot(IGtable, aes(x = heavyidxs + 0.5, y = lightidxs + 0.5, fill = organ, width = 1, height = 1)) + 
  geom_tile()

if (dim(duptable2)[1] != 0) {
  p = p + geom_rect(data = duptable2, aes(xmin = heavyidxs, xmax = heavyidxs + 1, ymin = lightidxs, 
                                         ymax = lightidxsceil, fill = organ), colour = "black", size = 0.1)
}

p = p + scale_x_continuous("Heavy Chain", breaks = 1:(length(uniqueheavy) + 1), labels = append(uniqueheavy, "")) + 
  scale_y_continuous("Light Chain", breaks = 1:(length(uniquelight) + 1), labels = append(uniquelight, "")) + 
  #scale_fill_manual(breaks = c("Brain",'Peripheral'), values = cbPalette[1:2]) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, margin = margin(r = 0)), 
        axis.text.y = element_text(vjust = 0, margin = margin(r = 0)), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.grid = element_line(size = 0.1, colour = "darkgrey"), 
        panel.ontop = T)
ggsave("data/MBPM/IG_analysis/plots_IG_combination_with_organ.pdf")
