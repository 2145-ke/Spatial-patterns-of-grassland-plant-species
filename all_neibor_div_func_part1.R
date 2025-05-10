# load packages ----
library(geoR)
library(spatstat)
library(codyn)
library(vegan)
library(betapart)
library(readxl)
library(dplyr)
library(tidyr)
library(writexl)
library(deldir)
library(psych)
# load data 2023 ----
setwd("D:/OneDrive/Rdata/7 spatial patterns of grassland plant species")
dot.comp_all <- read_excel("all_points_across_six_sites_2023.xlsx",sheet=1,na="NA")
head(dot.comp_all)
headTail(dot.comp_all)
str(dot.comp_all)
summary(dot.comp_all)

dot.comp_all$X_mm <- as.numeric(dot.comp_all$X_mm)
dot.comp_all$Y_mm <- as.numeric(dot.comp_all$Y_mm)

# 1. Spatial distribution pattern of grassland plant populations ----
## old species segregation index - relative neighborhood density Ωr ----
Ω_all_site_0_500r <- data.frame()

for (i in unique(dot.comp_all$site)) {
  for(j in unique(dot.comp_all$plot)) {
    subset_data <- dot.comp_all[dot.comp_all$site==i&dot.comp_all$plot==j,]
    
    subset_data <- subset_data %>%
      dplyr::select(!year) %>%
      filter(!duplicated(cbind(X_mm, Y_mm)))
    
    subset_data$species_English <- as.factor(subset_data$species_English) # make species code a factor
    
    Ω0_500r <- data.frame()
    bands <- seq(0, 1900, 100)
    
    for (k in 1:nrow(subset_data)) {
      for (m in bands) {
        
        spe_target <- subset_data$species_English[k]
        
        sub_Ω0_500r <- subset_data[-k,] %>%
          mutate(x_center = subset_data$X_mm[k],
                 y_center = subset_data$Y_mm[k]) %>%
          mutate(distance = sqrt((X_mm - x_center)^2 + (Y_mm - y_center)^2)) %>%
          filter(distance > m & distance <= m + 100) %>%
          mutate(spe_identity = if_else(species_English %in% spe_target, 'cD', 'hD')) %>%
          group_by(spe_identity) %>%
          dplyr::summarise(D_value = sum(abun)) %>%
          pivot_wider(names_from = spe_identity, values_from = D_value) %>%
          mutate(species_English = spe_target,
                 x_center = subset_data$X_mm[k],
                 y_center = subset_data$Y_mm[k],
                 cover_cm2 = subset_data$cover_cm2[k],
                 high_mm = subset_data$high_mm[k],
                 band = m,
                 area = case_when(
                   x_center >= m+100&x_center <= 2000-(m+100)&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2, 
                   x_center < m+100&x_center >= m&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2),
                   y_center < m+100&y_center >= m&x_center >= m+100&x_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center >= m+100&x_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2),
                   x_center < m&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   x_center > 2000-m&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   y_center < m&x_center >= m+100&x_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2),
                   y_center > 2000-m&x_center >= m+100&x_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   x_center < m+100&x_center >= m&y_center < m+100&y_center >= m&sqrt(x_center^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2),
                   x_center < m+100&x_center >= m&y_center < m+100&y_center >= m&sqrt(x_center^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 - (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 + ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2,
                   x_center < m+100&x_center >= m&y_center > 2000-(m+100)&y_center <= 2000-m&sqrt(x_center^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2),
                   x_center < m+100&x_center >= m&y_center > 2000-(m+100)&y_center <= 2000-m&sqrt(x_center^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 -(sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 + ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2,
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center < m+100&y_center >= m&sqrt((2000-x_center)^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center < m+100&y_center >= m&sqrt((2000-x_center)^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 - (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 + ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))- atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2,
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center > 2000-(m+100)&y_center <= 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center > 2000-(m+100)&y_center <= 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 - (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 + ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2,
                   x_center < m+100&x_center >= m&y_center < m&sqrt(x_center^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2),
                   x_center < m+100&x_center >= m&y_center < m&sqrt(x_center^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + + (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 + (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 - ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2 + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2), 
                   y_center < m+100&y_center >= m&x_center < m&sqrt(x_center^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   y_center < m+100&y_center >= m&x_center < m&sqrt(x_center^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + + (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 + (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 - ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2 + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   x_center < m+100&x_center >= m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   x_center < m+100&x_center >= m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + + (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 + (sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 - ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2 + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   y_center < m+100&y_center >= m&x_center > 2000-m&sqrt((2000-x_center)^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   y_center < m+100&y_center >= m&x_center > 2000-m&sqrt((2000-x_center)^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + + (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 + (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 - ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))-atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + + (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 + (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 - ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))-atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center < m&sqrt(x_center^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center < m&sqrt(x_center^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + + (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 + (sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 - ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2 + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + + (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 + (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 - ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + + (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 + (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 - ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   x_center < m&y_center < m&sqrt(x_center^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2), 
                   x_center < m&y_center < m&sqrt(x_center^2 + y_center^2) < m + 100&sqrt(x_center^2 + y_center^2) >= m ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 - (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 + ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2, 
                   x_center < m&y_center < m&sqrt(x_center^2 + y_center^2) < m ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 - (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 + ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2 + (sqrt(m^2-y_center^2)-x_center)*y_center/2 + (sqrt(m^2-x_center^2)-y_center)*x_center/2 - ((acos(y_center/m)-atan(x_center/y_center)) + (acos(x_center/m)-atan(y_center/x_center)))/(2*pi)*pi*m^2, 
                   x_center > 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2), 
                   x_center > 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) < m + 100&sqrt((2000-x_center)^2 + y_center^2) >= m ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 - (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 + ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))-atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2, 
                   x_center > 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) < m ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 - (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 + ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))-atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + (sqrt(m^2-y_center^2)-(2000-x_center))*y_center/2 + (sqrt(m^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 - ((acos(y_center/m)-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/m)-atan(y_center/(2000-x_center))))/(2*pi)*pi*m^2, 
                   x_center < m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2), 
                   x_center < m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) < m + 100&sqrt(x_center^2 + (2000-y_center)^2) >= m ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 - (sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 + ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2, 
                   x_center < m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) < m ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 - (sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 + ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2 + (sqrt(m^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 + (sqrt(m^2-x_center^2)-(2000-y_center))*x_center/2 - ((acos((2000-y_center)/m)-atan(x_center/(2000-y_center))) + (acos(x_center/m)-atan((2000-y_center)/x_center)))/(2*pi)*pi*m^2, 
                   x_center > 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2), 
                   x_center > 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m + 100&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 - (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 + ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2, 
                   x_center > 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 - (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 + ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + (sqrt(m^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 + (sqrt(m^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 - ((acos((2000-y_center)/m)-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/m)-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*m^2
                 ))
                                                                                                                                                                                                                                                                                                                                        
        Ω0_500r <- bind_rows(Ω0_500r, sub_Ω0_500r)
        
      }
      
    }    
    sel_species <- subset_data %>% 
      dplyr::select(species_English) %>%
      group_by(species_English) %>% 
      mutate(indi_num = length(species_English)) %>%
      filter(indi_num >= 20) %>%
      filter(!duplicated(species_English)) %>%
      mutate(mean_cden = indi_num/(2000*2000),
             mean_hden = (nrow(subset_data)-indi_num)/(2000*2000))
    
    Ω0_500r <- Ω0_500r %>% 
      filter(species_English %in% sel_species$species_English) %>% 
      merge(sel_species) %>%
      group_by(species_English, band) %>%
      dplyr::summarise(mean_cden = mean(mean_cden),
                       mean_hden = mean(mean_hden),
                       tcD = sum(cD),
                       thD = sum(hD),
                       tarea = sum(area),
                       spe_abun = mean(indi_num)) %>%
      mutate(rtcD = tcD/tarea,
             rthD = thD/tarea) %>%
      mutate(cΩ = rtcD/mean_cden,
             hΩ = rthD/mean_hden,
             #anpp_target = cover_cm2*high_mm, 
             site = i, plot = j)
    
    Ω_all_site_0_500r <- bind_rows(Ω_all_site_0_500r, Ω0_500r)
    
  }
  
}

write_xlsx(Ω_all_site_0_500r, 'Ω_all_site_0_500r.xlsx')

## old bootstrap method estimating confidence limits of Ωr ----
Ω_all_site_boot_0_500r <- data.frame()

for (i in unique(dot.comp_all$site)) {
  for(j in unique(dot.comp_all$plot)) {
    subset_data <- dot.comp_all[dot.comp_all$site==i&dot.comp_all$plot==j,]
    
    subset_data <- subset_data %>%
      dplyr::select(!year) %>%
      filter(!duplicated(cbind(X_mm, Y_mm)))
    
    subset_data$species_English <- as.factor(subset_data$species_English) # make species code a factor
    
    Ω0_500r <- data.frame()
    bands <- seq(0, 1900, 100)
    
    for (k in 1:nrow(subset_data)) {
      for (m in bands) {
        
        spe_target <- subset_data$species_English[k]
        
        sub_Ω0_500r <- subset_data[-k,] %>%
          mutate(x_center = subset_data$X_mm[k],
                 y_center = subset_data$Y_mm[k]) %>%
          mutate(distance = sqrt((X_mm - x_center)^2 + (Y_mm - y_center)^2)) %>%
          filter(distance > m & distance <= m + 100) %>%
          mutate(spe_identity = if_else(species_English %in% spe_target, 'cD', 'hD')) %>%
          group_by(spe_identity) %>%
          dplyr::summarise(D_value = sum(abun)) %>%
          pivot_wider(names_from = spe_identity, values_from = D_value) %>%
          mutate(species_English = spe_target,
                 x_center = subset_data$X_mm[k],
                 y_center = subset_data$Y_mm[k],
                 cover_cm2 = subset_data$cover_cm2[k],
                 high_mm = subset_data$high_mm[k],
                 band = m,
                 area = case_when(
                   x_center >= m+100&x_center <= 2000-(m+100)&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2, 
                   x_center < m+100&x_center >= m&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2),
                   y_center < m+100&y_center >= m&x_center >= m+100&x_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center >= m+100&x_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2),
                   x_center < m&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   x_center > 2000-m&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   y_center < m&x_center >= m+100&x_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2),
                   y_center > 2000-m&x_center >= m+100&x_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   x_center < m+100&x_center >= m&y_center < m+100&y_center >= m&sqrt(x_center^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2),
                   x_center < m+100&x_center >= m&y_center < m+100&y_center >= m&sqrt(x_center^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 - (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 + ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2,
                   x_center < m+100&x_center >= m&y_center > 2000-(m+100)&y_center <= 2000-m&sqrt(x_center^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2),
                   x_center < m+100&x_center >= m&y_center > 2000-(m+100)&y_center <= 2000-m&sqrt(x_center^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 -(sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 + ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2,
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center < m+100&y_center >= m&sqrt((2000-x_center)^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center < m+100&y_center >= m&sqrt((2000-x_center)^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 - (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 + ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))- atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2,
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center > 2000-(m+100)&y_center <= 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center > 2000-(m+100)&y_center <= 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 - (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 + ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2,
                   x_center < m+100&x_center >= m&y_center < m&sqrt(x_center^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2),
                   x_center < m+100&x_center >= m&y_center < m&sqrt(x_center^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + + (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 + (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 - ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2 + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2), 
                   y_center < m+100&y_center >= m&x_center < m&sqrt(x_center^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   y_center < m+100&y_center >= m&x_center < m&sqrt(x_center^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + + (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 + (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 - ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2 + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   x_center < m+100&x_center >= m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   x_center < m+100&x_center >= m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + + (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 + (sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 - ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2 + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   y_center < m+100&y_center >= m&x_center > 2000-m&sqrt((2000-x_center)^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   y_center < m+100&y_center >= m&x_center > 2000-m&sqrt((2000-x_center)^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + + (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 + (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 - ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))-atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + + (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 + (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 - ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))-atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center < m&sqrt(x_center^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center < m&sqrt(x_center^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + + (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 + (sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 - ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2 + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + + (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 + (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 - ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + + (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 + (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 - ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   x_center < m&y_center < m&sqrt(x_center^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2), 
                   x_center < m&y_center < m&sqrt(x_center^2 + y_center^2) < m + 100&sqrt(x_center^2 + y_center^2) >= m ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 - (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 + ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2, 
                   x_center < m&y_center < m&sqrt(x_center^2 + y_center^2) < m ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 - (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 + ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2 + (sqrt(m^2-y_center^2)-x_center)*y_center/2 + (sqrt(m^2-x_center^2)-y_center)*x_center/2 - ((acos(y_center/m)-atan(x_center/y_center)) + (acos(x_center/m)-atan(y_center/x_center)))/(2*pi)*pi*m^2, 
                   x_center > 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2), 
                   x_center > 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) < m + 100&sqrt((2000-x_center)^2 + y_center^2) >= m ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 - (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 + ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))-atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2, 
                   x_center > 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) < m ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 - (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 + ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))-atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + (sqrt(m^2-y_center^2)-(2000-x_center))*y_center/2 + (sqrt(m^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 - ((acos(y_center/m)-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/m)-atan(y_center/(2000-x_center))))/(2*pi)*pi*m^2, 
                   x_center < m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2), 
                   x_center < m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) < m + 100&sqrt(x_center^2 + (2000-y_center)^2) >= m ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 - (sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 + ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2, 
                   x_center < m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) < m ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 - (sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 + ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2 + (sqrt(m^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 + (sqrt(m^2-x_center^2)-(2000-y_center))*x_center/2 - ((acos((2000-y_center)/m)-atan(x_center/(2000-y_center))) + (acos(x_center/m)-atan((2000-y_center)/x_center)))/(2*pi)*pi*m^2, 
                   x_center > 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2), 
                   x_center > 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m + 100&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 - (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 + ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2, 
                   x_center > 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 - (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 + ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + (sqrt(m^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 + (sqrt(m^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 - ((acos((2000-y_center)/m)-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/m)-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*m^2
                 ))
        
        Ω0_500r <- bind_rows(Ω0_500r, sub_Ω0_500r)
        
      }
      
    } 
    
    sel_species <- subset_data %>% 
      dplyr::select(species_English) %>%
      group_by(species_English) %>% 
      mutate(indi_num = length(species_English)) %>%
      filter(indi_num >= 20) %>%
      filter(!duplicated(species_English)) %>%
      mutate(mean_cden = indi_num/(2000*2000),
             mean_hden = (nrow(subset_data)-indi_num)/(2000*2000))
    
    Ω0_500r <- Ω0_500r %>% 
      filter(species_English %in% sel_species$species_English) %>% 
      merge(sel_species) 
    
    Ω0_500r_boot <- data.frame()
    
    for (w in unique(sel_species$species_English)) {
      for(q in 1:15) {
        
        Ω0_500r_spe <- Ω0_500r[Ω0_500r$species_English==w,]
        
        n <- as.integer(nrow(Ω0_500r_spe)/2)
        
        random_rows <- sample(1:nrow(Ω0_500r_spe), n)
        Ω0_500r_sampled_data <- Ω0_500r_spe[random_rows, ]
        
        Ω0_boot_500r_spe <- Ω0_500r_sampled_data %>%
          group_by(band) %>%
          dplyr::summarise(mean_cden = mean(mean_cden),
                           mean_hden = mean(mean_hden),
                           tcD = sum(cD),
                           thD = sum(hD),
                           tarea = sum(area),
                           spe_abun = mean(indi_num)) %>%
          mutate(rtcD = tcD/tarea,
                 rthD = thD/tarea) %>%
          mutate(cΩ = rtcD/mean_cden,
                 hΩ = rthD/mean_hden,
                 #anpp_target = cover_cm2*high_mm, 
                 time = q)
        
        Ω0_500r_boot <- bind_rows(Ω0_500r_boot, Ω0_boot_500r_spe)
        
      }
    }
    
    Ω0_500r_boot <- Ω0_500r_boot %>%
      group_by(species_English, band) %>%
      dplyr::summarise(mean_cΩ = mean(cΩ),
                       #sd_cΩ = sd(cΩ),
                       conf_low_cΩ = t.test(cΩ)$conf.int[1]/2,
                       conf_high_cΩ = t.test(cΩ)$conf.int[2]/2,
                       mean_hΩ = mean(hΩ),
                       #sd_hΩ = sd(hΩ),
                       conf_low_hΩ = t.test(hΩ)$conf.int[1]/2,
                       conf_high_hΩ = t.test(hΩ)$conf.int[2]/2,
                       spe_abun = mean(indi_num))
      mutate(site = i, plot = j)
    
    Ω_all_site_boot_0_500r <- bind_rows(Ω_all_site_boot_0_500r, Ω0_500r_boot)
    
  }
  
}

write_xlsx(Ω_all_site_boot_0_500r, 'Ω_all_site_boot_0_500r.xlsx')


## species segregation index Ms ----
Ms_all_site <- data.frame() 
Com_mass_mean_msi_all_site_alpha <- data.frame()
Com_mass_mean_msi_all_site_beta <- data.frame()
scale_alpha <- c(100,150,200,250,300,350,400,450,500,600,700,800,900,1000)
scale_beta <- c(100,125,150,175,200,225,250,300,350,400,450,500)

for (i in unique(dot.comp_all$site)) {
  for(j in unique(dot.comp_all$plot)) {
    subset_data <- dot.comp_all[dot.comp_all$site==i&dot.comp_all$plot==j,]
    
    subset_data <- subset_data %>%
      dplyr::select(!year) %>%
      filter(!duplicated(cbind(X_mm, Y_mm))) %>%
      mutate(label = paste(X_mm, Y_mm, sep="-"))
    
    subset_data$species_English <- as.factor(subset_data$species_English) # make species code a factor
    
    subset_data_ppp <- ppp(subset_data$X_mm,
                           subset_data$Y_mm,
                           xrange=c(0,2000),yrange=c(0,2000),
                           marks=subset_data$species_English)
    
    voronoi <- deldir(subset_data_ppp)
    
    #write_xlsx(voronoi$summary, 'species segregation index_polygons.xlsx')
    #write_xlsx(voronoi$delsgs, 'species segregation index_delsgs.xlsx')
    #write_xlsx(voronoi$dirsgs, 'species segregation index_dirsgs.xlsx')
    
    delsgs_a <- voronoi$delsgs %>%
      dplyr::select(x_center = x1, y_center = y1, x_neigh = x2, y_neigh = y2, id_center = ind1)
    delsgs_b <- voronoi$delsgs %>%
      dplyr::select(x_center = x2, y_center = y2, x_neigh = x1, y_neigh = y1, id_center = ind2)
    
    delsgs <- bind_rows(delsgs_a, delsgs_b)
    
    #write_xlsx(delsgs, 'species segregation index_delsgs_all.xlsx')
    
    delsgs$label = paste(delsgs$x_neigh, delsgs$y_neigh, sep="-")
    #delsgs$labelcenter = paste(delsgs$x_center, delsgs$y_center, sep="-")
    #length(unique(delsgs$labelcenter)) # = nrow(subset_data) ensure all individuals as focal individuals
    
    neigh <- delsgs %>% 
      merge(subset_data %>% dplyr::select(label, species_English_neigh = species_English, cover_cm2_neigh = cover_cm2, high_mm_neigh = high_mm, abun_neigh = abun)) %>%
      dplyr::select(X_mm = x_center, Y_mm = y_center, id_center, x_neigh, y_neigh, species_English_neigh, cover_cm2_neigh, high_mm_neigh, abun_neigh) %>%
      merge(subset_data)
    #length(unique(neigh$id_center))   
    
    out_plant_id <- voronoi$dirsgs %>%
      mutate(outn = thirdv1*thirdv2) %>%
      filter(outn < 0) %>%
      dplyr::select(ind1, ind2) %>%
      pivot_longer(1:2, names_to = 'ind', values_to = 'id_center') %>%
      filter(!duplicated(id_center))
    
    neigh_sel <- neigh %>%
      mutate(outy = if_else(id_center %in% out_plant_id$id_center, 0, 1)) %>%
      filter(outy == 1) %>%
      mutate(spe_identity = if_else(species_English_neigh == species_English, '0', '1')) %>%
      mutate(spe_identity = as.numeric(spe_identity))
    
    sub_Msi <- neigh_sel %>%
      group_by(X_mm, Y_mm, species_English) %>%
      dplyr::summarise(Msi = if_else(sum(spe_identity) == length(species_English_neigh), 
                                     sum(spe_identity)*(length(unique(species_English_neigh))+1)/length(species_English_neigh)/(length(species_English_neigh)+1),
                                     sum(spe_identity)*length(unique(species_English_neigh))/length(species_English_neigh)/(length(species_English_neigh)+1)),
                       neigh_abun = length(species_English_neigh),
                       neigh_rich = length(unique(species_English_neigh)),
                       focal_Ind_mass = log(mean(cover_cm2*high_mm)+1)) 
      
    com_Msi = mean(sub_Msi$Msi)
    
    sub_Ms_all_site <- sub_Msi %>% merge(neigh_sel) %>% mutate(site = i, plot = j, com_Msi = com_Msi) 
    
    Ms_all_site <- bind_rows(Ms_all_site, sub_Ms_all_site)
    write_xlsx(Ms_all_site, 'Ms_all_site_2023.xlsx')
    
    Com_mass_mean_msi_alpha <- data.frame()

    for (k in scale_alpha) {
    sub_Ms_all_site$xt <- cut(sub_Ms_all_site$X_mm,seq(min(sub_Ms_all_site$X_mm),round(max(sub_Ms_all_site$X_mm),min(sub_Ms_all_site$X_mm)),k)) 
    sub_Ms_all_site$yt <- cut(sub_Ms_all_site$Y_mm,seq(min(sub_Ms_all_site$Y_mm),round(max(sub_Ms_all_site$Y_mm),min(sub_Ms_all_site$Y_mm)),k))
    
    sub_Com_mass_mean_msi <- sub_Ms_all_site %>%
      filter(xt != 'NA') %>%
      filter(yt != 'NA') %>%
      group_by(site, plot, xt, yt, X_mm, Y_mm, species_English) %>%
      dplyr::summarise(Msi=mean(Msi),
                       abun=mean(abun),
                       neigh_abun1=mean(neigh_abun),neigh_abun2=sum(abun_neigh),
                       heter_neigh_abun1=sum(spe_identity),heter_neigh_abun2=sum(abun_neigh*spe_identity),
                       neigh_rich=mean(neigh_rich),cover_cm2=mean(cover_cm2),high_mm=mean(high_mm)) %>%
      mutate(con_neigh_abun1=neigh_abun1-heter_neigh_abun1,con_neigh_abun2=neigh_abun2-heter_neigh_abun2) %>%
      ungroup() %>%
      group_by(site, plot, xt, yt) %>%
      dplyr::summarise(rich = length(unique(species_English)),
                       com_mass=log(sum(cover_cm2*high_mm)+1),
                       com_Msi=mean(Msi),
                       mean_neigh_abun1=mean(neigh_abun1),mean_neigh_abun2=mean(neigh_abun2),
                       mean_heter_neigh_abun1=mean(heter_neigh_abun1),mean_heter_neigh_abun2=mean(heter_neigh_abun2),
                       mean_con_neigh_abun1=mean(con_neigh_abun1),mean_con_neigh_abun2=mean(con_neigh_abun2),
                       mean_neigh_rich=mean(neigh_rich),com_abun=sum(abun)) %>%
      mutate(scale = k)
    
    Com_mass_mean_msi_alpha <- bind_rows(Com_mass_mean_msi_alpha, sub_Com_mass_mean_msi)
    
    }
    
    Com_mass_mean_msi_all_site_alpha <- bind_rows(Com_mass_mean_msi_all_site_alpha,Com_mass_mean_msi_alpha)
    
    write_xlsx(Com_mass_mean_msi_all_site_alpha, 'Com_mass_mean_msi_all_site_alpha_2023.xlsx')
    
    Com_mass_mean_msi_beta <- data.frame()
    
    for (s in scale_beta) {
      
      sub_Ms_all_site$xt1 <- cut(sub_Ms_all_site$X_mm,seq(min(sub_Ms_all_site$X_mm),round(max(sub_Ms_all_site$X_mm),min(sub_Ms_all_site$X_mm)),s)) 
      sub_Ms_all_site$yt1 <- cut(sub_Ms_all_site$Y_mm,seq(min(sub_Ms_all_site$Y_mm),round(max(sub_Ms_all_site$Y_mm),min(sub_Ms_all_site$Y_mm)),s))
      
      sub_Ms_all_site$xt2 <- cut(sub_Ms_all_site$X_mm,seq(min(sub_Ms_all_site$X_mm),round(max(sub_Ms_all_site$X_mm),min(sub_Ms_all_site$X_mm)),2*s)) 
      sub_Ms_all_site$yt2 <- cut(sub_Ms_all_site$Y_mm,seq(min(sub_Ms_all_site$Y_mm),round(max(sub_Ms_all_site$Y_mm),min(sub_Ms_all_site$Y_mm)),2*s))
      
      
      sub_Ms_all_site2 <- sub_Ms_all_site %>%
        filter(xt2 != 'NA') %>%
        filter(yt2 != 'NA')
      
      rich1 <- sub_Ms_all_site2 %>%
        group_by(site, plot, xt2, yt2, xt1, yt1) %>%
        dplyr::summarise(rich = length(unique(species_English))) %>%
        group_by(site, plot, xt2, yt2) %>%
        dplyr::summarise(rich1 = mean(rich))
      
      rich2 <- sub_Ms_all_site2 %>%
        group_by(site, plot, xt2, yt2) %>%
        dplyr::summarise(rich2 = length(unique(species_English)))
      
      beta_W <- rich1 %>%
        merge(rich2) %>%
        mutate(beta_add = rich2-rich1,
               beta_mul = rich2/rich1)
      
      spe_abu <- sub_Ms_all_site2 %>% 
        group_by(xt2, yt2, xt1, yt1, X_mm, Y_mm, species_English) %>%
        dplyr::summarise(abun = mean(abun)) %>%
        ungroup() %>%
        dplyr::select(xt2, yt2, xt1, yt1, species_English, abun) %>% 
        group_by(xt2, yt2, xt1, yt1, species_English) %>% 
        dplyr::summarise(abun = sum(abun)) %>%
        mutate(xyt2 = paste(xt2, yt2, sep="-"))
      
      beta_jac_bray <- data.frame()
      
      for (m in unique(spe_abu$xyt2)) {
        subdata <- spe_abu[spe_abu$xyt2 == m, ]
        
        subdata_jac <- subdata %>%
          ungroup() %>%
          dplyr::select(xt1, yt1, species_English, abun) %>% 
          mutate(abun = ifelse(abun > 0, 1, 0)) %>%
          pivot_wider(values_from = abun, names_from = species_English)
        
        subdata_jac <- replace(subdata_jac, is.na(subdata_jac), 0)
        sub_beta_jac <- data.frame(beta.multi(subdata_jac[, -c(1:2)], index.family = 'jaccard'))
        sub_beta_jac <- sub_beta_jac %>% mutate(xyt2 = m)
        
        subdata_bray <- subdata %>%
          ungroup() %>%
          dplyr::select(xt1, yt1, species_English, abun) %>% 
          pivot_wider(values_from = abun, names_from = species_English)
        
        subdata_bray <- replace(subdata_bray, is.na(subdata_bray), 0)
        sub_beta_bray <- data.frame(beta.multi.abund(subdata_bray[, -c(1:2)], index.family = 'bray'))
        sub_beta_bray <- sub_beta_bray %>% mutate(xyt2 = m)
        
        sub_beta_jac_bray <- sub_beta_jac %>% merge(sub_beta_bray)
        
        sub_beta_jac_bray <- separate(sub_beta_jac_bray,xyt2,into = c("xt2","yt2"),sep = "-")
        
        
        beta_jac_bray <- bind_rows(beta_jac_bray, sub_beta_jac_bray)
        
      }
      
      sub_Com_mass_mean_msi2 <- sub_Ms_all_site %>%
        filter(xt2 != 'NA') %>%
        filter(yt2 != 'NA') %>%
        group_by(site, plot, xt2, yt2, X_mm, Y_mm, species_English) %>%
        dplyr::summarise(Msi=mean(Msi),
                         abun=mean(abun),
                         neigh_abun1=mean(neigh_abun),neigh_abun2=sum(abun_neigh),
                         heter_neigh_abun1=sum(spe_identity),heter_neigh_abun2=sum(abun_neigh*spe_identity),
                         neigh_rich=mean(neigh_rich),cover_cm2=mean(cover_cm2),high_mm=mean(high_mm)) %>%
        mutate(con_neigh_abun1=neigh_abun1-heter_neigh_abun1,con_neigh_abun2=neigh_abun2-heter_neigh_abun2) %>%
        ungroup() %>%
        group_by(site, plot, xt2, yt2) %>%
        dplyr::summarise(rich = length(unique(species_English)),
                         com_mass=log(sum(cover_cm2*high_mm)+1),
                         com_Msi=mean(Msi),
                         mean_neigh_abun1=mean(neigh_abun1),mean_neigh_abun2=mean(neigh_abun2),
                         mean_heter_neigh_abun1=mean(heter_neigh_abun1),mean_heter_neigh_abun2=mean(heter_neigh_abun2),
                         mean_con_neigh_abun1=mean(con_neigh_abun1),mean_con_neigh_abun2=mean(con_neigh_abun2),
                         mean_neigh_rich=mean(neigh_rich),com_abun=sum(abun)) %>%
        mutate(scale = 2*s) %>%
        merge(beta_W) %>% merge(beta_jac_bray)
      
      Com_mass_mean_msi_beta <- bind_rows(Com_mass_mean_msi_beta, sub_Com_mass_mean_msi2)
      
    }
    
    Com_mass_mean_msi_all_site_beta <- bind_rows(Com_mass_mean_msi_all_site_beta,Com_mass_mean_msi_beta)
    write_xlsx(Com_mass_mean_msi_all_site_beta, 'Com_mass_mean_msi_all_site_beta_2023.xlsx')
  }
}

## expected value based on CSR of species segregation index Ms ----
Com_mass_mean_msi_RS_exp_all_sites <- data.frame()
scale <- c(100,150,200,250,300,350,400,450,500,600,700,800,900,1000)

for (i in unique(dot.comp_all$site)) {
  for(j in unique(dot.comp_all$plot)) {
    subset_data <- dot.comp_all[dot.comp_all$site==i&dot.comp_all$plot==j,]
    
    subset_data <- subset_data %>%
      dplyr::select(!year) %>%
      filter(!duplicated(cbind(X_mm, Y_mm)))
    
    subset_data <- subset_data %>%
      mutate(id = 1:nrow(subset_data))
    
    subset_data$species_English <- as.factor(subset_data$species_English) # make species code a factor
    
    Com_mass_mean_msi_all_time <- data.frame()
    
    for (m in 1:199) {
      #for (w in scale) {
      
      n = nrow(subset_data)
      
      random_spe <- subset_data %>%
        dplyr::select(species_English, cover_cm2, high_mm, abun) %>%
        mutate(id = sample(1:n, n))
      
      subset_data2 <- subset_data %>%
        dplyr::select(-c(species_English, cover_cm2, high_mm, abun))
      
      subset_data_csr <- subset_data2 %>%
        merge(random_spe) %>%
        dplyr::select(-id) %>%
        mutate(label = paste(X_mm, Y_mm, sep="-"))
      
      subset_data_csr_ppp <- ppp(subset_data_csr$X_mm,
                                 subset_data_csr$Y_mm,
                                 xrange=c(0,2000),yrange=c(0,2000),
                                 marks=subset_data_csr$species_English)
      
      voronoi <- deldir(subset_data_csr_ppp)
      
      delsgs_a <- voronoi$delsgs %>%
        dplyr::select(x_center = x1, y_center = y1, x_neigh = x2, y_neigh = y2, id_center = ind1)
      delsgs_b <- voronoi$delsgs %>%
        dplyr::select(x_center = x2, y_center = y2, x_neigh = x1, y_neigh = y1, id_center = ind2)
      
      delsgs <- bind_rows(delsgs_a, delsgs_b)
      
      delsgs$label = paste(delsgs$x_neigh, delsgs$y_neigh, sep="-")
      #delsgs$labelcenter = paste(delsgs$x_center, delsgs$y_center, sep="-")
      #length(unique(delsgs$labelcenter)) # = nrow(subset_data) ensure all individuals as focal individuals
      
      neigh <- delsgs %>% 
        merge(subset_data_csr %>% dplyr::select(label, species_English_neigh = species_English, cover_cm2_neigh = cover_cm2, high_mm_neigh = high_mm , abun_neigh = abun)) %>%
        dplyr::select(X_mm = x_center, Y_mm = y_center, id_center, x_neigh, y_neigh, species_English_neigh, cover_cm2_neigh, high_mm_neigh, abun_neigh) %>%
        merge(subset_data_csr)
      #length(unique(neigh$id_center))   
      
      out_plant_id <- voronoi$dirsgs %>%
        mutate(outn = thirdv1*thirdv2) %>%
        filter(outn < 0) %>%
        dplyr::select(ind1, ind2) %>%
        pivot_longer(1:2, names_to = 'ind', values_to = 'id_center') %>%
        filter(!duplicated(id_center))
      
      neigh_sel <- neigh %>%
        mutate(outy = if_else(id_center %in% out_plant_id$id_center, 0, 1)) %>%
        filter(outy == 1) %>%
        mutate(spe_identity = if_else(species_English_neigh == species_English, '0', '1')) %>%
        mutate(spe_identity = as.numeric(spe_identity))
      
      sub_Msi <- neigh_sel %>%
        group_by(X_mm, Y_mm, species_English) %>%
        dplyr::summarise(Msi = if_else(sum(spe_identity) == length(species_English_neigh), 
                                       sum(spe_identity)*(length(unique(species_English_neigh))+1)/length(species_English_neigh)/(length(species_English_neigh)+1),
                                       sum(spe_identity)*length(unique(species_English_neigh))/length(species_English_neigh)/(length(species_English_neigh)+1)),
                         neigh_abun = length(species_English_neigh),
                         neigh_rich = length(unique(species_English_neigh))) 
      
      com_Msi = mean(sub_Msi$Msi)
      
      sub_Ms_all_site <- sub_Msi %>% merge(neigh_sel) #%>% mutate(site = i, plot = j, com_Msi = com_Msi) 
      
      #Ms_all_site <- bind_rows(Ms_all_site, sub_Ms_all_site)
      
      Com_mass_mean_msi <- data.frame()
      
      for (k in scale) {
        sub_Ms_all_site$xt <- cut(sub_Ms_all_site$X_mm,seq(min(sub_Ms_all_site$X_mm),round(max(sub_Ms_all_site$X_mm),min(sub_Ms_all_site$X_mm)),k)) 
        sub_Ms_all_site$yt <- cut(sub_Ms_all_site$Y_mm,seq(min(sub_Ms_all_site$Y_mm),round(max(sub_Ms_all_site$Y_mm),min(sub_Ms_all_site$Y_mm)),k))
        
        sub_Com_mass_mean_msi <- sub_Ms_all_site %>%
          filter(xt != 'NA') %>%
          filter(yt != 'NA') %>%
          dplyr::group_by(xt, yt, X_mm, Y_mm, species_English) %>%
          dplyr::summarise(Msi=mean(Msi),
                           neigh_abun1=mean(neigh_abun),neigh_abun2=sum(abun_neigh),
                           heter_neigh_abun1=sum(spe_identity),heter_neigh_abun2=sum(abun_neigh*spe_identity),
                           neigh_rich=mean(neigh_rich),cover_cm2=mean(cover_cm2),high_mm=mean(high_mm)) %>%
          dplyr::mutate(con_neigh_abun1=neigh_abun1-heter_neigh_abun1,con_neigh_abun2=neigh_abun2-heter_neigh_abun2) %>%
          dplyr::mutate(struc_rich=if_else(con_neigh_abun1 == 0, neigh_rich+1, neigh_rich)) %>%
          dplyr::ungroup() %>%
          dplyr::group_by(xt, yt) %>%
          dplyr::summarise(rich = length(unique(species_English)),
                           com_mass=log(sum(cover_cm2*high_mm)+1),
                           com_Msi=mean(Msi),
                           mean_neigh_abun1=mean(neigh_abun1),mean_neigh_abun2=mean(neigh_abun2),
                           mean_heter_neigh_abun1=mean(heter_neigh_abun1),mean_heter_neigh_abun2=mean(heter_neigh_abun2),
                           mean_con_neigh_abun1=mean(con_neigh_abun1),mean_con_neigh_abun2=mean(con_neigh_abun2),
                           mean_neigh_rich=mean(neigh_rich),
                           mean_struc_rich=mean(struc_rich)) %>%
          dplyr::mutate(scale = k)
        
        Com_mass_mean_msi <- bind_rows(Com_mass_mean_msi, sub_Com_mass_mean_msi)
      }
      
      Com_mass_mean_msi <- Com_mass_mean_msi %>% mutate(time = m)
      
      Com_mass_mean_msi_all_time <- bind_rows(Com_mass_mean_msi_all_time,Com_mass_mean_msi)
      
    }
    
    Com_mass_mean_msi_exp <- Com_mass_mean_msi_all_time %>%
      group_by(xt, yt, scale) %>%
      dplyr::summarise(mean_rich = mean(rich),
                       sd_rich = sd(rich),
                       #conf_low_rich = t.test(rich)$conf.int[1],
                       #conf_high_rich = t.test(rich)$conf.int[2],
                       mean_com_mass = mean(com_mass),
                       sd_com_mass = sd(com_mass),
                       conf_low_com_mass = t.test(com_mass)$conf.int[1],
                       conf_high_com_mass = t.test(com_mass)$conf.int[2],
                       mean_com_Msi = mean(com_Msi),
                       sd_com_Msi = sd(com_Msi),
                       conf_low_com_Msi = t.test(com_Msi)$conf.int[1],
                       conf_high_com_Msi = t.test(com_Msi)$conf.int[2],
                       mean_neigh_abun1 = mean(mean_neigh_abun1),
                       sd_neigh_abun1 = sd(mean_neigh_abun1),
                       #conf_low_neigh_abun1 = t.test(mean_neigh_abun1)$conf.int[1],
                       #conf_high_neigh_abun1 = t.test(mean_neigh_abun1)$conf.int[2],
                       aver_neigh_abun2 = mean(mean_neigh_abun2),
                       sd_neigh_abun2 = sd(mean_neigh_abun2),
                       conf_low_neigh_abun2 = t.test(mean_neigh_abun2)$conf.int[1],
                       conf_high_neigh_abun2 = t.test(mean_neigh_abun2)$conf.int[2],
                       aver_heter_neigh_abun1 = mean(mean_heter_neigh_abun1),
                       sd_heter_neigh_abun1 = sd(mean_heter_neigh_abun1),
                       conf_low_heter_neigh_abun1 = t.test(mean_heter_neigh_abun1)$conf.int[1],
                       conf_high_heter_neigh_abun1 = t.test(mean_heter_neigh_abun1)$conf.int[2],
                       aver_heter_neigh_abun2 = mean(mean_heter_neigh_abun2),
                       sd_heter_neigh_abun2 = sd(mean_heter_neigh_abun2),
                       conf_low_heter_neigh_abun2 = t.test(mean_heter_neigh_abun2)$conf.int[1],
                       conf_high_heter_neigh_abun2 = t.test(mean_heter_neigh_abun2)$conf.int[2],
                       aver_con_neigh_abun1 = mean(mean_con_neigh_abun1),
                       sd_con_neigh_abun1 = sd(mean_con_neigh_abun1),
                       conf_low_con_neigh_abun1 = t.test(mean_con_neigh_abun1)$conf.int[1],
                       conf_high_con_neigh_abun1 = t.test(mean_con_neigh_abun1)$conf.int[2],
                       aver_con_neigh_abun2 = mean(mean_con_neigh_abun2),
                       sd_con_neigh_abun2 = sd(mean_con_neigh_abun2),
                       conf_low_con_neigh_abun2 = t.test(mean_con_neigh_abun2)$conf.int[1],
                       conf_high_con_neigh_abun2 = t.test(mean_con_neigh_abun2)$conf.int[2],
                       aver_neigh_rich = mean(mean_neigh_rich),
                       sd_neigh_rich = sd(mean_neigh_rich),
                       conf_low_neigh_rich = t.test(mean_neigh_rich)$conf.int[1],
                       conf_high_neigh_rich = t.test(mean_neigh_rich)$conf.int[2],
                       aver_struc_rich = mean(mean_struc_rich),
                       sd_struc_rich = sd(mean_struc_rich),
                       conf_low_struc_rich = t.test(mean_struc_rich)$conf.int[1],
                       conf_high_struc_rich = t.test(mean_struc_rich)$conf.int[2]) %>%
      mutate(site = i, plot = j)
    
    Com_mass_mean_msi_RS_exp_all_sites <- bind_rows(Com_mass_mean_msi_RS_exp_all_sites, Com_mass_mean_msi_exp)
    
    write_xlsx(Com_mass_mean_msi_RS_exp_all_sites, 'Com_mass_mean_msi_RS_exp_all_sites_2023.xlsx')
    
  }
}

## species segregation index Ms distributed ----
sel_species <- read_excel("sel_spe_2023.xlsx",sheet=1,na="NA")
Ms_distribution_all_site <- data.frame()

for (i in unique(dot.comp_all$site)) {
  for(j in unique(dot.comp_all$plot)) {
    subset_data <- dot.comp_all[dot.comp_all$site==i&dot.comp_all$plot==j,]
    
    subset_data <- subset_data %>%
      dplyr::select(!year) %>%
      filter(!duplicated(cbind(X_mm, Y_mm)))
    
    subset_data <- subset_data %>%
      mutate(id = 1:nrow(subset_data))
    
    subset_data$species_English <- as.factor(subset_data$species_English) # make species code a factor
    
    sub_Msi_sel <- sel_species[sel_species$site==i&sel_species$plot==j,]
    
    Ms_distribution_plot <- data.frame()
    
    for (m in 1:99) {
      #for (w in scale) {
      
      n = nrow(subset_data)
      
      random_spe <- subset_data %>%
        dplyr::select(species_English, cover_cm2, high_mm, abun) %>%
        mutate(id = sample(1:n, n))
      
      subset_data2 <- subset_data %>%
        dplyr::select(-c(species_English, cover_cm2, high_mm, abun))
      
      subset_data_csr <- subset_data2 %>%
        merge(random_spe) %>%
        dplyr::select(-id) %>%
        mutate(label = paste(X_mm, Y_mm, sep="-"))
      
      subset_data_csr_ppp <- ppp(subset_data_csr$X_mm,
                                 subset_data_csr$Y_mm,
                                 xrange=c(0,2000),yrange=c(0,2000),
                                 marks=subset_data_csr$species_English)
      
      voronoi <- deldir(subset_data_csr_ppp)
      
      delsgs_a <- voronoi$delsgs %>%
        dplyr::select(x_center = x1, y_center = y1, x_neigh = x2, y_neigh = y2, id_center = ind1)
      delsgs_b <- voronoi$delsgs %>%
        dplyr::select(x_center = x2, y_center = y2, x_neigh = x1, y_neigh = y1, id_center = ind2)
      
      delsgs <- bind_rows(delsgs_a, delsgs_b)
      
      delsgs$label = paste(delsgs$x_neigh, delsgs$y_neigh, sep="-")
      #delsgs$labelcenter = paste(delsgs$x_center, delsgs$y_center, sep="-")
      #length(unique(delsgs$labelcenter)) # = nrow(subset_data) ensure all individuals as focal individuals
      
      neigh <- delsgs %>% 
        merge(subset_data_csr %>% dplyr::select(label, species_English_neigh = species_English, cover_cm2_neigh = cover_cm2, high_mm_neigh = high_mm , abun_neigh = abun)) %>%
        dplyr::select(X_mm = x_center, Y_mm = y_center, id_center, x_neigh, y_neigh, species_English_neigh, cover_cm2_neigh, high_mm_neigh, abun_neigh) %>%
        merge(subset_data_csr)
      #length(unique(neigh$id_center))   
      
      out_plant_id <- voronoi$dirsgs %>%
        mutate(outn = thirdv1*thirdv2) %>%
        filter(outn < 0) %>%
        dplyr::select(ind1, ind2) %>%
        pivot_longer(1:2, names_to = 'ind', values_to = 'id_center') %>%
        filter(!duplicated(id_center))
      
      neigh_sel <- neigh %>%
        mutate(outy = if_else(id_center %in% out_plant_id$id_center, 0, 1)) %>%
        filter(outy == 1) %>%
        mutate(spe_identity = if_else(species_English_neigh == species_English, '0', '1')) %>%
        mutate(spe_identity = as.numeric(spe_identity))
      
      sub_Msi <- neigh_sel %>%
        group_by(X_mm, Y_mm, species_English) %>%
        dplyr::summarise(Msi = if_else(sum(spe_identity) == length(species_English_neigh), 
                                       sum(spe_identity)*(length(unique(species_English_neigh))+1)/length(species_English_neigh)/(length(species_English_neigh)+1),
                                       sum(spe_identity)*length(unique(species_English_neigh))/length(species_English_neigh)/(length(species_English_neigh)+1)),
                         neigh_abun = length(species_English_neigh),
                         neigh_rich = length(unique(species_English_neigh))) 
      
      #com_Msi = mean(sub_Msi$Msi)
      #sub_Ms_all_site <- sub_Msi %>% merge(neigh_sel)
      
      sub_Msi <- sub_Msi %>%
        filter(species_English %in% sub_Msi_sel$species_English) %>%
        mutate(time = as.factor(m)) 
      
      Ms_distribution_plot <- bind_rows(Ms_distribution_plot, sub_Msi)
    }
    
    sub_Ms_distribution <- data.frame()
    
    for (s in unique(sub_Msi_sel$species_English)) {
      
      sub_Ms_distribution_plot <- Ms_distribution_plot[Ms_distribution_plot$species_English==s,]
      
      global_min <- min(sub_Ms_distribution_plot$Msi)
      global_max <- max(sub_Ms_distribution_plot$Msi)
      
      x_grid <- seq(global_min, global_max, length.out = 100)
      
      sub_Ms_distribution_plot <- sub_Ms_distribution_plot %>% dplyr::select(time, Msi)
  
      list_sub_Ms_distribution_plot <- split(sub_Ms_distribution_plot$Msi, sub_Ms_distribution_plot$time)
      
      densities_spe <- lapply(list_sub_Ms_distribution_plot, function(sim) {
        dens <- density(sim, n = 100)
        approx(dens$x, dens$y, xout = x_grid)$y
      })
      
      density_matrix_spe <- do.call(rbind, densities_spe)
      density_matrix_spe[is.na(density_matrix_spe)] <- 0
      
      mean_density <- colMeans(density_matrix_spe)
      lower_ci <- apply(density_matrix_spe, 2, quantile, probs = 0.025)
      upper_ci <- apply(density_matrix_spe, 2, quantile, probs = 0.975)
      
      sub_Ms_distribution_spe <- data.frame(x = x_grid, mean = mean_density, lower = lower_ci, upper = upper_ci, 
                                            species_English = s)
      
      sub_Ms_distribution <- bind_rows(sub_Ms_distribution, sub_Ms_distribution_spe)
    }
    
    sub_Ms_distribution <- sub_Ms_distribution %>% mutate(site = i, plot = j)
    
    Ms_distribution_all_site <- bind_rows(Ms_distribution_all_site, sub_Ms_distribution)
    
    write_xlsx(Ms_distribution_all_site, 'Ms_distribution_all_site_2023.xlsx')
  }
}
## maybe for species --- expected value based on CSR of species segregation index Ms ----
Ms_all_CSR_exp_site <- data.frame() 
#Ms_all_CSR_exp_com_scaling <- data.frame()
scale <- c(100,150,200,250,300,350,400,450,500,600,700,800,900,1000)

for (i in unique(dot.comp_all$site)) {
  for(j in unique(dot.comp_all$plot)) {
    subset_data <- dot.comp_all[dot.comp_all$site==i&dot.comp_all$plot==j,]
    
    subset_data <- subset_data %>%
      dplyr::select(!year) %>%
      filter(!duplicated(cbind(X_mm, Y_mm)))
    
    subset_data <- subset_data %>%
      mutate(id = 1:nrow(subset_data))
    
    subset_data$species_English <- as.factor(subset_data$species_English) # make species code a factor
    
    #sub_Msi_sel <- sel_species[sel_species$site==i&sel_species$plot==j,]
    
    sub_Ms_csr <- data.frame()
    
    for (k in 1:99) {
      #for (w in scale) {
      
      n = nrow(subset_data)
      
      random_spe <- subset_data %>%
        dplyr::select(species_English, cover_cm2, high_mm, abun) %>%
        mutate(id = sample(1:n, n))
      
      subset_data2 <- subset_data %>%
        dplyr::select(-c(species_English, cover_cm2, high_mm, abun))
      
      subset_data_csr <- subset_data2 %>%
        merge(random_spe) %>%
        dplyr::select(-id) %>%
        mutate(label = paste(X_mm, Y_mm, sep="-"))
      
      subset_data_csr_ppp <- ppp(subset_data_csr$X_mm,
                                 subset_data_csr$Y_mm,
                                 xrange=c(0,2000),yrange=c(0,2000),
                                 marks=subset_data_csr$species_English)
      
      voronoi <- deldir(subset_data_csr_ppp)
      
      delsgs_a <- voronoi$delsgs %>%
        dplyr::select(x_center = x1, y_center = y1, x_neigh = x2, y_neigh = y2, id_center = ind1)
      delsgs_b <- voronoi$delsgs %>%
        dplyr::select(x_center = x2, y_center = y2, x_neigh = x1, y_neigh = y1, id_center = ind2)
      
      delsgs <- bind_rows(delsgs_a, delsgs_b)
      
      delsgs$label = paste(delsgs$x_neigh, delsgs$y_neigh, sep="-")
      #delsgs$labelcenter = paste(delsgs$x_center, delsgs$y_center, sep="-")
      #length(unique(delsgs$labelcenter)) # = nrow(subset_data) ensure all individuals as focal individuals
      
      neigh <- delsgs %>% 
        merge(subset_data_csr %>% dplyr::select(label, species_English_neigh = species_English, cover_cm2_neigh = cover_cm2, high_mm_neigh = high_mm , abun_neigh = abun)) %>%
        dplyr::select(X_mm = x_center, Y_mm = y_center, id_center, x_neigh, y_neigh, species_English_neigh, cover_cm2_neigh, high_mm_neigh, abun_neigh) %>%
        merge(subset_data_csr)
      #length(unique(neigh$id_center))   
      
      out_plant_id <- voronoi$dirsgs %>%
        mutate(outn = thirdv1*thirdv2) %>%
        filter(outn < 0) %>%
        dplyr::select(ind1, ind2) %>%
        pivot_longer(1:2, names_to = 'ind', values_to = 'id_center') %>%
        filter(!duplicated(id_center))
      
      neigh_sel <- neigh %>%
        mutate(outy = if_else(id_center %in% out_plant_id$id_center, 0, 1)) %>%
        filter(outy == 1) %>%
        mutate(spe_identity = if_else(species_English_neigh == species_English, '0', '1')) %>%
        mutate(spe_identity = as.numeric(spe_identity))
      
      sub_Msi <- neigh_sel %>%
        group_by(X_mm, Y_mm, species_English) %>%
        dplyr::summarise(Msi = if_else(sum(spe_identity) == length(species_English_neigh), 
                                       sum(spe_identity)*(length(unique(species_English_neigh))+1)/length(species_English_neigh)/(length(species_English_neigh)+1),
                                       sum(spe_identity)*length(unique(species_English_neigh))/length(species_English_neigh)/(length(species_English_neigh)+1))) %>%
        mutate(site = i, plot = j, time = k)
      
      com_Msi = mean(sub_Msi$Msi)
      
      sub_Ms_csr_time <- sub_Msi %>%
        mutate(com_Msi = com_Msi) %>%
        group_by(species_English) %>%
        dplyr::summarise(spe_Msi = mean(Msi),
                         com_Msi = mean(com_Msi)) %>%
        mutate(site = i, plot = j, time = k)
      
      
      sub_Ms_csr <- bind_rows(sub_Ms_csr, sub_Ms_csr_time)
      
      #sub_Msi$xt <- cut(sub_Msi$X_mm,seq(min(sub_Msi$X_mm),round(max(sub_Msi$X_mm),min(sub_Msi$X_mm)),w)) 
      #sub_Msi$yt <- cut(sub_Msi$Y_mm,seq(min(sub_Msi$Y_mm),round(max(sub_Msi$Y_mm),min(sub_Msi$Y_mm)),w))
        
      #sub2_com_mass_mean_msi_csr <- sub_Msi %>%
      #    filter(xt != 'NA') %>%
      #    filter(yt != 'NA') %>%
      #    group_by(site, plot, time, xt, yt, species_English) %>%
      #    dplyr::summarise(Msi=mean(Msi)) %>%
      #    ungroup() %>%
      #    group_by(site, plot, time, xt, yt) %>%
      #    dplyr::summarise(mean_Msi=mean(Msi)) %>%
      #    mutate(scale = w)

      #sub_com_mass_mean_msi_csr_scaling <- bind_rows(sub_com_mass_mean_msi_csr_scaling, sub2_com_mass_mean_msi_csr)
      }
    #}
    
    #sub_com_mass_mean_msi_csr_scaling_exp <- sub_com_mass_mean_msi_csr_scaling %>%
    #  group_by(site, plot, scale, xt, yt) %>%
    #  dplyr::summarise(mean_com_Msi=mean(mean_Msi),
    #                   sd_com_Msi=sd(mean_Msi),
    #                   conf_low_com_Msi = t.test(mean_Msi)$conf.int[1],
    #                   conf_high_com_Msi = t.test(mean_Msi)$conf.int[2])
    
    #Ms_all_CSR_exp_com_scaling <- bind_rows(Ms_all_CSR_exp_com_scaling, sub_com_mass_mean_msi_csr_scaling)
    
    sub_Ms_spe_csr_exp <- sub_Ms_csr %>%
      group_by(site, plot, species_English) %>%
      dplyr::summarise(mean_spe_Msi = mean(spe_Msi),
                       sd_spe_Msi = sd(spe_Msi),
                       conf_low_spe_Msi = t.test(spe_Msi)$conf.int[1],
                       conf_high_spe_Msi = t.test(spe_Msi)$conf.int[2])
    
    sub_Ms_com_csr_exp <- sub_Ms_csr %>%
      group_by(site, plot, time) %>%
      dplyr::summarise(com_Msi = mean(com_Msi)) 
    
    sub_Ms_com_csr_exp <- sub_Ms_com_csr_exp %>%
      group_by(site, plot) %>%   
      dplyr::summarise(mean_com_Msi = mean(com_Msi),
                       sd_com_Msi = sd(com_Msi),                 
                       conf_low_com_Msi = t.test(com_Msi)$conf.int[1],
                       conf_high_com_Msi = t.test(com_Msi)$conf.int[2])
    
    sub_Ms_csr_exp <- sub_Ms_spe_csr_exp %>%
      merge(sub_Ms_com_csr_exp)
    
    Ms_all_CSR_exp_site <- bind_rows(Ms_all_CSR_exp_site, sub_Ms_csr_exp)
    
    write_xlsx(Ms_all_CSR_exp_site, 'Ms_all_CSR_exp_site_2023.xlsx')
    #write_xlsx(Ms_all_CSR_exp_com_scaling, 'Ms_all_CSR_exp_com_scaling.xlsx')
  }
}

## expected species richness range based on CSR ----
setwd("D:/OneDrive/Rdata/7 spatial patterns of grassland plant species")
dot.comp_all <- read_excel("all_points_across_six_sites_2023.xlsx",sheet=1,na="NA")
head(dot.comp_all)
headTail(dot.comp_all)
str(dot.comp_all)
summary(dot.comp_all)

dot.comp_all$X_mm <- as.numeric(dot.comp_all$X_mm)
dot.comp_all$Y_mm <- as.numeric(dot.comp_all$Y_mm)

Com_mass_mean_msi_RS_exp_all_sites <- data.frame()
scale <- c(100,150,200,250,300,350,400,450,500,600,700,800,900,1000)

for (i in unique(dot.comp_all$site)) {
  
  Com_mass_mean_msi_all_time_site <- data.frame()
  
  for(j in unique(dot.comp_all$plot)) {
    subset_data <- dot.comp_all[dot.comp_all$site==i&dot.comp_all$plot==j,]
    
    subset_data <- subset_data %>%
      dplyr::select(!year) %>%
      filter(!duplicated(cbind(X_mm, Y_mm)))
    
    subset_data <- subset_data %>%
      mutate(id = 1:nrow(subset_data))
    
    subset_data$species_English <- as.factor(subset_data$species_English) # make species code a factor
    
    Com_mass_mean_msi_all_time <- data.frame()
    
    for (m in 1:199) {
      #for (w in scale) {
      
      n = nrow(subset_data)
      
      random_spe <- subset_data %>%
        dplyr::select(species_English, cover_cm2, high_mm, abun) %>%
        mutate(id = sample(1:n, n))
      
      subset_data2 <- subset_data %>%
        dplyr::select(-c(species_English, cover_cm2, high_mm, abun))
      
      subset_data_csr <- subset_data2 %>%
        merge(random_spe) %>%
        dplyr::select(-id) %>%
        mutate(label = paste(X_mm, Y_mm, sep="-"))
      
      subset_data_csr_ppp <- ppp(subset_data_csr$X_mm,
                                 subset_data_csr$Y_mm,
                                 xrange=c(0,2000),yrange=c(0,2000),
                                 marks=subset_data_csr$species_English)
      
      voronoi <- deldir(subset_data_csr_ppp)
      
      delsgs_a <- voronoi$delsgs %>%
        dplyr::select(x_center = x1, y_center = y1, x_neigh = x2, y_neigh = y2, id_center = ind1)
      delsgs_b <- voronoi$delsgs %>%
        dplyr::select(x_center = x2, y_center = y2, x_neigh = x1, y_neigh = y1, id_center = ind2)
      
      delsgs <- bind_rows(delsgs_a, delsgs_b)
      
      delsgs$label = paste(delsgs$x_neigh, delsgs$y_neigh, sep="-")
      #delsgs$labelcenter = paste(delsgs$x_center, delsgs$y_center, sep="-")
      #length(unique(delsgs$labelcenter)) # = nrow(subset_data) ensure all individuals as focal individuals
      
      neigh <- delsgs %>% 
        merge(subset_data_csr %>% dplyr::select(label, species_English_neigh = species_English, cover_cm2_neigh = cover_cm2, high_mm_neigh = high_mm , abun_neigh = abun)) %>%
        dplyr::select(X_mm = x_center, Y_mm = y_center, id_center, x_neigh, y_neigh, species_English_neigh, cover_cm2_neigh, high_mm_neigh, abun_neigh) %>%
        merge(subset_data_csr)
      #length(unique(neigh$id_center))   
      
      out_plant_id <- voronoi$dirsgs %>%
        mutate(outn = thirdv1*thirdv2) %>%
        filter(outn < 0) %>%
        dplyr::select(ind1, ind2) %>%
        pivot_longer(1:2, names_to = 'ind', values_to = 'id_center') %>%
        filter(!duplicated(id_center))
      
      neigh_sel <- neigh %>%
        mutate(outy = if_else(id_center %in% out_plant_id$id_center, 0, 1)) %>%
        filter(outy == 1) %>%
        mutate(spe_identity = if_else(species_English_neigh == species_English, '0', '1')) %>%
        mutate(spe_identity = as.numeric(spe_identity))
      
      sub_Msi <- neigh_sel %>%
        group_by(X_mm, Y_mm, species_English) %>%
        dplyr::summarise(Msi = if_else(sum(spe_identity) == length(species_English_neigh), 
                                       sum(spe_identity)*(length(unique(species_English_neigh))+1)/length(species_English_neigh)/(length(species_English_neigh)+1),
                                       sum(spe_identity)*length(unique(species_English_neigh))/length(species_English_neigh)/(length(species_English_neigh)+1)),
                         neigh_abun = length(species_English_neigh),
                         neigh_rich = length(unique(species_English_neigh))) 
      
      sub_Ms_all_site <- sub_Msi %>% merge(neigh_sel) #%>% mutate(site = i, plot = j, com_Msi = com_Msi) 
      
      #Ms_all_site <- bind_rows(Ms_all_site, sub_Ms_all_site)
      
      Com_mass_mean_msi <- data.frame()
      
      for (k in scale) {
        sub_Ms_all_site$xt <- cut(sub_Ms_all_site$X_mm,seq(min(sub_Ms_all_site$X_mm),round(max(sub_Ms_all_site$X_mm),min(sub_Ms_all_site$X_mm)),k)) 
        sub_Ms_all_site$yt <- cut(sub_Ms_all_site$Y_mm,seq(min(sub_Ms_all_site$Y_mm),round(max(sub_Ms_all_site$Y_mm),min(sub_Ms_all_site$Y_mm)),k))
        
        sub_Com_mass_mean_msi <- sub_Ms_all_site %>%
          filter(xt != 'NA') %>%
          filter(yt != 'NA') %>%
          dplyr::group_by(xt, yt, X_mm, Y_mm, species_English) %>%
          dplyr::summarise(Msi=mean(Msi),cover_cm2=mean(cover_cm2),high_mm=mean(high_mm)) %>%
          dplyr::ungroup() %>%
          dplyr::group_by(xt, yt) %>%
          dplyr::summarise(rich = length(unique(species_English)),
                           com_mass=log(sum(cover_cm2*high_mm/10)+1),
                           com_Msi=mean(Msi)) %>%
          ungroup() %>%
          dplyr::mutate(scale = k)
        
        Com_mass_mean_msi <- bind_rows(Com_mass_mean_msi, sub_Com_mass_mean_msi)
      }
      
      Com_mass_mean_msi <- Com_mass_mean_msi %>% mutate(time = m)
      Com_mass_mean_msi_all_time <- bind_rows(Com_mass_mean_msi_all_time,Com_mass_mean_msi)
      
    }
    
    Com_mass_mean_msi_all_time <- Com_mass_mean_msi_all_time %>% mutate(plot = j)
    Com_mass_mean_msi_all_time_site <- bind_rows(Com_mass_mean_msi_all_time_site, Com_mass_mean_msi_all_time)
    
  }
  
   Com_mass_mean_msi_exp <- Com_mass_mean_msi_all_time_site %>%
      dplyr::group_by(time,scale) %>%
      dplyr::summarise(mean_rich = mean(rich),
                       mean_com_Msi = mean(com_Msi),
                       mean_com_mass = mean(com_mass),
                       max_rich = max(rich),
                       max_com_Msi = max(com_Msi),
                       max_com_mass = max(com_mass),
                       min_rich = min(rich),
                       min_com_Msi = min(com_Msi),
                       min_com_mass = min(com_mass),
                       med_rich = median(rich),
                       med_com_Msi = median(com_Msi),
                       med_com_mass = median(com_mass)) %>%
      ungroup() %>%
      mutate(D_rich = max_rich-min_rich,
             D_com_Msi = max_com_Msi-min_com_Msi,
             D_com_mass = max_com_mass-min_com_mass) %>%
      dplyr::group_by(scale) %>% 
      dplyr::summarise(mean_D_rich = mean(D_rich),
                       sd_D_rich = sd(D_rich),
                       mean_D_com_Msi = mean(D_com_Msi),
                       sd_D_com_Msi = sd(D_com_Msi),
                       mean_D_com_mass = mean(D_com_mass),
                       sd_D_com_mass = sd(D_com_mass),
                       
                       mean_rich0 = mean(mean_rich),
                       sd_mean_rich = sd(mean_rich),
                       #conf_low_mean_rich = t.test(mean_rich)$conf.int[1],
                       #conf_high_mean_rich = t.test(mean_rich)$conf.int[2],
                       mean_com_Msi0 = mean(mean_com_Msi),
                       sd_mean_com_Msi = sd(mean_com_Msi),
                       #conf_low_mean_com_Msi = t.test(mean_com_Msi)$conf.int[1],
                       #conf_high_mean_com_Msi = t.test(mean_com_Msi)$conf.int[2],
                       mean_com_mass0 = mean(mean_com_mass),
                       sd_mean_com_mass = sd(mean_com_mass),
                       #conf_low_mean_com_mass = t.test(mean_com_mass)$conf.int[1],
                       #conf_high_mean_com_mass = t.test(mean_com_mass)$conf.int[2],
                       
                       mean_max_rich0 = mean(max_rich),
                       sd_max_rich = sd(max_rich),
                       #conf_low_max_rich = t.test(max_rich)$conf.int[1],
                       #conf_high_max_rich = t.test(max_rich)$conf.int[2],
                       mean_max_com_Msi0 = mean(max_com_Msi),
                       sd_max_com_Msi = sd(max_com_Msi),
                       #conf_low_max_com_Msi = t.test(max_com_Msi)$conf.int[1],
                       #conf_high_max_com_Msi = t.test(max_com_Msi)$conf.int[2],
                       mean_max_com_mass0 = mean(max_com_mass),
                       sd_max_com_mass = sd(max_com_mass),
                       #conf_low_max_com_mass = t.test(max_com_mass)$conf.int[1],
                       #conf_high_max_com_mass = t.test(max_com_mass)$conf.int[2],
                       
                       mean_min_rich0 = mean(min_rich),
                       sd_min_rich = sd(min_rich),
                       #conf_low_min_rich = t.test(min_rich)$conf.int[1],
                       #conf_high_min_rich = t.test(min_rich)$conf.int[2],
                       mean_min_com_Msi0 = mean(min_com_Msi),
                       sd_min_com_Msi = sd(min_com_Msi),
                       #conf_low_min_com_Msi = t.test(min_com_Msi)$conf.int[1],
                       #conf_high_min_com_Msi = t.test(min_com_Msi)$conf.int[2],
                       mean_min_com_mass0 = mean(min_com_mass),
                       sd_min_com_mass = sd(min_com_mass),
                       #conf_low_min_com_mass = t.test(min_com_mass)$conf.int[1],
                       #conf_high_min_com_mass = t.test(min_com_mass)$conf.int[2],
                       
                       mean_med_rich0 = mean(med_rich),
                       sd_med_rich = sd(med_rich),
                       #conf_low_med_rich = t.test(med_rich)$conf.int[1],
                       #conf_high_med_rich = t.test(med_rich)$conf.int[2],
                       mean_med_com_Msi0 = mean(med_com_Msi),
                       sd_med_com_Msi = sd(med_com_Msi),
                       #conf_low_med_com_Msi = t.test(med_com_Msi)$conf.int[1],
                       #conf_high_med_com_Msi = t.test(med_com_Msi)$conf.int[2],
                       mean_med_com_mass0 = mean(med_com_mass),
                       sd_med_com_mass = sd(med_com_mass),
                       #conf_low_med_com_mass = t.test(med_com_mass)$conf.int[1],
                       #conf_high_med_com_mass = t.test(med_com_mass)$conf.int[2]
                       ) %>%
      mutate(site = i)
    
    Com_mass_mean_msi_RS_exp_all_sites <- bind_rows(Com_mass_mean_msi_RS_exp_all_sites, Com_mass_mean_msi_exp)
    
    write_xlsx(Com_mass_mean_msi_RS_exp_all_sites, 'Rich_Msi_Mass_RS_exp_all_sites_2023.xlsx')
  
}


# load data 2024 ----
setwd("D:/OneDrive/Rdata/7 spatial patterns of grassland plant species")
dot.comp_all <- read_excel("all_points_across_six_sites_2024.xlsx",sheet=1,na="NA")
head(dot.comp_all)
headTail(dot.comp_all)
str(dot.comp_all)
summary(dot.comp_all)

dot.comp_all$X_mm <- as.numeric(dot.comp_all$X_mm)
dot.comp_all$Y_mm <- as.numeric(dot.comp_all$Y_mm)

# 
i = 'WLT'
j = '2'
subset_data <- dot.comp_all[dot.comp_all$site==i&dot.comp_all$plot==j,]

subset_data <- subset_data %>%
  dplyr::select(!year) %>%
  filter(!duplicated(cbind(X_mm, Y_mm))) %>%
  mutate(label = paste(X_mm, Y_mm, sep="-"))

subset_data$species_English <- as.factor(subset_data$species_English) # make species code a factor

subset_data_ppp <- ppp(subset_data$X_mm,
                       subset_data$Y_mm,
                       xrange=c(0,2000),yrange=c(0,2000),
                       marks=subset_data$species_English)

voronoi <- deldir(subset_data_ppp)

plot(voronoi)

# 1. Spatial distribution pattern of grassland plant populations ----
## old species segregation index - relative neighborhood density Ωr ----
Ω_all_site_0_500r <- data.frame()

for (i in unique(dot.comp_all$site)) {
  for(j in unique(dot.comp_all$plot)) {
    subset_data <- dot.comp_all[dot.comp_all$site==i&dot.comp_all$plot==j,]
    
    subset_data <- subset_data %>%
      dplyr::select(!year) %>%
      filter(!duplicated(cbind(X_mm, Y_mm)))
    
    subset_data$species_English <- as.factor(subset_data$species_English) # make species code a factor
    
    Ω0_500r <- data.frame()
    bands <- seq(0, 1900, 100)
    
    for (k in 1:nrow(subset_data)) {
      for (m in bands) {
        
        spe_target <- subset_data$species_English[k]
        
        sub_Ω0_500r <- subset_data[-k,] %>%
          mutate(x_center = subset_data$X_mm[k],
                 y_center = subset_data$Y_mm[k]) %>%
          mutate(distance = sqrt((X_mm - x_center)^2 + (Y_mm - y_center)^2)) %>%
          filter(distance > m & distance <= m + 100) %>%
          mutate(spe_identity = if_else(species_English %in% spe_target, 'cD', 'hD')) %>%
          group_by(spe_identity) %>%
          dplyr::summarise(D_value = sum(abun)) %>%
          pivot_wider(names_from = spe_identity, values_from = D_value) %>%
          mutate(species_English = spe_target,
                 x_center = subset_data$X_mm[k],
                 y_center = subset_data$Y_mm[k],
                 cover_cm2 = subset_data$cover_cm2[k],
                 high_mm = subset_data$high_mm[k],
                 band = m,
                 area = case_when(
                   x_center >= m+100&x_center <= 2000-(m+100)&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2, 
                   x_center < m+100&x_center >= m&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2),
                   y_center < m+100&y_center >= m&x_center >= m+100&x_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center >= m+100&x_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2),
                   x_center < m&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   x_center > 2000-m&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   y_center < m&x_center >= m+100&x_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2),
                   y_center > 2000-m&x_center >= m+100&x_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   x_center < m+100&x_center >= m&y_center < m+100&y_center >= m&sqrt(x_center^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2),
                   x_center < m+100&x_center >= m&y_center < m+100&y_center >= m&sqrt(x_center^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 - (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 + ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2,
                   x_center < m+100&x_center >= m&y_center > 2000-(m+100)&y_center <= 2000-m&sqrt(x_center^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2),
                   x_center < m+100&x_center >= m&y_center > 2000-(m+100)&y_center <= 2000-m&sqrt(x_center^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 -(sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 + ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2,
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center < m+100&y_center >= m&sqrt((2000-x_center)^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center < m+100&y_center >= m&sqrt((2000-x_center)^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 - (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 + ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))- atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2,
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center > 2000-(m+100)&y_center <= 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center > 2000-(m+100)&y_center <= 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 - (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 + ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2,
                   x_center < m+100&x_center >= m&y_center < m&sqrt(x_center^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2),
                   x_center < m+100&x_center >= m&y_center < m&sqrt(x_center^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + + (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 + (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 - ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2 + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2), 
                   y_center < m+100&y_center >= m&x_center < m&sqrt(x_center^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   y_center < m+100&y_center >= m&x_center < m&sqrt(x_center^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + + (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 + (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 - ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2 + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   x_center < m+100&x_center >= m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   x_center < m+100&x_center >= m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + + (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 + (sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 - ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2 + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   y_center < m+100&y_center >= m&x_center > 2000-m&sqrt((2000-x_center)^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   y_center < m+100&y_center >= m&x_center > 2000-m&sqrt((2000-x_center)^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + + (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 + (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 - ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))-atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + + (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 + (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 - ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))-atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center < m&sqrt(x_center^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center < m&sqrt(x_center^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + + (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 + (sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 - ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2 + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + + (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 + (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 - ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + + (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 + (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 - ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   x_center < m&y_center < m&sqrt(x_center^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2), 
                   x_center < m&y_center < m&sqrt(x_center^2 + y_center^2) < m + 100&sqrt(x_center^2 + y_center^2) >= m ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 - (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 + ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2, 
                   x_center < m&y_center < m&sqrt(x_center^2 + y_center^2) < m ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 - (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 + ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2 + (sqrt(m^2-y_center^2)-x_center)*y_center/2 + (sqrt(m^2-x_center^2)-y_center)*x_center/2 - ((acos(y_center/m)-atan(x_center/y_center)) + (acos(x_center/m)-atan(y_center/x_center)))/(2*pi)*pi*m^2, 
                   x_center > 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2), 
                   x_center > 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) < m + 100&sqrt((2000-x_center)^2 + y_center^2) >= m ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 - (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 + ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))-atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2, 
                   x_center > 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) < m ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 - (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 + ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))-atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + (sqrt(m^2-y_center^2)-(2000-x_center))*y_center/2 + (sqrt(m^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 - ((acos(y_center/m)-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/m)-atan(y_center/(2000-x_center))))/(2*pi)*pi*m^2, 
                   x_center < m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2), 
                   x_center < m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) < m + 100&sqrt(x_center^2 + (2000-y_center)^2) >= m ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 - (sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 + ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2, 
                   x_center < m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) < m ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 - (sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 + ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2 + (sqrt(m^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 + (sqrt(m^2-x_center^2)-(2000-y_center))*x_center/2 - ((acos((2000-y_center)/m)-atan(x_center/(2000-y_center))) + (acos(x_center/m)-atan((2000-y_center)/x_center)))/(2*pi)*pi*m^2, 
                   x_center > 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2), 
                   x_center > 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m + 100&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 - (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 + ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2, 
                   x_center > 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 - (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 + ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + (sqrt(m^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 + (sqrt(m^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 - ((acos((2000-y_center)/m)-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/m)-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*m^2
                 ))
        
        Ω0_500r <- bind_rows(Ω0_500r, sub_Ω0_500r)
        
      }
      
    }    
    sel_species <- subset_data %>% 
      dplyr::select(species_English) %>%
      group_by(species_English) %>% 
      mutate(indi_num = length(species_English)) %>%
      filter(indi_num >= 20) %>%
      filter(!duplicated(species_English)) %>%
      mutate(mean_cden = indi_num/(2000*2000),
             mean_hden = (nrow(subset_data)-indi_num)/(2000*2000))
    
    Ω0_500r <- Ω0_500r %>% 
      filter(species_English %in% sel_species$species_English) %>% 
      merge(sel_species) %>%
      group_by(species_English, band) %>%
      dplyr::summarise(mean_cden = mean(mean_cden),
                       mean_hden = mean(mean_hden),
                       tcD = sum(cD),
                       thD = sum(hD),
                       tarea = sum(area),
                       spe_abun = mean(indi_num)) %>%
      mutate(rtcD = tcD/tarea,
             rthD = thD/tarea) %>%
      mutate(cΩ = rtcD/mean_cden,
             hΩ = rthD/mean_hden,
             #anpp_target = cover_cm2*high_mm, 
             site = i, plot = j)
    
    Ω_all_site_0_500r <- bind_rows(Ω_all_site_0_500r, Ω0_500r)
    
  }
  
}

write_xlsx(Ω_all_site_0_500r, 'Ω_all_site_0_500r.xlsx')

## old bootstrap method estimating confidence limits of Ωr ----
Ω_all_site_boot_0_500r <- data.frame()

for (i in unique(dot.comp_all$site)) {
  for(j in unique(dot.comp_all$plot)) {
    subset_data <- dot.comp_all[dot.comp_all$site==i&dot.comp_all$plot==j,]
    
    subset_data <- subset_data %>%
      dplyr::select(!year) %>%
      filter(!duplicated(cbind(X_mm, Y_mm)))
    
    subset_data$species_English <- as.factor(subset_data$species_English) # make species code a factor
    
    Ω0_500r <- data.frame()
    bands <- seq(0, 1900, 100)
    
    for (k in 1:nrow(subset_data)) {
      for (m in bands) {
        
        spe_target <- subset_data$species_English[k]
        
        sub_Ω0_500r <- subset_data[-k,] %>%
          mutate(x_center = subset_data$X_mm[k],
                 y_center = subset_data$Y_mm[k]) %>%
          mutate(distance = sqrt((X_mm - x_center)^2 + (Y_mm - y_center)^2)) %>%
          filter(distance > m & distance <= m + 100) %>%
          mutate(spe_identity = if_else(species_English %in% spe_target, 'cD', 'hD')) %>%
          group_by(spe_identity) %>%
          dplyr::summarise(D_value = sum(abun)) %>%
          pivot_wider(names_from = spe_identity, values_from = D_value) %>%
          mutate(species_English = spe_target,
                 x_center = subset_data$X_mm[k],
                 y_center = subset_data$Y_mm[k],
                 cover_cm2 = subset_data$cover_cm2[k],
                 high_mm = subset_data$high_mm[k],
                 band = m,
                 area = case_when(
                   x_center >= m+100&x_center <= 2000-(m+100)&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2, 
                   x_center < m+100&x_center >= m&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2),
                   y_center < m+100&y_center >= m&x_center >= m+100&x_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center >= m+100&x_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2),
                   x_center < m&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   x_center > 2000-m&y_center >= m+100&y_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   y_center < m&x_center >= m+100&x_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2),
                   y_center > 2000-m&x_center >= m+100&x_center <= 2000-(m+100) ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   x_center < m+100&x_center >= m&y_center < m+100&y_center >= m&sqrt(x_center^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2),
                   x_center < m+100&x_center >= m&y_center < m+100&y_center >= m&sqrt(x_center^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 - (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 + ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2,
                   x_center < m+100&x_center >= m&y_center > 2000-(m+100)&y_center <= 2000-m&sqrt(x_center^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2),
                   x_center < m+100&x_center >= m&y_center > 2000-(m+100)&y_center <= 2000-m&sqrt(x_center^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 -(sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 + ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2,
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center < m+100&y_center >= m&sqrt((2000-x_center)^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center < m+100&y_center >= m&sqrt((2000-x_center)^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 - (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 + ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))- atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2,
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center > 2000-(m+100)&y_center <= 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center > 2000-(m+100)&y_center <= 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 - (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 + ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2,
                   x_center < m+100&x_center >= m&y_center < m&sqrt(x_center^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2),
                   x_center < m+100&x_center >= m&y_center < m&sqrt(x_center^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + + (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 + (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 - ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2 + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2), 
                   y_center < m+100&y_center >= m&x_center < m&sqrt(x_center^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   y_center < m+100&y_center >= m&x_center < m&sqrt(x_center^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + + (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 + (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 - ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2 + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   x_center < m+100&x_center >= m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   x_center < m+100&x_center >= m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + + (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 + (sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 - ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2 + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   y_center < m+100&y_center >= m&x_center > 2000-m&sqrt((2000-x_center)^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   y_center < m+100&y_center >= m&x_center > 2000-m&sqrt((2000-x_center)^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + + (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 + (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 - ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))-atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + + (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 + (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 - ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))-atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center < m&sqrt(x_center^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center < m&sqrt(x_center^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + + (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 + (sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 - ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2 + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   x_center > 2000-(m+100)&x_center <= 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + + (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 + (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 - ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   y_center > 2000-(m+100)&y_center <= 2000-m&x_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) - acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + + (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 + (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 - ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2),
                   x_center < m&y_center < m&sqrt(x_center^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2), 
                   x_center < m&y_center < m&sqrt(x_center^2 + y_center^2) < m + 100&sqrt(x_center^2 + y_center^2) >= m ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 - (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 + ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2, 
                   x_center < m&y_center < m&sqrt(x_center^2 + y_center^2) < m ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-x_center)*y_center/2 - (sqrt((m+100)^2-x_center^2)-y_center)*x_center/2 + ((acos(y_center/(m+100))-atan(x_center/y_center)) + (acos(x_center/(m+100))-atan(y_center/x_center)))/(2*pi)*pi*(m+100)^2 + (sqrt(m^2-y_center^2)-x_center)*y_center/2 + (sqrt(m^2-x_center^2)-y_center)*x_center/2 - ((acos(y_center/m)-atan(x_center/y_center)) + (acos(x_center/m)-atan(y_center/x_center)))/(2*pi)*pi*m^2, 
                   x_center > 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2), 
                   x_center > 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) < m + 100&sqrt((2000-x_center)^2 + y_center^2) >= m ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 - (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 + ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))-atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2, 
                   x_center > 2000-m&y_center < m&sqrt((2000-x_center)^2 + y_center^2) < m ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos(y_center/(m+100))/pi*pi*(m + 100)^2 + y_center*sqrt((m + 100)^2-y_center^2) + acos(y_center/m)/pi*pi*m^2 - y_center*sqrt(m^2-y_center^2) - (sqrt((m+100)^2-y_center^2)-(2000-x_center))*y_center/2 - (sqrt((m+100)^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 + ((acos(y_center/(m+100))-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/(m+100))-atan(y_center/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + (sqrt(m^2-y_center^2)-(2000-x_center))*y_center/2 + (sqrt(m^2-(2000-x_center)^2)-y_center)*(2000-x_center)/2 - ((acos(y_center/m)-atan((2000-x_center)/y_center)) + (acos((2000-x_center)/m)-atan(y_center/(2000-x_center))))/(2*pi)*pi*m^2, 
                   x_center < m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2), 
                   x_center < m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) < m + 100&sqrt(x_center^2 + (2000-y_center)^2) >= m ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 - (sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 + ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2, 
                   x_center < m&y_center > 2000-m&sqrt(x_center^2 + (2000-y_center)^2) < m ~ pi*(m + 100)^2 - pi*m^2 - acos(x_center/(m+100))/pi*pi*(m + 100)^2 + x_center*sqrt((m + 100)^2-x_center^2) + acos(x_center/m)/pi*pi*m^2 - x_center*sqrt(m^2-x_center^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 - (sqrt((m+100)^2-x_center^2)-(2000-y_center))*x_center/2 + ((acos((2000-y_center)/(m+100))-atan(x_center/(2000-y_center))) + (acos(x_center/(m+100))-atan((2000-y_center)/x_center)))/(2*pi)*pi*(m+100)^2 + (sqrt(m^2-(2000-y_center)^2)-x_center)*(2000-y_center)/2 + (sqrt(m^2-x_center^2)-(2000-y_center))*x_center/2 - ((acos((2000-y_center)/m)-atan(x_center/(2000-y_center))) + (acos(x_center/m)-atan((2000-y_center)/x_center)))/(2*pi)*pi*m^2, 
                   x_center > 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m + 100 ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2), 
                   x_center > 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m + 100&sqrt((2000-x_center)^2 + (2000-y_center)^2) >= m ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 - (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 + ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2, 
                   x_center > 2000-m&y_center > 2000-m&sqrt((2000-x_center)^2 + (2000-y_center)^2) < m ~ pi*(m + 100)^2 - pi*m^2 - acos((2000-x_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-x_center)*sqrt((m + 100)^2-(2000-x_center)^2) + acos((2000-x_center)/m)/pi*pi*m^2 - (2000-x_center)*sqrt(m^2-(2000-x_center)^2)- acos((2000-y_center)/(m+100))/pi*pi*(m + 100)^2 + (2000-y_center)*sqrt((m + 100)^2-(2000-y_center)^2) + acos((2000-y_center)/m)/pi*pi*m^2 - (2000-y_center)*sqrt(m^2-(2000-y_center)^2) - (sqrt((m+100)^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 - (sqrt((m+100)^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 + ((acos((2000-y_center)/(m+100))-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/(m+100))-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*(m+100)^2 + (sqrt(m^2-(2000-y_center)^2)-(2000-x_center))*(2000-y_center)/2 + (sqrt(m^2-(2000-x_center)^2)-(2000-y_center))*(2000-x_center)/2 - ((acos((2000-y_center)/m)-atan((2000-x_center)/(2000-y_center))) + (acos((2000-x_center)/m)-atan((2000-y_center)/(2000-x_center))))/(2*pi)*pi*m^2
                 ))
        
        Ω0_500r <- bind_rows(Ω0_500r, sub_Ω0_500r)
        
      }
      
    } 
    
    sel_species <- subset_data %>% 
      dplyr::select(species_English) %>%
      group_by(species_English) %>% 
      mutate(indi_num = length(species_English)) %>%
      filter(indi_num >= 20) %>%
      filter(!duplicated(species_English)) %>%
      mutate(mean_cden = indi_num/(2000*2000),
             mean_hden = (nrow(subset_data)-indi_num)/(2000*2000))
    
    Ω0_500r <- Ω0_500r %>% 
      filter(species_English %in% sel_species$species_English) %>% 
      merge(sel_species) 
    
    Ω0_500r_boot <- data.frame()
    
    for (w in unique(sel_species$species_English)) {
      for(q in 1:15) {
        
        Ω0_500r_spe <- Ω0_500r[Ω0_500r$species_English==w,]
        
        n <- as.integer(nrow(Ω0_500r_spe)/2)
        
        random_rows <- sample(1:nrow(Ω0_500r_spe), n)
        Ω0_500r_sampled_data <- Ω0_500r_spe[random_rows, ]
        
        Ω0_boot_500r_spe <- Ω0_500r_sampled_data %>%
          group_by(band) %>%
          dplyr::summarise(mean_cden = mean(mean_cden),
                           mean_hden = mean(mean_hden),
                           tcD = sum(cD),
                           thD = sum(hD),
                           tarea = sum(area),
                           spe_abun = mean(indi_num)) %>%
          mutate(rtcD = tcD/tarea,
                 rthD = thD/tarea) %>%
          mutate(cΩ = rtcD/mean_cden,
                 hΩ = rthD/mean_hden,
                 #anpp_target = cover_cm2*high_mm, 
                 time = q)
        
        Ω0_500r_boot <- bind_rows(Ω0_500r_boot, Ω0_boot_500r_spe)
        
      }
    }
    
    Ω0_500r_boot <- Ω0_500r_boot %>%
      group_by(species_English, band) %>%
      dplyr::summarise(mean_cΩ = mean(cΩ),
                       #sd_cΩ = sd(cΩ),
                       conf_low_cΩ = t.test(cΩ)$conf.int[1]/2,
                       conf_high_cΩ = t.test(cΩ)$conf.int[2]/2,
                       mean_hΩ = mean(hΩ),
                       #sd_hΩ = sd(hΩ),
                       conf_low_hΩ = t.test(hΩ)$conf.int[1]/2,
                       conf_high_hΩ = t.test(hΩ)$conf.int[2]/2,
                       spe_abun = mean(indi_num))
    mutate(site = i, plot = j)
    
    Ω_all_site_boot_0_500r <- bind_rows(Ω_all_site_boot_0_500r, Ω0_500r_boot)
    
  }
  
}

write_xlsx(Ω_all_site_boot_0_500r, 'Ω_all_site_boot_0_500r.xlsx')


## species segregation index Ms ----
Ms_all_site <- data.frame() 
Com_mass_mean_msi_all_site_alpha <- data.frame()
Com_mass_mean_msi_all_site_beta <- data.frame()
scale_alpha <- c(100,150,200,250,300,350,400,450,500,600,700,800,900,1000)
scale_beta <- c(100,125,150,175,200,225,250,300,350,400,450,500)

for (i in unique(dot.comp_all$site)) {
  for(j in unique(dot.comp_all$plot)) {
    subset_data <- dot.comp_all[dot.comp_all$site==i&dot.comp_all$plot==j,]
    
    subset_data <- subset_data %>%
      dplyr::select(!year) %>%
      filter(!duplicated(cbind(X_mm, Y_mm))) %>%
      mutate(label = paste(X_mm, Y_mm, sep="-"))
    
    subset_data$species_English <- as.factor(subset_data$species_English) # make species code a factor
    
    subset_data_ppp <- ppp(subset_data$X_mm,
                           subset_data$Y_mm,
                           xrange=c(0,2000),yrange=c(0,2000),
                           marks=subset_data$species_English)
    
    voronoi <- deldir(subset_data_ppp)
    
    #write_xlsx(voronoi$summary, 'species segregation index_polygons.xlsx')
    #write_xlsx(voronoi$delsgs, 'species segregation index_delsgs.xlsx')
    #write_xlsx(voronoi$dirsgs, 'species segregation index_dirsgs.xlsx')
    
    delsgs_a <- voronoi$delsgs %>%
      dplyr::select(x_center = x1, y_center = y1, x_neigh = x2, y_neigh = y2, id_center = ind1)
    delsgs_b <- voronoi$delsgs %>%
      dplyr::select(x_center = x2, y_center = y2, x_neigh = x1, y_neigh = y1, id_center = ind2)
    
    delsgs <- bind_rows(delsgs_a, delsgs_b)
    
    #write_xlsx(delsgs, 'species segregation index_delsgs_all.xlsx')
    
    delsgs$label = paste(delsgs$x_neigh, delsgs$y_neigh, sep="-")
    #delsgs$labelcenter = paste(delsgs$x_center, delsgs$y_center, sep="-")
    #length(unique(delsgs$labelcenter)) # = nrow(subset_data) ensure all individuals as focal individuals
    
    neigh <- delsgs %>% 
      merge(subset_data %>% dplyr::select(label, species_English_neigh = species_English, cover_cm2_neigh = cover_cm2, high_mm_neigh = high_mm, abun_neigh = abun)) %>%
      dplyr::select(X_mm = x_center, Y_mm = y_center, id_center, x_neigh, y_neigh, species_English_neigh, cover_cm2_neigh, high_mm_neigh, abun_neigh) %>%
      merge(subset_data)
    #length(unique(neigh$id_center))   
    
    out_plant_id <- voronoi$dirsgs %>%
      mutate(outn = thirdv1*thirdv2) %>%
      filter(outn < 0) %>%
      dplyr::select(ind1, ind2) %>%
      pivot_longer(1:2, names_to = 'ind', values_to = 'id_center') %>%
      filter(!duplicated(id_center))
    
    neigh_sel <- neigh %>%
      mutate(outy = if_else(id_center %in% out_plant_id$id_center, 0, 1)) %>%
      filter(outy == 1) %>%
      mutate(spe_identity = if_else(species_English_neigh == species_English, '0', '1')) %>%
      mutate(spe_identity = as.numeric(spe_identity))
    
    sub_Msi <- neigh_sel %>%
      group_by(X_mm, Y_mm, species_English) %>%
      dplyr::summarise(Msi = if_else(sum(spe_identity) == length(species_English_neigh), 
                                     sum(spe_identity)*(length(unique(species_English_neigh))+1)/length(species_English_neigh)/(length(species_English_neigh)+1),
                                     sum(spe_identity)*length(unique(species_English_neigh))/length(species_English_neigh)/(length(species_English_neigh)+1)),
                       neigh_abun = length(species_English_neigh),
                       neigh_rich = length(unique(species_English_neigh)),
                       focal_Ind_mass = log(mean(cover_cm2*high_mm)+1)) 
    
    com_Msi = mean(sub_Msi$Msi)
    
    sub_Ms_all_site <- sub_Msi %>% merge(neigh_sel) %>% mutate(site = i, plot = j, com_Msi = com_Msi) 
    
    Ms_all_site <- bind_rows(Ms_all_site, sub_Ms_all_site)
    write_xlsx(Ms_all_site, 'Ms_all_site_2024.xlsx')
    
    Com_mass_mean_msi_alpha <- data.frame()
    
    for (k in scale_alpha) {
      sub_Ms_all_site$xt <- cut(sub_Ms_all_site$X_mm,seq(min(sub_Ms_all_site$X_mm),round(max(sub_Ms_all_site$X_mm),min(sub_Ms_all_site$X_mm)),k)) 
      sub_Ms_all_site$yt <- cut(sub_Ms_all_site$Y_mm,seq(min(sub_Ms_all_site$Y_mm),round(max(sub_Ms_all_site$Y_mm),min(sub_Ms_all_site$Y_mm)),k))
      
      sub_Com_mass_mean_msi <- sub_Ms_all_site %>%
        filter(xt != 'NA') %>%
        filter(yt != 'NA') %>%
        group_by(site, plot, xt, yt, X_mm, Y_mm, species_English) %>%
        dplyr::summarise(Msi=mean(Msi),
                         abun=mean(abun),
                         neigh_abun1=mean(neigh_abun),neigh_abun2=sum(abun_neigh),
                         heter_neigh_abun1=sum(spe_identity),heter_neigh_abun2=sum(abun_neigh*spe_identity),
                         neigh_rich=mean(neigh_rich),cover_cm2=mean(cover_cm2),high_mm=mean(high_mm)) %>%
        mutate(con_neigh_abun1=neigh_abun1-heter_neigh_abun1,con_neigh_abun2=neigh_abun2-heter_neigh_abun2) %>%
        ungroup() %>%
        group_by(site, plot, xt, yt) %>%
        dplyr::summarise(rich = length(unique(species_English)),
                         com_mass=log(sum(cover_cm2*high_mm)+1),
                         com_Msi=mean(Msi),
                         mean_neigh_abun1=mean(neigh_abun1),mean_neigh_abun2=mean(neigh_abun2),
                         mean_heter_neigh_abun1=mean(heter_neigh_abun1),mean_heter_neigh_abun2=mean(heter_neigh_abun2),
                         mean_con_neigh_abun1=mean(con_neigh_abun1),mean_con_neigh_abun2=mean(con_neigh_abun2),
                         mean_neigh_rich=mean(neigh_rich),com_abun=sum(abun)) %>%
        mutate(scale = k)
      
      Com_mass_mean_msi_alpha <- bind_rows(Com_mass_mean_msi_alpha, sub_Com_mass_mean_msi)
      
    }
    
    Com_mass_mean_msi_all_site_alpha <- bind_rows(Com_mass_mean_msi_all_site_alpha,Com_mass_mean_msi_alpha)
    
    write_xlsx(Com_mass_mean_msi_all_site_alpha, 'Com_mass_mean_msi_all_site_alpha_2024.xlsx')
    
    Com_mass_mean_msi_beta <- data.frame()
    
    for (s in scale_beta) {
      
      sub_Ms_all_site$xt1 <- cut(sub_Ms_all_site$X_mm,seq(min(sub_Ms_all_site$X_mm),round(max(sub_Ms_all_site$X_mm),min(sub_Ms_all_site$X_mm)),s)) 
      sub_Ms_all_site$yt1 <- cut(sub_Ms_all_site$Y_mm,seq(min(sub_Ms_all_site$Y_mm),round(max(sub_Ms_all_site$Y_mm),min(sub_Ms_all_site$Y_mm)),s))
      
      sub_Ms_all_site$xt2 <- cut(sub_Ms_all_site$X_mm,seq(min(sub_Ms_all_site$X_mm),round(max(sub_Ms_all_site$X_mm),min(sub_Ms_all_site$X_mm)),2*s)) 
      sub_Ms_all_site$yt2 <- cut(sub_Ms_all_site$Y_mm,seq(min(sub_Ms_all_site$Y_mm),round(max(sub_Ms_all_site$Y_mm),min(sub_Ms_all_site$Y_mm)),2*s))
      
      
      sub_Ms_all_site2 <- sub_Ms_all_site %>%
        filter(xt2 != 'NA') %>%
        filter(yt2 != 'NA')
      
      rich1 <- sub_Ms_all_site2 %>%
        group_by(site, plot, xt2, yt2, xt1, yt1) %>%
        dplyr::summarise(rich = length(unique(species_English))) %>%
        group_by(site, plot, xt2, yt2) %>%
        dplyr::summarise(rich1 = mean(rich))
      
      rich2 <- sub_Ms_all_site2 %>%
        group_by(site, plot, xt2, yt2) %>%
        dplyr::summarise(rich2 = length(unique(species_English)))
      
      beta_W <- rich1 %>%
        merge(rich2) %>%
        mutate(beta_add = rich2-rich1,
               beta_mul = rich2/rich1)
      
      spe_abu <- sub_Ms_all_site2 %>% 
        group_by(xt2, yt2, xt1, yt1, X_mm, Y_mm, species_English) %>%
        dplyr::summarise(abun = mean(abun)) %>%
        ungroup() %>%
        dplyr::select(xt2, yt2, xt1, yt1, species_English, abun) %>% 
        group_by(xt2, yt2, xt1, yt1, species_English) %>% 
        dplyr::summarise(abun = sum(abun)) %>%
        mutate(xyt2 = paste(xt2, yt2, sep="-"))
      
      beta_jac_bray <- data.frame()
      
      for (m in unique(spe_abu$xyt2)) {
        subdata <- spe_abu[spe_abu$xyt2 == m, ]
        
        subdata_jac <- subdata %>%
          ungroup() %>%
          dplyr::select(xt1, yt1, species_English, abun) %>% 
          mutate(abun = ifelse(abun > 0, 1, 0)) %>%
          pivot_wider(values_from = abun, names_from = species_English)
        
        subdata_jac <- replace(subdata_jac, is.na(subdata_jac), 0)
        sub_beta_jac <- data.frame(beta.multi(subdata_jac[, -c(1:2)], index.family = 'jaccard'))
        sub_beta_jac <- sub_beta_jac %>% mutate(xyt2 = m)
        
        subdata_bray <- subdata %>%
          ungroup() %>%
          dplyr::select(xt1, yt1, species_English, abun) %>% 
          pivot_wider(values_from = abun, names_from = species_English)
        
        subdata_bray <- replace(subdata_bray, is.na(subdata_bray), 0)
        sub_beta_bray <- data.frame(beta.multi.abund(subdata_bray[, -c(1:2)], index.family = 'bray'))
        sub_beta_bray <- sub_beta_bray %>% mutate(xyt2 = m)
        
        sub_beta_jac_bray <- sub_beta_jac %>% merge(sub_beta_bray)
        
        sub_beta_jac_bray <- separate(sub_beta_jac_bray,xyt2,into = c("xt2","yt2"),sep = "-")
        
        
        beta_jac_bray <- bind_rows(beta_jac_bray, sub_beta_jac_bray)
        
      }
      
      sub_Com_mass_mean_msi2 <- sub_Ms_all_site %>%
        filter(xt2 != 'NA') %>%
        filter(yt2 != 'NA') %>%
        group_by(site, plot, xt2, yt2, X_mm, Y_mm, species_English) %>%
        dplyr::summarise(Msi=mean(Msi),
                         abun=mean(abun),
                         neigh_abun1=mean(neigh_abun),neigh_abun2=sum(abun_neigh),
                         heter_neigh_abun1=sum(spe_identity),heter_neigh_abun2=sum(abun_neigh*spe_identity),
                         neigh_rich=mean(neigh_rich),cover_cm2=mean(cover_cm2),high_mm=mean(high_mm)) %>%
        mutate(con_neigh_abun1=neigh_abun1-heter_neigh_abun1,con_neigh_abun2=neigh_abun2-heter_neigh_abun2) %>%
        ungroup() %>%
        group_by(site, plot, xt2, yt2) %>%
        dplyr::summarise(rich = length(unique(species_English)),
                         com_mass=log(sum(cover_cm2*high_mm)+1),
                         com_Msi=mean(Msi),
                         mean_neigh_abun1=mean(neigh_abun1),mean_neigh_abun2=mean(neigh_abun2),
                         mean_heter_neigh_abun1=mean(heter_neigh_abun1),mean_heter_neigh_abun2=mean(heter_neigh_abun2),
                         mean_con_neigh_abun1=mean(con_neigh_abun1),mean_con_neigh_abun2=mean(con_neigh_abun2),
                         mean_neigh_rich=mean(neigh_rich),com_abun=sum(abun)) %>%
        mutate(scale = 2*s) %>%
        merge(beta_W) %>% merge(beta_jac_bray)
      
      Com_mass_mean_msi_beta <- bind_rows(Com_mass_mean_msi_beta, sub_Com_mass_mean_msi2)
      
    }
    
    Com_mass_mean_msi_all_site_beta <- bind_rows(Com_mass_mean_msi_all_site_beta,Com_mass_mean_msi_beta)
    write_xlsx(Com_mass_mean_msi_all_site_beta, 'Com_mass_mean_msi_all_site_beta_2024.xlsx')
  }
}

## expected value based on CSR of species segregation index Ms ----
Com_mass_mean_msi_RS_exp_all_sites <- data.frame()
scale <- c(100,150,200,250,300,350,400,450,500,600,700,800,900,1000)

for (i in unique(dot.comp_all$site)) {
  for(j in unique(dot.comp_all$plot)) {
    subset_data <- dot.comp_all[dot.comp_all$site==i&dot.comp_all$plot==j,]
    
    subset_data <- subset_data %>%
      dplyr::select(!year) %>%
      filter(!duplicated(cbind(X_mm, Y_mm)))
    
    subset_data <- subset_data %>%
      mutate(id = 1:nrow(subset_data))
    
    subset_data$species_English <- as.factor(subset_data$species_English) # make species code a factor
    
    Com_mass_mean_msi_all_time <- data.frame()
    
    for (m in 1:199) {
      #for (w in scale) {
      
      n = nrow(subset_data)
      
      random_spe <- subset_data %>%
        dplyr::select(species_English, cover_cm2, high_mm, abun) %>%
        mutate(id = sample(1:n, n))
      
      subset_data2 <- subset_data %>%
        dplyr::select(-c(species_English, cover_cm2, high_mm, abun))
      
      subset_data_csr <- subset_data2 %>%
        merge(random_spe) %>%
        dplyr::select(-id) %>%
        mutate(label = paste(X_mm, Y_mm, sep="-"))
      
      subset_data_csr_ppp <- ppp(subset_data_csr$X_mm,
                                 subset_data_csr$Y_mm,
                                 xrange=c(0,2000),yrange=c(0,2000),
                                 marks=subset_data_csr$species_English)
      
      voronoi <- deldir(subset_data_csr_ppp)
      
      delsgs_a <- voronoi$delsgs %>%
        dplyr::select(x_center = x1, y_center = y1, x_neigh = x2, y_neigh = y2, id_center = ind1)
      delsgs_b <- voronoi$delsgs %>%
        dplyr::select(x_center = x2, y_center = y2, x_neigh = x1, y_neigh = y1, id_center = ind2)
      
      delsgs <- bind_rows(delsgs_a, delsgs_b)
      
      delsgs$label = paste(delsgs$x_neigh, delsgs$y_neigh, sep="-")
      #delsgs$labelcenter = paste(delsgs$x_center, delsgs$y_center, sep="-")
      #length(unique(delsgs$labelcenter)) # = nrow(subset_data) ensure all individuals as focal individuals
      
      neigh <- delsgs %>% 
        merge(subset_data_csr %>% dplyr::select(label, species_English_neigh = species_English, cover_cm2_neigh = cover_cm2, high_mm_neigh = high_mm , abun_neigh = abun)) %>%
        dplyr::select(X_mm = x_center, Y_mm = y_center, id_center, x_neigh, y_neigh, species_English_neigh, cover_cm2_neigh, high_mm_neigh, abun_neigh) %>%
        merge(subset_data_csr)
      #length(unique(neigh$id_center))   
      
      out_plant_id <- voronoi$dirsgs %>%
        mutate(outn = thirdv1*thirdv2) %>%
        filter(outn < 0) %>%
        dplyr::select(ind1, ind2) %>%
        pivot_longer(1:2, names_to = 'ind', values_to = 'id_center') %>%
        filter(!duplicated(id_center))
      
      neigh_sel <- neigh %>%
        mutate(outy = if_else(id_center %in% out_plant_id$id_center, 0, 1)) %>%
        filter(outy == 1) %>%
        mutate(spe_identity = if_else(species_English_neigh == species_English, '0', '1')) %>%
        mutate(spe_identity = as.numeric(spe_identity))
      
      sub_Msi <- neigh_sel %>%
        group_by(X_mm, Y_mm, species_English) %>%
        dplyr::summarise(Msi = if_else(sum(spe_identity) == length(species_English_neigh), 
                                       sum(spe_identity)*(length(unique(species_English_neigh))+1)/length(species_English_neigh)/(length(species_English_neigh)+1),
                                       sum(spe_identity)*length(unique(species_English_neigh))/length(species_English_neigh)/(length(species_English_neigh)+1)),
                         neigh_abun = length(species_English_neigh),
                         neigh_rich = length(unique(species_English_neigh))) 
      
      com_Msi = mean(sub_Msi$Msi)
      
      sub_Ms_all_site <- sub_Msi %>% merge(neigh_sel) #%>% mutate(site = i, plot = j, com_Msi = com_Msi) 
      
      #Ms_all_site <- bind_rows(Ms_all_site, sub_Ms_all_site)
      
      Com_mass_mean_msi <- data.frame()
      
      for (k in scale) {
        sub_Ms_all_site$xt <- cut(sub_Ms_all_site$X_mm,seq(min(sub_Ms_all_site$X_mm),round(max(sub_Ms_all_site$X_mm),min(sub_Ms_all_site$X_mm)),k)) 
        sub_Ms_all_site$yt <- cut(sub_Ms_all_site$Y_mm,seq(min(sub_Ms_all_site$Y_mm),round(max(sub_Ms_all_site$Y_mm),min(sub_Ms_all_site$Y_mm)),k))
        
        sub_Com_mass_mean_msi <- sub_Ms_all_site %>%
          filter(xt != 'NA') %>%
          filter(yt != 'NA') %>%
          dplyr::group_by(xt, yt, X_mm, Y_mm, species_English) %>%
          dplyr::summarise(Msi=mean(Msi),
                           neigh_abun1=mean(neigh_abun),neigh_abun2=sum(abun_neigh),
                           heter_neigh_abun1=sum(spe_identity),heter_neigh_abun2=sum(abun_neigh*spe_identity),
                           neigh_rich=mean(neigh_rich),cover_cm2=mean(cover_cm2),high_mm=mean(high_mm)) %>%
          dplyr::mutate(con_neigh_abun1=neigh_abun1-heter_neigh_abun1,con_neigh_abun2=neigh_abun2-heter_neigh_abun2) %>%
          dplyr::mutate(struc_rich=if_else(con_neigh_abun1 == 0, neigh_rich+1, neigh_rich)) %>%
          dplyr::ungroup() %>%
          dplyr::group_by(xt, yt) %>%
          dplyr::summarise(rich = length(unique(species_English)),
                           com_mass=log(sum(cover_cm2*high_mm)+1),
                           com_Msi=mean(Msi),
                           mean_neigh_abun1=mean(neigh_abun1),mean_neigh_abun2=mean(neigh_abun2),
                           mean_heter_neigh_abun1=mean(heter_neigh_abun1),mean_heter_neigh_abun2=mean(heter_neigh_abun2),
                           mean_con_neigh_abun1=mean(con_neigh_abun1),mean_con_neigh_abun2=mean(con_neigh_abun2),
                           mean_neigh_rich=mean(neigh_rich),
                           mean_struc_rich=mean(struc_rich)) %>%
          dplyr::mutate(scale = k)
        
        Com_mass_mean_msi <- bind_rows(Com_mass_mean_msi, sub_Com_mass_mean_msi)
      }
      
      Com_mass_mean_msi <- Com_mass_mean_msi %>% mutate(time = m)
      
      Com_mass_mean_msi_all_time <- bind_rows(Com_mass_mean_msi_all_time,Com_mass_mean_msi)
      
    }
    
    Com_mass_mean_msi_exp <- Com_mass_mean_msi_all_time %>%
      group_by(xt, yt, scale) %>%
      dplyr::summarise(mean_rich = mean(rich),
                       sd_rich = sd(rich),
                       #conf_low_rich = t.test(rich)$conf.int[1],
                       #conf_high_rich = t.test(rich)$conf.int[2],
                       mean_com_mass = mean(com_mass),
                       sd_com_mass = sd(com_mass),
                       conf_low_com_mass = t.test(com_mass)$conf.int[1],
                       conf_high_com_mass = t.test(com_mass)$conf.int[2],
                       mean_com_Msi = mean(com_Msi),
                       sd_com_Msi = sd(com_Msi),
                       conf_low_com_Msi = t.test(com_Msi)$conf.int[1],
                       conf_high_com_Msi = t.test(com_Msi)$conf.int[2],
                       mean_neigh_abun1 = mean(mean_neigh_abun1),
                       sd_neigh_abun1 = sd(mean_neigh_abun1),
                       #conf_low_neigh_abun1 = t.test(mean_neigh_abun1)$conf.int[1],
                       #conf_high_neigh_abun1 = t.test(mean_neigh_abun1)$conf.int[2],
                       aver_neigh_abun2 = mean(mean_neigh_abun2),
                       sd_neigh_abun2 = sd(mean_neigh_abun2),
                       conf_low_neigh_abun2 = t.test(mean_neigh_abun2)$conf.int[1],
                       conf_high_neigh_abun2 = t.test(mean_neigh_abun2)$conf.int[2],
                       aver_heter_neigh_abun1 = mean(mean_heter_neigh_abun1),
                       sd_heter_neigh_abun1 = sd(mean_heter_neigh_abun1),
                       conf_low_heter_neigh_abun1 = t.test(mean_heter_neigh_abun1)$conf.int[1],
                       conf_high_heter_neigh_abun1 = t.test(mean_heter_neigh_abun1)$conf.int[2],
                       aver_heter_neigh_abun2 = mean(mean_heter_neigh_abun2),
                       sd_heter_neigh_abun2 = sd(mean_heter_neigh_abun2),
                       conf_low_heter_neigh_abun2 = t.test(mean_heter_neigh_abun2)$conf.int[1],
                       conf_high_heter_neigh_abun2 = t.test(mean_heter_neigh_abun2)$conf.int[2],
                       aver_con_neigh_abun1 = mean(mean_con_neigh_abun1),
                       sd_con_neigh_abun1 = sd(mean_con_neigh_abun1),
                       conf_low_con_neigh_abun1 = t.test(mean_con_neigh_abun1)$conf.int[1],
                       conf_high_con_neigh_abun1 = t.test(mean_con_neigh_abun1)$conf.int[2],
                       aver_con_neigh_abun2 = mean(mean_con_neigh_abun2),
                       sd_con_neigh_abun2 = sd(mean_con_neigh_abun2),
                       conf_low_con_neigh_abun2 = t.test(mean_con_neigh_abun2)$conf.int[1],
                       conf_high_con_neigh_abun2 = t.test(mean_con_neigh_abun2)$conf.int[2],
                       aver_neigh_rich = mean(mean_neigh_rich),
                       sd_neigh_rich = sd(mean_neigh_rich),
                       conf_low_neigh_rich = t.test(mean_neigh_rich)$conf.int[1],
                       conf_high_neigh_rich = t.test(mean_neigh_rich)$conf.int[2],
                       aver_struc_rich = mean(mean_struc_rich),
                       sd_struc_rich = sd(mean_struc_rich),
                       conf_low_struc_rich = t.test(mean_struc_rich)$conf.int[1],
                       conf_high_struc_rich = t.test(mean_struc_rich)$conf.int[2]) %>%
      mutate(site = i, plot = j)
    
    Com_mass_mean_msi_RS_exp_all_sites <- bind_rows(Com_mass_mean_msi_RS_exp_all_sites, Com_mass_mean_msi_exp)
    
    write_xlsx(Com_mass_mean_msi_RS_exp_all_sites, 'Com_mass_mean_msi_RS_exp_all_sites_2024.xlsx')
    
  }
}

## species segregation index Ms distributed ----
sel_species <- read_excel("sel_spe_2024.xlsx",sheet=1,na="NA")
Ms_distribution_all_site <- data.frame()

for (i in unique(dot.comp_all$site)) {
  for(j in unique(dot.comp_all$plot)) {
    subset_data <- dot.comp_all[dot.comp_all$site==i&dot.comp_all$plot==j,]
    
    subset_data <- subset_data %>%
      dplyr::select(!year) %>%
      filter(!duplicated(cbind(X_mm, Y_mm)))
    
    subset_data <- subset_data %>%
      mutate(id = 1:nrow(subset_data))
    
    subset_data$species_English <- as.factor(subset_data$species_English) # make species code a factor
    
    sub_Msi_sel <- sel_species[sel_species$site==i&sel_species$plot==j,]
    
    Ms_distribution_plot <- data.frame()
    
    for (m in 1:99) {
      #for (w in scale) {
      
      n = nrow(subset_data)
      
      random_spe <- subset_data %>%
        dplyr::select(species_English, cover_cm2, high_mm, abun) %>%
        mutate(id = sample(1:n, n))
      
      subset_data2 <- subset_data %>%
        dplyr::select(-c(species_English, cover_cm2, high_mm, abun))
      
      subset_data_csr <- subset_data2 %>%
        merge(random_spe) %>%
        dplyr::select(-id) %>%
        mutate(label = paste(X_mm, Y_mm, sep="-"))
      
      subset_data_csr_ppp <- ppp(subset_data_csr$X_mm,
                                 subset_data_csr$Y_mm,
                                 xrange=c(0,2000),yrange=c(0,2000),
                                 marks=subset_data_csr$species_English)
      
      voronoi <- deldir(subset_data_csr_ppp)
      
      delsgs_a <- voronoi$delsgs %>%
        dplyr::select(x_center = x1, y_center = y1, x_neigh = x2, y_neigh = y2, id_center = ind1)
      delsgs_b <- voronoi$delsgs %>%
        dplyr::select(x_center = x2, y_center = y2, x_neigh = x1, y_neigh = y1, id_center = ind2)
      
      delsgs <- bind_rows(delsgs_a, delsgs_b)
      
      delsgs$label = paste(delsgs$x_neigh, delsgs$y_neigh, sep="-")
      #delsgs$labelcenter = paste(delsgs$x_center, delsgs$y_center, sep="-")
      #length(unique(delsgs$labelcenter)) # = nrow(subset_data) ensure all individuals as focal individuals
      
      neigh <- delsgs %>% 
        merge(subset_data_csr %>% dplyr::select(label, species_English_neigh = species_English, cover_cm2_neigh = cover_cm2, high_mm_neigh = high_mm , abun_neigh = abun)) %>%
        dplyr::select(X_mm = x_center, Y_mm = y_center, id_center, x_neigh, y_neigh, species_English_neigh, cover_cm2_neigh, high_mm_neigh, abun_neigh) %>%
        merge(subset_data_csr)
      #length(unique(neigh$id_center))   
      
      out_plant_id <- voronoi$dirsgs %>%
        mutate(outn = thirdv1*thirdv2) %>%
        filter(outn < 0) %>%
        dplyr::select(ind1, ind2) %>%
        pivot_longer(1:2, names_to = 'ind', values_to = 'id_center') %>%
        filter(!duplicated(id_center))
      
      neigh_sel <- neigh %>%
        mutate(outy = if_else(id_center %in% out_plant_id$id_center, 0, 1)) %>%
        filter(outy == 1) %>%
        mutate(spe_identity = if_else(species_English_neigh == species_English, '0', '1')) %>%
        mutate(spe_identity = as.numeric(spe_identity))
      
      sub_Msi <- neigh_sel %>%
        group_by(X_mm, Y_mm, species_English) %>%
        dplyr::summarise(Msi = if_else(sum(spe_identity) == length(species_English_neigh), 
                                       sum(spe_identity)*(length(unique(species_English_neigh))+1)/length(species_English_neigh)/(length(species_English_neigh)+1),
                                       sum(spe_identity)*length(unique(species_English_neigh))/length(species_English_neigh)/(length(species_English_neigh)+1)),
                         neigh_abun = length(species_English_neigh),
                         neigh_rich = length(unique(species_English_neigh))) 
      
      #com_Msi = mean(sub_Msi$Msi)
      #sub_Ms_all_site <- sub_Msi %>% merge(neigh_sel)
      
      sub_Msi <- sub_Msi %>%
        filter(species_English %in% sub_Msi_sel$species_English) %>%
        mutate(time = as.factor(m)) 
      
      Ms_distribution_plot <- bind_rows(Ms_distribution_plot, sub_Msi)
    }
    
    sub_Ms_distribution <- data.frame()
    
    for (s in unique(sub_Msi_sel$species_English)) {
      
      sub_Ms_distribution_plot <- Ms_distribution_plot[Ms_distribution_plot$species_English==s,]
      
      global_min <- min(sub_Ms_distribution_plot$Msi)
      global_max <- max(sub_Ms_distribution_plot$Msi)
      
      x_grid <- seq(global_min, global_max, length.out = 100)
      
      sub_Ms_distribution_plot <- sub_Ms_distribution_plot %>% dplyr::select(time, Msi)
      
      list_sub_Ms_distribution_plot <- split(sub_Ms_distribution_plot$Msi, sub_Ms_distribution_plot$time)
      
      densities_spe <- lapply(list_sub_Ms_distribution_plot, function(sim) {
        dens <- density(sim, n = 100)
        approx(dens$x, dens$y, xout = x_grid)$y
      })
      
      density_matrix_spe <- do.call(rbind, densities_spe)
      density_matrix_spe[is.na(density_matrix_spe)] <- 0
      
      mean_density <- colMeans(density_matrix_spe)
      lower_ci <- apply(density_matrix_spe, 2, quantile, probs = 0.025)
      upper_ci <- apply(density_matrix_spe, 2, quantile, probs = 0.975)
      
      sub_Ms_distribution_spe <- data.frame(x = x_grid, mean = mean_density, lower = lower_ci, upper = upper_ci, 
                                            species_English = s)
      
      sub_Ms_distribution <- bind_rows(sub_Ms_distribution, sub_Ms_distribution_spe)
    }
    
    sub_Ms_distribution <- sub_Ms_distribution %>% mutate(site = i, plot = j)
    
    Ms_distribution_all_site <- bind_rows(Ms_distribution_all_site, sub_Ms_distribution)
    
    write_xlsx(Ms_distribution_all_site, 'Ms_distribution_all_site_2024.xlsx')
  }
}

## maybe for species --- expected value based on CSR of species segregation index Ms ----
Ms_all_CSR_exp_site <- data.frame() 
#Ms_all_CSR_exp_com_scaling <- data.frame()
scale <- c(100,150,200,250,300,350,400,450,500,600,700,800,900,1000)

for (i in unique(dot.comp_all$site)) {
  for(j in unique(dot.comp_all$plot)) {
    subset_data <- dot.comp_all[dot.comp_all$site==i&dot.comp_all$plot==j,]
    
    subset_data <- subset_data %>%
      dplyr::select(!year) %>%
      filter(!duplicated(cbind(X_mm, Y_mm)))
    
    subset_data <- subset_data %>%
      mutate(id = 1:nrow(subset_data))
    
    subset_data$species_English <- as.factor(subset_data$species_English) # make species code a factor
    
    #sub_Msi_sel <- sel_species[sel_species$site==i&sel_species$plot==j,]
    
    sub_Ms_csr <- data.frame()
    
    for (k in 1:99) {
      #for (w in scale) {
      
      n = nrow(subset_data)
      
      random_spe <- subset_data %>%
        dplyr::select(species_English, cover_cm2, high_mm, abun) %>%
        mutate(id = sample(1:n, n))
      
      subset_data2 <- subset_data %>%
        dplyr::select(-c(species_English, cover_cm2, high_mm, abun))
      
      subset_data_csr <- subset_data2 %>%
        merge(random_spe) %>%
        dplyr::select(-id) %>%
        mutate(label = paste(X_mm, Y_mm, sep="-"))
      
      subset_data_csr_ppp <- ppp(subset_data_csr$X_mm,
                                 subset_data_csr$Y_mm,
                                 xrange=c(0,2000),yrange=c(0,2000),
                                 marks=subset_data_csr$species_English)
      
      voronoi <- deldir(subset_data_csr_ppp)
      
      delsgs_a <- voronoi$delsgs %>%
        dplyr::select(x_center = x1, y_center = y1, x_neigh = x2, y_neigh = y2, id_center = ind1)
      delsgs_b <- voronoi$delsgs %>%
        dplyr::select(x_center = x2, y_center = y2, x_neigh = x1, y_neigh = y1, id_center = ind2)
      
      delsgs <- bind_rows(delsgs_a, delsgs_b)
      
      delsgs$label = paste(delsgs$x_neigh, delsgs$y_neigh, sep="-")
      #delsgs$labelcenter = paste(delsgs$x_center, delsgs$y_center, sep="-")
      #length(unique(delsgs$labelcenter)) # = nrow(subset_data) ensure all individuals as focal individuals
      
      neigh <- delsgs %>% 
        merge(subset_data_csr %>% dplyr::select(label, species_English_neigh = species_English, cover_cm2_neigh = cover_cm2, high_mm_neigh = high_mm , abun_neigh = abun)) %>%
        dplyr::select(X_mm = x_center, Y_mm = y_center, id_center, x_neigh, y_neigh, species_English_neigh, cover_cm2_neigh, high_mm_neigh, abun_neigh) %>%
        merge(subset_data_csr)
      #length(unique(neigh$id_center))   
      
      out_plant_id <- voronoi$dirsgs %>%
        mutate(outn = thirdv1*thirdv2) %>%
        filter(outn < 0) %>%
        dplyr::select(ind1, ind2) %>%
        pivot_longer(1:2, names_to = 'ind', values_to = 'id_center') %>%
        filter(!duplicated(id_center))
      
      neigh_sel <- neigh %>%
        mutate(outy = if_else(id_center %in% out_plant_id$id_center, 0, 1)) %>%
        filter(outy == 1) %>%
        mutate(spe_identity = if_else(species_English_neigh == species_English, '0', '1')) %>%
        mutate(spe_identity = as.numeric(spe_identity))
      
      sub_Msi <- neigh_sel %>%
        group_by(X_mm, Y_mm, species_English) %>%
        dplyr::summarise(Msi = if_else(sum(spe_identity) == length(species_English_neigh), 
                                       sum(spe_identity)*(length(unique(species_English_neigh))+1)/length(species_English_neigh)/(length(species_English_neigh)+1),
                                       sum(spe_identity)*length(unique(species_English_neigh))/length(species_English_neigh)/(length(species_English_neigh)+1))) %>%
        mutate(site = i, plot = j, time = k)
      
      com_Msi = mean(sub_Msi$Msi)
      
      sub_Ms_csr_time <- sub_Msi %>%
        mutate(com_Msi = com_Msi) %>%
        group_by(species_English) %>%
        dplyr::summarise(spe_Msi = mean(Msi),
                         com_Msi = mean(com_Msi)) %>%
        mutate(site = i, plot = j, time = k)
      
      
      sub_Ms_csr <- bind_rows(sub_Ms_csr, sub_Ms_csr_time)
      
      #sub_Msi$xt <- cut(sub_Msi$X_mm,seq(min(sub_Msi$X_mm),round(max(sub_Msi$X_mm),min(sub_Msi$X_mm)),w)) 
      #sub_Msi$yt <- cut(sub_Msi$Y_mm,seq(min(sub_Msi$Y_mm),round(max(sub_Msi$Y_mm),min(sub_Msi$Y_mm)),w))
      
      #sub2_com_mass_mean_msi_csr <- sub_Msi %>%
      #    filter(xt != 'NA') %>%
      #    filter(yt != 'NA') %>%
      #    group_by(site, plot, time, xt, yt, species_English) %>%
      #    dplyr::summarise(Msi=mean(Msi)) %>%
      #    ungroup() %>%
      #    group_by(site, plot, time, xt, yt) %>%
      #    dplyr::summarise(mean_Msi=mean(Msi)) %>%
      #    mutate(scale = w)
      
      #sub_com_mass_mean_msi_csr_scaling <- bind_rows(sub_com_mass_mean_msi_csr_scaling, sub2_com_mass_mean_msi_csr)
    }
    #}
    
    #sub_com_mass_mean_msi_csr_scaling_exp <- sub_com_mass_mean_msi_csr_scaling %>%
    #  group_by(site, plot, scale, xt, yt) %>%
    #  dplyr::summarise(mean_com_Msi=mean(mean_Msi),
    #                   sd_com_Msi=sd(mean_Msi),
    #                   conf_low_com_Msi = t.test(mean_Msi)$conf.int[1],
    #                   conf_high_com_Msi = t.test(mean_Msi)$conf.int[2])
    
    #Ms_all_CSR_exp_com_scaling <- bind_rows(Ms_all_CSR_exp_com_scaling, sub_com_mass_mean_msi_csr_scaling)
    
    sub_Ms_spe_csr_exp <- sub_Ms_csr %>%
      group_by(site, plot, species_English) %>%
      dplyr::summarise(mean_spe_Msi = mean(spe_Msi),
                       sd_spe_Msi = sd(spe_Msi),
                       conf_low_spe_Msi = t.test(spe_Msi)$conf.int[1],
                       conf_high_spe_Msi = t.test(spe_Msi)$conf.int[2])
    
    sub_Ms_com_csr_exp <- sub_Ms_csr %>%
      group_by(site, plot, time) %>%
      dplyr::summarise(com_Msi = mean(com_Msi)) 
    
    sub_Ms_com_csr_exp <- sub_Ms_com_csr_exp %>%
      group_by(site, plot) %>%   
      dplyr::summarise(mean_com_Msi = mean(com_Msi),
                       sd_com_Msi = sd(com_Msi),                 
                       conf_low_com_Msi = t.test(com_Msi)$conf.int[1],
                       conf_high_com_Msi = t.test(com_Msi)$conf.int[2])
    
    sub_Ms_csr_exp <- sub_Ms_spe_csr_exp %>%
      merge(sub_Ms_com_csr_exp)
    
    Ms_all_CSR_exp_site <- bind_rows(Ms_all_CSR_exp_site, sub_Ms_csr_exp)
    
    write_xlsx(Ms_all_CSR_exp_site, 'Ms_all_CSR_exp_site_2024.xlsx')
    #write_xlsx(Ms_all_CSR_exp_com_scaling, 'Ms_all_CSR_exp_com_scaling.xlsx')
  }
}

## expected species richness range based on CSR ----
setwd("D:/OneDrive/Rdata/7 spatial patterns of grassland plant species")
dot.comp_all <- read_excel("all_points_across_six_sites_2024.xlsx",sheet=1,na="NA")
head(dot.comp_all)
headTail(dot.comp_all)
str(dot.comp_all)
summary(dot.comp_all)

dot.comp_all$X_mm <- as.numeric(dot.comp_all$X_mm)
dot.comp_all$Y_mm <- as.numeric(dot.comp_all$Y_mm)

Com_mass_mean_msi_RS_exp_all_sites <- data.frame()
scale <- c(100,150,200,250,300,350,400,450,500,600,700,800,900,1000)

for (i in unique(dot.comp_all$site)) {
  
  Com_mass_mean_msi_all_time_site <- data.frame()
  
  for(j in unique(dot.comp_all$plot)) {
    subset_data <- dot.comp_all[dot.comp_all$site==i&dot.comp_all$plot==j,]
    
    subset_data <- subset_data %>%
      dplyr::select(!year) %>%
      filter(!duplicated(cbind(X_mm, Y_mm)))
    
    subset_data <- subset_data %>%
      mutate(id = 1:nrow(subset_data))
    
    subset_data$species_English <- as.factor(subset_data$species_English) # make species code a factor
    
    Com_mass_mean_msi_all_time <- data.frame()
    
    for (m in 1:199) {
      #for (w in scale) {
      
      n = nrow(subset_data)
      
      random_spe <- subset_data %>%
        dplyr::select(species_English, cover_cm2, high_mm, abun) %>%
        mutate(id = sample(1:n, n))
      
      subset_data2 <- subset_data %>%
        dplyr::select(-c(species_English, cover_cm2, high_mm, abun))
      
      subset_data_csr <- subset_data2 %>%
        merge(random_spe) %>%
        dplyr::select(-id) %>%
        mutate(label = paste(X_mm, Y_mm, sep="-"))
      
      subset_data_csr_ppp <- ppp(subset_data_csr$X_mm,
                                 subset_data_csr$Y_mm,
                                 xrange=c(0,2000),yrange=c(0,2000),
                                 marks=subset_data_csr$species_English)
      
      voronoi <- deldir(subset_data_csr_ppp)
      
      delsgs_a <- voronoi$delsgs %>%
        dplyr::select(x_center = x1, y_center = y1, x_neigh = x2, y_neigh = y2, id_center = ind1)
      delsgs_b <- voronoi$delsgs %>%
        dplyr::select(x_center = x2, y_center = y2, x_neigh = x1, y_neigh = y1, id_center = ind2)
      
      delsgs <- bind_rows(delsgs_a, delsgs_b)
      
      delsgs$label = paste(delsgs$x_neigh, delsgs$y_neigh, sep="-")
      #delsgs$labelcenter = paste(delsgs$x_center, delsgs$y_center, sep="-")
      #length(unique(delsgs$labelcenter)) # = nrow(subset_data) ensure all individuals as focal individuals
      
      neigh <- delsgs %>% 
        merge(subset_data_csr %>% dplyr::select(label, species_English_neigh = species_English, cover_cm2_neigh = cover_cm2, high_mm_neigh = high_mm , abun_neigh = abun)) %>%
        dplyr::select(X_mm = x_center, Y_mm = y_center, id_center, x_neigh, y_neigh, species_English_neigh, cover_cm2_neigh, high_mm_neigh, abun_neigh) %>%
        merge(subset_data_csr)
      #length(unique(neigh$id_center))   
      
      out_plant_id <- voronoi$dirsgs %>%
        mutate(outn = thirdv1*thirdv2) %>%
        filter(outn < 0) %>%
        dplyr::select(ind1, ind2) %>%
        pivot_longer(1:2, names_to = 'ind', values_to = 'id_center') %>%
        filter(!duplicated(id_center))
      
      neigh_sel <- neigh %>%
        mutate(outy = if_else(id_center %in% out_plant_id$id_center, 0, 1)) %>%
        filter(outy == 1) %>%
        mutate(spe_identity = if_else(species_English_neigh == species_English, '0', '1')) %>%
        mutate(spe_identity = as.numeric(spe_identity))
      
      sub_Msi <- neigh_sel %>%
        group_by(X_mm, Y_mm, species_English) %>%
        dplyr::summarise(Msi = if_else(sum(spe_identity) == length(species_English_neigh), 
                                       sum(spe_identity)*(length(unique(species_English_neigh))+1)/length(species_English_neigh)/(length(species_English_neigh)+1),
                                       sum(spe_identity)*length(unique(species_English_neigh))/length(species_English_neigh)/(length(species_English_neigh)+1)),
                         neigh_abun = length(species_English_neigh),
                         neigh_rich = length(unique(species_English_neigh))) 
      
      sub_Ms_all_site <- sub_Msi %>% merge(neigh_sel) #%>% mutate(site = i, plot = j, com_Msi = com_Msi) 
      
      #Ms_all_site <- bind_rows(Ms_all_site, sub_Ms_all_site)
      
      Com_mass_mean_msi <- data.frame()
      
      for (k in scale) {
        sub_Ms_all_site$xt <- cut(sub_Ms_all_site$X_mm,seq(min(sub_Ms_all_site$X_mm),round(max(sub_Ms_all_site$X_mm),min(sub_Ms_all_site$X_mm)),k)) 
        sub_Ms_all_site$yt <- cut(sub_Ms_all_site$Y_mm,seq(min(sub_Ms_all_site$Y_mm),round(max(sub_Ms_all_site$Y_mm),min(sub_Ms_all_site$Y_mm)),k))
        
        sub_Com_mass_mean_msi <- sub_Ms_all_site %>%
          filter(xt != 'NA') %>%
          filter(yt != 'NA') %>%
          dplyr::group_by(xt, yt, X_mm, Y_mm, species_English) %>%
          dplyr::summarise(Msi=mean(Msi),cover_cm2=mean(cover_cm2),high_mm=mean(high_mm)) %>%
          dplyr::ungroup() %>%
          dplyr::group_by(xt, yt) %>%
          dplyr::summarise(rich = length(unique(species_English)),
                           com_mass=log(sum(cover_cm2*high_mm/10)+1),
                           com_Msi=mean(Msi)) %>%
          ungroup() %>%
          dplyr::mutate(scale = k)
        
        Com_mass_mean_msi <- bind_rows(Com_mass_mean_msi, sub_Com_mass_mean_msi)
      }
      
      Com_mass_mean_msi <- Com_mass_mean_msi %>% mutate(time = m)
      Com_mass_mean_msi_all_time <- bind_rows(Com_mass_mean_msi_all_time,Com_mass_mean_msi)
      
    }
    
    Com_mass_mean_msi_all_time <- Com_mass_mean_msi_all_time %>% mutate(plot = j)
    Com_mass_mean_msi_all_time_site <- bind_rows(Com_mass_mean_msi_all_time_site, Com_mass_mean_msi_all_time)
    
  }
  
  Com_mass_mean_msi_exp <- Com_mass_mean_msi_all_time_site %>%
    dplyr::group_by(time,scale) %>%
    dplyr::summarise(mean_rich = mean(rich),
                     mean_com_Msi = mean(com_Msi),
                     mean_com_mass = mean(com_mass),
                     max_rich = max(rich),
                     max_com_Msi = max(com_Msi),
                     max_com_mass = max(com_mass),
                     min_rich = min(rich),
                     min_com_Msi = min(com_Msi),
                     min_com_mass = min(com_mass),
                     med_rich = median(rich),
                     med_com_Msi = median(com_Msi),
                     med_com_mass = median(com_mass)) %>%
    ungroup() %>%
    mutate(D_rich = max_rich-min_rich,
           D_com_Msi = max_com_Msi-min_com_Msi,
           D_com_mass = max_com_mass-min_com_mass) %>%
    dplyr::group_by(scale) %>% 
    dplyr::summarise(mean_D_rich = mean(D_rich),
                     sd_D_rich = sd(D_rich),
                     mean_D_com_Msi = mean(D_com_Msi),
                     sd_D_com_Msi = sd(D_com_Msi),
                     mean_D_com_mass = mean(D_com_mass),
                     sd_D_com_mass = sd(D_com_mass),
                     
                     mean_rich0 = mean(mean_rich),
                     sd_mean_rich = sd(mean_rich),
                     #conf_low_mean_rich = t.test(mean_rich)$conf.int[1],
                     #conf_high_mean_rich = t.test(mean_rich)$conf.int[2],
                     mean_com_Msi0 = mean(mean_com_Msi),
                     sd_mean_com_Msi = sd(mean_com_Msi),
                     #conf_low_mean_com_Msi = t.test(mean_com_Msi)$conf.int[1],
                     #conf_high_mean_com_Msi = t.test(mean_com_Msi)$conf.int[2],
                     mean_com_mass0 = mean(mean_com_mass),
                     sd_mean_com_mass = sd(mean_com_mass),
                     #conf_low_mean_com_mass = t.test(mean_com_mass)$conf.int[1],
                     #conf_high_mean_com_mass = t.test(mean_com_mass)$conf.int[2],
                     
                     mean_max_rich0 = mean(max_rich),
                     sd_max_rich = sd(max_rich),
                     #conf_low_max_rich = t.test(max_rich)$conf.int[1],
                     #conf_high_max_rich = t.test(max_rich)$conf.int[2],
                     mean_max_com_Msi0 = mean(max_com_Msi),
                     sd_max_com_Msi = sd(max_com_Msi),
                     #conf_low_max_com_Msi = t.test(max_com_Msi)$conf.int[1],
                     #conf_high_max_com_Msi = t.test(max_com_Msi)$conf.int[2],
                     mean_max_com_mass0 = mean(max_com_mass),
                     sd_max_com_mass = sd(max_com_mass),
                     #conf_low_max_com_mass = t.test(max_com_mass)$conf.int[1],
                     #conf_high_max_com_mass = t.test(max_com_mass)$conf.int[2],
                     
                     mean_min_rich0 = mean(min_rich),
                     sd_min_rich = sd(min_rich),
                     #conf_low_min_rich = t.test(min_rich)$conf.int[1],
                     #conf_high_min_rich = t.test(min_rich)$conf.int[2],
                     mean_min_com_Msi0 = mean(min_com_Msi),
                     sd_min_com_Msi = sd(min_com_Msi),
                     #conf_low_min_com_Msi = t.test(min_com_Msi)$conf.int[1],
                     #conf_high_min_com_Msi = t.test(min_com_Msi)$conf.int[2],
                     mean_min_com_mass0 = mean(min_com_mass),
                     sd_min_com_mass = sd(min_com_mass),
                     #conf_low_min_com_mass = t.test(min_com_mass)$conf.int[1],
                     #conf_high_min_com_mass = t.test(min_com_mass)$conf.int[2],
                     
                     mean_med_rich0 = mean(med_rich),
                     sd_med_rich = sd(med_rich),
                     #conf_low_med_rich = t.test(med_rich)$conf.int[1],
                     #conf_high_med_rich = t.test(med_rich)$conf.int[2],
                     mean_med_com_Msi0 = mean(med_com_Msi),
                     sd_med_com_Msi = sd(med_com_Msi),
                     #conf_low_med_com_Msi = t.test(med_com_Msi)$conf.int[1],
                     #conf_high_med_com_Msi = t.test(med_com_Msi)$conf.int[2],
                     mean_med_com_mass0 = mean(med_com_mass),
                     sd_med_com_mass = sd(med_com_mass),
                     #conf_low_med_com_mass = t.test(med_com_mass)$conf.int[1],
                     #conf_high_med_com_mass = t.test(med_com_mass)$conf.int[2]
    ) %>%
    mutate(site = i)
  
  Com_mass_mean_msi_RS_exp_all_sites <- bind_rows(Com_mass_mean_msi_RS_exp_all_sites, Com_mass_mean_msi_exp)
  
  write_xlsx(Com_mass_mean_msi_RS_exp_all_sites, 'Rich_Msi_Mass_RS_exp_all_sites_2024.xlsx')
  
}
getwd()
