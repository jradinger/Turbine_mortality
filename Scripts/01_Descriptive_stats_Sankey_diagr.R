####################################################
### Mortality of fish passing turbines
### Descriptive statistics of turbine mortality dataset
### Author: Johannes Radinger
### October 2021
####################################################

## Load libraries
library(devtools)
library(ggplot2)
library(randomcoloR)
library(plyr)
library(dplyr)
library(forcats)
#install_github("Displayr/flipPlots")
library(flipPlots)
#install_github("schmidtchristoph/js2graphic")
library(js2graphic)
library(webshot)
library(htmlwidgets)

############################################### 
# Load data
turb_mort_df <- read.csv("./Data/turb_mort_df.csv",stringsAsFactors=T)
turb_df <- read.csv("./Data/turb_df.csv",stringsAsFactors=T)

source("./Scripts/00_Helper_functions.R")

###############################################
## Kaplan, LHP, Salmoniform +3 missing
# Aggregate individuals of less studied families into family "others"
turb_mort_df$family_agg <- turb_mort_df$family
sum_family <- aggregate(turb_mort_df$N,by=list(family = turb_mort_df$family),FUN=sum)
turb_mort_df$family_agg <- fct_collapse(turb_mort_df$family_agg,
                            "Others" = as.character(sum_family$family[sum_family$x<500]))

# Aggregate individuals of less studied orders into order "others"
turb_mort_df$taxorder_agg <- turb_mort_df$taxorder
sum_taxorder <- aggregate(turb_mort_df$N,by=list(taxorder = turb_mort_df$taxorder),FUN=sum)
turb_mort_df$taxorder_agg <- fct_collapse(turb_mort_df$taxorder_agg,
                                        "Others" = as.character(sum_taxorder$taxorder[sum_taxorder$x<5000]))


# Group into length classes
summary(turb_mort_df$avg_length)
turb_mort_df$length_class <- cut(turb_mort_df$avg_length,
                                        breaks=c(0,10,30,50,Inf),
                                        labels=c("<10","10-30","30-50",">50"))
turb_mort_df$length_class <- fct_explicit_na(turb_mort_df$length_class,
                                             na_level="Length not reported")

#########################
# Descriptive stats
#########################
# N assessments
nrow(turb_mort_df)
# N studies
length(unique(turb_mort_df$study))
# N separate assessments
length(unique(turb_mort_df$study_ID))
# N locations
length(unique(turb_mort_df$location))
# N countries
length(unique(turb_mort_df$country))
# N per global region
table(turb_mort_df$global_region) 
# N per sampling method
table(turb_mort_df$sampling_method)
# N per turbine type
table(turb_mort_df$turbine_type)

#########################
##### Descriptive stats across studies (not cases)
#########################
# Distribution capacity class
table(turb_df$capacity_class)
# Distribution capacity class
sum(table(turb_df$turbine_type))
# Distribution global region
table(turb_df$global_region)
# Turbine type
table(turb_df$turbine_type)
# sampling methods
table(turb_df$sampling_method)


# Total N fish
sum(turb_mort_df$N)
  sum(turb_mort_df$n_mort)
(sum(turb_mort_df$n_mort)/sum(turb_mort_df$N))*100
sum(turb_mort_df$n_not_mort)
(sum(turb_mort_df$n_not_mort)/sum(turb_mort_df$N))*100

# N species/families
all_species <- unique(turb_mort_df$species)
all_species <- all_species[order(all_species)]
length(all_species[all_species!="Spp. spp."])-1 #-1 for Leuciscus spp.
length(unique(turb_mort_df$family[!is.na(turb_mort_df$family)]))-1 #-1 for Family = "Spp."
length(unique(turb_mort_df$taxorder[!is.na(turb_mort_df$taxorder)]))-1 #-1 for order = "Spp."

# Overview reported mortality rate across all species
summary(turb_mort_df$mort_rate)
sd(turb_mort_df$mort_rate)
weighted.mean(turb_mort_df$mort_rate,turb_mort_df$N)

# Mort rate per turbine type (no uncertainties)
aggregate(turb_mort_df$mort_rate,by=list(turb_mort_df$turbine_type),FUN=summary)

## Summary by spp (unweighted mean, no uncertainties)
summary_by_spp <- turb_mort_df %>% 
  group_by(species_code) %>% 
  summarize(mean = mean(mort_rate),
            sd = sd(mort_rate),
            n = length(mort_rate))
summary_by_spp <- as.data.frame(summary_by_spp)



#########################################################
## Alluvial diagramm / Sankey diagramm (flipPlots:Displayr)
#########################################################
# Expand df to relect 1 row = 1 fish
df.expanded_mort <- turb_mort_df[rep(row.names(turb_mort_df), turb_mort_df$n_mort), ]
df.expanded_mort$mort <- "mort"
df.expanded_not_mort <- turb_mort_df[rep(row.names(turb_mort_df), turb_mort_df$n_not_mort), ]
df.expanded_not_mort$mort <- "not_mort"
turb_mort_df_expand <- rbind(df.expanded_mort,df.expanded_not_mort)
turb_mort_df_expand$mort <- as.factor(turb_mort_df_expand$mort)
colnames(turb_mort_df_expand)
nrow(turb_mort_df_expand)

# Relative frequency per order
taxorder_prop <- prop.table(table(turb_mort_df_expand$taxorder))*100
taxorder_prop[order(-taxorder_prop)]

# Relative frequency per family
fam_prop <- prop.table(table(turb_mort_df_expand$family))*100
fam_prop[order(-fam_prop)]

# Relative frequency per species
spp_prop <- prop.table(table(turb_mort_df_expand$species_code))*100
spp_prop <- spp_prop[order(-spp_prop)]
round(spp_prop,2)


## Define variables to plot in alluvial diagram
#Sankey_variables <- c("family_agg","capacity_class","turbine_type","mort")
Sankey_variables <- c("taxorder_agg","capacity_class","turbine_type","mort")

apply(turb_mort_df_expand[,Sankey_variables],2,FUN=function(x){levels(x)})

# Create alluvial diagram
S <- SankeyDiagram(data = turb_mort_df_expand[,Sankey_variables], 
              links.and.nodes = NULL,
              output.data.only = F, 
              max.categories = sum(unlist(lapply(sapply(turb_mort_df_expand[,Sankey_variables], levels),length))), 
              subset = NULL,
              weights = NULL, font.size = 12, font.family = "Sans",
              font.unit = "px", 
              colors = c(randomColor(count = length(levels(turb_mort_df_expand[,Sankey_variables[1]])),
                                     hue = "blue",luminosity="bright"),
                         randomColor(count = length(levels(turb_mort_df_expand[,Sankey_variables[2]])),
                                     hue = "red",luminosity="bright"),
                         randomColor(count = length(levels(turb_mort_df_expand[,Sankey_variables[3]])),
                                     hue = "yellow",luminosity="bright"),
                         "red","green"), 
              link.color = c("None", "Source","Target", "First variable", "Last variable")[2],
              variables.share.values = FALSE, label.show.varname = FALSE,
              label.max.length = 100, label.show.counts = T,
              label.show.percentages = FALSE, node.width = 30, node.padding = 7,
              sinks.right = TRUE, hovertext.show.percentages = T)
S
# Save Alluvial diagram as html
#saveWidget(S, paste(getwd(),"/Results_Figs_Tabs/turbine_alluvial_diagram_taxorder.html",sep=""))
### --> Extract svg/pdf from html using the chrome extension SVG Crowbar 2



