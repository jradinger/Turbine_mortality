##############################################
### Mortality of fish passing turbines
### Calculate bootstrapped confidence interval of turbine-related mortality based on individual counts
### Author: Johannes Radinger#
### October 2021
##############################################

# load libaries
library(boot)
library(broom)
library(xtable)
library(dplyr)

###############################################
# Load data
turb_mort_df <- read.csv("./Data/turb_mort_df.csv",stringsAsFactors=T)
source("./Scripts/00_Helper_functions.R")

###############################################
### Descriptive stats on uncertainty classes
###############################################
turb_mort_df$combined_uncertainty <- as.factor(paste(turb_mort_df$post_monitoring_uncertainty,
                                           turb_mort_df$catch_related_uncertainty,sep="_"))
prop.table(table(turb_mort_df$post_monitoring_uncertainty))
prop.table(table(turb_mort_df$catch_related_uncertainty))
nrow(turb_mort_df)

###############################################
##### Bootstrap average mortality and extract CI
###############################################
# Add triangle weigth for resampling
turb_mort_df$triangle_weight_mort <- 1

# Generate bootstrapped replicate dataframes
B <- 10000
rep.data <- list()
i <- 1
while(i <= B){
  print(i)
  rep.data_i <- try(.cases.within.cases.resamp(dat=turb_mort_df,
                                               cluster = c("experiment",".id"),
                                               resample = c(TRUE,TRUE),
                                               resample.mort = list(n_mort = "n_mort",
                                                                    n_not_mort = "n_not_mort",
                                                                    uncertainty_lwr = "catch_related_uncertainty",
                                                                    uncertainty_upr = "post_monitoring_uncertainty",
                                                                    triangle_weight_mort = "triangle_weight_mort")),
                    silent=TRUE)
  if(isTRUE(class(rep.data_i)=="try-error")) {
    next}else{
      # Update mort_rate from resampled data
      rep.data_i$mort_rate <- rep.data_i$n_mort/(rep.data_i$n_mort+rep.data_i$n_not_mort)
        
      # Collect results
      rep.data[[i]] <- rep.data_i
      i <- i+1}
}

#saveRDS(rep.data,"./Results_Figs_Tabs/bootstrapped_turb_mort_df.rds")
#rep.data <- loadRDS("./Results_Figs_Tabs/bootstrapped_turb_mort_df.rds")

###############################
## Bootstrapped replicate dataframes subset of handling corrected data with postmonitoring (i.e. gold standard)
###############################
turb_mort_df_subset <- 
  turb_mort_df%>%
  filter(catch_related_uncertainty==0 &
           post_monitoring_uncertainty==0)

rep.data_subset <- list()
for(i in 1:length(rep.data)){
  #i=1
  print(paste("Subset replicate: ",as.character(i),sep=""))
  rep.data_subset[[i]] <- 
    rep.data[[i]]%>%
    filter(catch_related_uncertainty==0 &
             post_monitoring_uncertainty==0)
}

###############################
# Mean and CI of bootstrapped dataframes
###############################
# Overall mean
bootstrap_avg_mort <- bootstrap_avg_mort_weighted(orig.data=turb_mort_df,
                                                                 rep.data=rep.data,group_var=NA)
as.data.frame(tidy(bootstrap_avg_mort,conf.int=TRUE,conf.method="perc"))
sd(turb_mort_df$mort_rate)

# Across classes of uncertainty
bootstrap_avg_mort_combined_uncertainty_classes <- bootstrap_avg_mort_weighted(orig.data=turb_mort_df,
                                                                               rep.data=rep.data,group_var="combined_uncertainty")
bootstrap_avg_mort_combined_uncertainty_classes_df <- multi.boot.ci(bootstrap_avg_mort_combined_uncertainty_classes)
bootstrap_avg_mort_combined_uncertainty_classes_df <- cbind(bootstrap_avg_mort_combined_uncertainty_classes_df,n=as.vector(table(turb_mort_df$combined_uncertainty)))
bootstrap_avg_mort_combined_uncertainty_classes_df$CI <- paste("(",round(bootstrap_avg_mort_combined_uncertainty_classes_df$conf.low,2),
                                                               ", ",round(bootstrap_avg_mort_combined_uncertainty_classes_df$conf.high,2),
                                                               ")",sep="")
bootstrap_avg_mort_combined_uncertainty_classes_df$post_monitoring_uncertainty <- unlist(purrr::map(strsplit(rownames(bootstrap_avg_mort_combined_uncertainty_classes_df),"_"),1))
bootstrap_avg_mort_combined_uncertainty_classes_df$catch_related_uncertainty <- unlist(purrr::map(strsplit(rownames(bootstrap_avg_mort_combined_uncertainty_classes_df),"_"),2))

print(xtable(bootstrap_avg_mort_combined_uncertainty_classes_df[,c("post_monitoring_uncertainty",
                                                                   "catch_related_uncertainty",
                                                                   "statistic",
                                                                   "CI",
                                                                   "n")], type = "latex"), file = "./Results_Figs_Tabs/bootstrap_avg_mort_combined_uncertainty_classes.tex")


##########################
# Capacity class
bootstrap_avg_mort_capacity_class <- bootstrap_avg_mort_weighted(orig.data=turb_mort_df,
                                         rep.data=rep.data,group_var="capacity_class")
bootstrap_avg_mort_capacity_class_df <- multi.boot.ci(bootstrap_avg_mort_capacity_class)
bootstrap_avg_mort_capacity_class_df <- cbind(bootstrap_avg_mort_capacity_class_df,n=as.vector(table(turb_mort_df$capacity_class)))
bootstrap_avg_mort_capacity_class_df$CI <- paste("(",round(bootstrap_avg_mort_capacity_class_df$conf.low,2),
                                                 ", ",round(bootstrap_avg_mort_capacity_class_df$conf.high,2),
                                                 ")",sep="")
print(xtable(bootstrap_avg_mort_capacity_class_df[,c(1,5,4)], type = "latex"), file = "./Results_Figs_Tabs/bootstrap_avg_mort_capacity_class.tex")

##########################
# Species
bootstrap_avg_mort_species <- bootstrap_avg_mort_weighted(orig.data=turb_mort_df,
                                                                 rep.data=rep.data,group_var="species")
bootstrap_avg_mort_species_df <- multi.boot.ci(bootstrap_avg_mort_species)
bootstrap_avg_mort_species_df <- cbind(bootstrap_avg_mort_species_df,n=as.vector(table(turb_mort_df$species)))
bootstrap_avg_mort_species_df$CI <- paste("(",round(bootstrap_avg_mort_species_df$conf.low,2),
                                          ", ",round(bootstrap_avg_mort_species_df$conf.high,2),
                                          ")",sep="")
bootstrap_avg_mort_species_subset <- bootstrap_avg_mort_weighted(orig.data=turb_mort_df_subset,
                                                                 rep.data=rep.data_subset,group_var="species")
bootstrap_avg_mort_species_subset_df <- multi.boot.ci(bootstrap_avg_mort_species_subset)
bootstrap_avg_mort_species_subset_df <- cbind(bootstrap_avg_mort_species_subset_df,n=as.vector(table(turb_mort_df_subset$species)))
bootstrap_avg_mort_species_subset_df$CI <- paste("(",round(bootstrap_avg_mort_species_subset_df$conf.low,2),
                                                 ", ",round(bootstrap_avg_mort_species_subset_df$conf.high,2),
                                                 ")",sep="")
colnames(bootstrap_avg_mort_species_subset_df) <- paste(colnames(bootstrap_avg_mort_species_subset_df),"_subset",sep="")
bootstrap_avg_mort_species_df_combined <-
  merge(bootstrap_avg_mort_species_df,
      bootstrap_avg_mort_species_subset_df[!is.na(bootstrap_avg_mort_species_subset_df$statistic_subset),],
      by="row.names",
      all.x=T)
print(xtable(bootstrap_avg_mort_species_df_combined[,c(1,2,6,5,7,11,10)], type = "latex"),
      file = "./Results_Figs_Tabs/bootstrap_avg_mort_species.tex",
      include.rownames=F)

##########################
# Family
bootstrap_avg_mort_family <- bootstrap_avg_mort_weighted(orig.data=turb_mort_df,
                                                          rep.data=rep.data,group_var="family")
bootstrap_avg_mort_family_df <- multi.boot.ci(bootstrap_avg_mort_family)
bootstrap_avg_mort_family_df <- cbind(bootstrap_avg_mort_family_df,n=as.vector(table(turb_mort_df$family)))
bootstrap_avg_mort_family_df$CI <- paste("(",round(bootstrap_avg_mort_family_df$conf.low,2),
                                          ", ",round(bootstrap_avg_mort_family_df$conf.high,2),
                                          ")",sep="")
bootstrap_avg_mort_family_subset <- bootstrap_avg_mort_weighted(orig.data=turb_mort_df_subset,
                                                                 rep.data=rep.data_subset,group_var="family")
bootstrap_avg_mort_family_subset_df <- multi.boot.ci(bootstrap_avg_mort_family_subset)
bootstrap_avg_mort_family_subset_df <- cbind(bootstrap_avg_mort_family_subset_df,n=as.vector(table(turb_mort_df_subset$family)))
bootstrap_avg_mort_family_subset_df$CI <- paste("(",round(bootstrap_avg_mort_family_subset_df$conf.low,2),
                                                 ", ",round(bootstrap_avg_mort_family_subset_df$conf.high,2),
                                                 ")",sep="")
colnames(bootstrap_avg_mort_family_subset_df) <- paste(colnames(bootstrap_avg_mort_family_subset_df),"_subset",sep="")
bootstrap_avg_mort_family_df_combined <-
  merge(bootstrap_avg_mort_family_df,
        bootstrap_avg_mort_family_subset_df[!is.na(bootstrap_avg_mort_family_subset_df$statistic_subset),],
        by="row.names",
        all.x=T)
print(xtable(bootstrap_avg_mort_family_df_combined[,c(1,2,6,5,7,11,10)], type = "latex"),
      file = "./Results_Figs_Tabs/bootstrap_avg_mort_family.tex",
      include.rownames=F)

##########################
# Turbine type
bootstrap_avg_mort_turbine_type <- bootstrap_avg_mort_weighted(orig.data=turb_mort_df,
                                                         rep.data=rep.data,group_var="turbine_type")
bootstrap_avg_mort_turbine_type_df <- multi.boot.ci(bootstrap_avg_mort_turbine_type)
bootstrap_avg_mort_turbine_type_df <- cbind(bootstrap_avg_mort_turbine_type_df,n=as.vector(table(turb_mort_df$turbine_type)))
bootstrap_avg_mort_turbine_type_df$CI <- paste("(",round(bootstrap_avg_mort_turbine_type_df$conf.low,2),
                                         ", ",round(bootstrap_avg_mort_turbine_type_df$conf.high,2),
                                         ")",sep="")
bootstrap_avg_mort_turbine_type_subset <- bootstrap_avg_mort_weighted(orig.data=turb_mort_df_subset,
                                                                rep.data=rep.data_subset,group_var="turbine_type")
bootstrap_avg_mort_turbine_type_subset_df <- multi.boot.ci(bootstrap_avg_mort_turbine_type_subset)
bootstrap_avg_mort_turbine_type_subset_df <- cbind(bootstrap_avg_mort_turbine_type_subset_df,n=as.vector(table(turb_mort_df_subset$turbine_type)))
bootstrap_avg_mort_turbine_type_subset_df$CI <- paste("(",round(bootstrap_avg_mort_turbine_type_subset_df$conf.low,2),
                                                ", ",round(bootstrap_avg_mort_turbine_type_subset_df$conf.high,2),
                                                ")",sep="")
colnames(bootstrap_avg_mort_turbine_type_subset_df) <- paste(colnames(bootstrap_avg_mort_turbine_type_subset_df),"_subset",sep="")
bootstrap_avg_mort_turbine_type_df_combined <-
  merge(bootstrap_avg_mort_turbine_type_df,
        bootstrap_avg_mort_turbine_type_subset_df[!is.na(bootstrap_avg_mort_turbine_type_subset_df$statistic_subset),],
        by="row.names",
        all.x=T)
print(xtable(bootstrap_avg_mort_turbine_type_df_combined[,c(1,2,6,5,7,11,10)], type = "latex"),
      file = "./Results_Figs_Tabs/bootstrap_avg_mort_turbine_type.tex",
      include.rownames=F)

##########################
# Across classes of uncertainty
bootstrap_avg_mort_combined_uncertainty_classes <- bootstrap_avg_mort_weighted(orig.data=turb_mort_df,
                                                                          rep.data=rep.data,group_var="combined_uncertainty")
bootstrap_avg_mort_combined_uncertainty_classes_df <- multi.boot.ci(bootstrap_avg_mort_combined_uncertainty_classes)
bootstrap_avg_mort_combined_uncertainty_classes_df <- cbind(bootstrap_avg_mort_combined_uncertainty_classes_df,n=as.vector(table(turb_mort_df$combined_uncertainty)))
bootstrap_avg_mort_combined_uncertainty_classes_df$CI <- paste("(",round(bootstrap_avg_mort_combined_uncertainty_classes_df$conf.low,2),
                                                          ", ",round(bootstrap_avg_mort_combined_uncertainty_classes_df$conf.high,2),
                                                          ")",sep="")
print(xtable(bootstrap_avg_mort_combined_uncertainty_classes_df[,c(1,5,4)], type = "latex"), file = "./Results_Figs_Tabs/bootstrap_avg_mort_combined_uncertainty_classes.tex")
