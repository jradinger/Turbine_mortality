##############################################
###  Mixed model to investigate turbine-related mortality
### Author: Johannes Radinger
##############################################
library(lme4)
library(car)
library(boot)
library(broom)
library(ggplot2)
library(blmeco)
library(lmtest)
library(tinter)
library(see)
library(xtable)

###############################################
# Load data
turb_mort_df <- read.csv("./Data/turb_mort_df.csv",stringsAsFactors=T)
summary(turb_mort_df)
source("./Scripts/00_Helper_functions.R")

#########################################################
### Prepare data for GLMMs
#########################################################
# Observation effect dummy variable (https://rdrr.io/github/markushuff/PsychHelperFunctions/man/overdisp_fun.html)
# to deal with overdispersion
turb_mort_df$obs_effect <- 1:nrow(turb_mort_df)
summary(turb_mort_df)
#Define reference levels
turb_mort_df <- within(turb_mort_df, species_code <- relevel(species_code, ref = "Sal_tru"))
turb_mort_df <- within(turb_mort_df, turbine_type <- relevel(turbine_type, ref = "Kaplan"))
turb_mort_df <- within(turb_mort_df, sampling_method <- relevel(sampling_method, ref = "net"))

# Standarize (mean/sd) single predictors
summary(turb_mort_df$avg_length)
length(turb_mort_df$avg_length[!is.na(turb_mort_df$avg_length)])
scaling_length = mean(turb_mort_df$avg_length,na.rm=T)
turb_mort_df$avg_length_scaled <- as.vector(scale(turb_mort_df$avg_length,scale=scaling_length,center=F))
turb_mort_df$min_length_scaled <- as.vector(scale(turb_mort_df$min_length,scale=scaling_length,center=F))
turb_mort_df$max_length_scaled <- as.vector(scale(turb_mort_df$max_length,scale=scaling_length,center=F))


# Get triangle-probability weights for resampling
turb_mort_df$triangle_weight_length <- 1
turb_mort_df$triangle_weight_mort <- 1

# Number of cases with reported length information
nrow(turb_mort_df[!is.na(turb_mort_df$avg_length),])
summary(turb_mort_df$avg_length)

# Make complete dataframe
mort_analysis_df <- completeFun(turb_mort_df, c("n_mort", "n_not_mort","mort_rate",
                                                 "post_monitoring_uncertainty","catch_related_uncertainty",
                                                "avg_length","capacity_class",
                                                "triangle_weight_length","triangle_weight_mort",
                                                "avg_length_scaled","min_length_scaled","max_length_scaled",
                                                "sampling_method",
                                                "turbine_type","family","location",
                                                "species_code","experiment","obs_effect"))
mort_analysis_df <- droplevels(mort_analysis_df[!mort_analysis_df$turbine_type %in% c("Other"),])
table(mort_analysis_df$turbine_type)
summary(mort_analysis_df)

########################
### Original model(s)
########################
glmm_A <- glmer(cbind(n_mort,n_not_mort)~avg_length_scaled+turbine_type+
                  (1|location)+(1|location:experiment)+
                  (1|family)+(1|family:species_code)+
                  (1|sampling_method)+
                  (1|obs_effect),
                family=binomial,
                data=mort_analysis_df,
                control=glmerControl(calc.derivs=FALSE))
glmm_B <- glmer(cbind(n_mort,n_not_mort)~avg_length_scaled*turbine_type+
                  (1|location)+(1|location:experiment)+
                  (1|family)+(1|family:species_code)+
                  (1|sampling_method)+
                  (1|obs_effect),
                family=binomial,
                data=mort_analysis_df,
                control=glmerControl(calc.derivs=FALSE))
summary(glmm_B)
vif(glmm_B) 
glmm_C <- glmer(cbind(n_mort,n_not_mort)~avg_length_scaled+
                  (1|location)+(1|location:experiment)+
                  (1|family)+(1|family:species_code)+
                  (1|sampling_method)+
                  (1|obs_effect),
                family=binomial,
                data=mort_analysis_df,
                control=glmerControl(calc.derivs=FALSE))
glmm_D <- glmer(cbind(n_mort,n_not_mort)~turbine_type+
                  (1|location)+(1|location:experiment)+
                  (1|family)+(1|family:species_code)+
                  (1|sampling_method)+
                  (1|obs_effect),
                family=binomial,
                data=mort_analysis_df,
                control=glmerControl(calc.derivs=FALSE))
glmm_E <- glmer(cbind(n_mort,n_not_mort)~avg_length_scaled*turbine_type+capacity_class+
                   (1|location)+(1|location:experiment)+
                   (1|family)+(1|family:species_code)+
                  (1|sampling_method)+
                   (1|obs_effect),
                 family=binomial,
                 data=mort_analysis_df,
                 control=glmerControl(calc.derivs=FALSE))

# Likelihoodratio-Test whether single variables and or interaction terms are relevant
lrt_results <- anova(glmm_A,glmm_B,glmm_C,glmm_D,glmm_E)
anova(glmm_B,glmm_E)
anova(glmm_E,glmm_F)
print(xtable(lrt_results, type = "latex"), file = "./Results_Figs_Tabs/LRT_anova_table.tex")
### --> Complete model glmm_F with length and turbine type and their interaction and sampling method is 
#### best model based on maximum-likelihood ratio test

### To assess overdispersion
dispersion_glmer(glmm_F) #should not be over 1.4 https://www.rdocumentation.org/packages/blmeco/versions/1.4/topics/dispersion_glmer 
overdisp_fun(glmm_F) #https://rdrr.io/github/markushuff/PsychHelperFunctions/man/overdisp_fun.html 

################################
#### Preparation for prediction 
################################
# Dataframe for predict
predict_df_families <- cbind(n_mort = 0,
                             n_not_mort = 0,
                             avg_length = 20,
                             avg_length_scaled = scale(20,scale=scaling_length,center=F),
                             expand.grid(family = levels(mort_analysis_df$family),
                                         turbine_type = levels(mort_analysis_df$turbine_type),
                                         sampling_method = levels(mort_analysis_df$sampling_method)),
                             species_code=NA,obs_effect=NA,location=NA,experiment=NA)
predict_df_families$predict_id <- paste("family_pred",seq(1,nrow(predict_df_families)),sep="_")

predict_df_length <- cbind(n_mort = 0,
                           n_not_mort = 0,
                           expand.grid(avg_length = seq(1,105,1),
                                       turbine_type = levels(mort_analysis_df$turbine_type)),
                           family = NA,species_code=NA,obs_effect=NA,location=NA,experiment=NA)
predict_df_length$avg_length_scaled <- scale(predict_df_length$avg_length,scale=scaling_length,center=F)
predict_df_length$predict_id <- paste("length_pred",seq(1,nrow(predict_df_length)),sep="_")


#### Function to collect results of model
predict_boot <- function(model){
  res <- c(fixef(model),
    predict(model,newdata=predict_df_families,
              re.form=~(1|family),type="response",
              allow.new.levels=TRUE),
    predict(model,newdata=predict_df_length,
              re.form=~0,type="response",
              allow.new.levels=TRUE))
  names(res) <- c(names(fixef(model)),
                  predict_df_families$predict_id,
                  predict_df_length$predict_id)
  return(res)
}

######################################
#### Bootstrapping model
#### Run repeated model runs on bootstrapped cases_within_cases
######################################
# Number of bootstrap replicates
B <- 10000
# Generate bootstrapped replicate dataframes
rep.data <- list()
i <- 1
#for(i in (1:B)){
while(i <= B){
  print(i)
  rep.data_i <- try(.cases.within.cases.resamp(dat=mort_analysis_df,
                                               #cluster = c(rev(names(lme4::getME(glmm_B, "flist"))), ".id"),
                                               #resample = c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE),
                                               cluster = c("experiment",".id"),
                                               resample = c(TRUE,TRUE),
                                               resample.mort = list(n_mort = "n_mort",
                                                                    n_not_mort = "n_not_mort",
                                                                    uncertainty_lwr = "catch_related_uncertainty",
                                                                    uncertainty_upr = "post_monitoring_uncertainty",
                                                                    triangle_weight_mort = "triangle_weight_mort"),
                                               resample.length = list(avg_length = "avg_length_scaled",
                                                                      min_length="min_length_scaled",
                                                                      max_length="max_length_scaled",
                                                                      triangle_weight_length = "triangle_weight_length")),
                    silent=TRUE)
  if(isTRUE(class(rep.data_i)=="try-error")) {
    next}else{
      # Update mort_rate from resampled data
      rep.data_i$mort_rate <- rep.data_i$n_mort/(rep.data_i$n_mort+rep.data_i$n_not_mort)
      
      # Collect results
      rep.data[[i]] <- rep.data_i
      i <- i+1}
}
saveRDS(rep.data,paste("./Results_Figs_Tabs/rep.data_","10000_IGB",".RDS",sep=""))

###############################
# Run replicated GLMM runs on bootstrapped dataframes
###############################
bootstrap_results_glmm <- .cases.completion(model=glmm_B, rep.data=rep.data, fn=predict_boot)
saveRDS(bootstrap_results_glmm,paste("./Results_Figs_Tabs/bootstrap_results_","10000_IGB",".RDS",sep=""))
#bootstrap_results_glmm <- readRDS("./Results_Figs_Tabs/bootstrap_results_1000_IGB.RDS")

###########################################
## Analyse bootstrapped model results (e.g. CI)
##########################################
## Get bootstrapped confidence intervals for model statistics and predicted values
bc_df <- as.data.frame(tidy(bootstrap_results_glmm,conf.int=TRUE,conf.method="perc"))
bootstrap_results_glmm$t0
bc_df_coef <- bc_df[!grepl(".*pred.*",bc_df$term),]
bc_df_length <- bc_df[grepl("^length_pred.*",bc_df$term),]
bc_df_length <- merge(bc_df_length,predict_df_length,by.x="term",by.y="predict_id")

# Moratility for 25 cm fish across turbine types
bc_df_length[bc_df_length$avg_length == 25,]

##########################################
### Plot predicted relationship length vs. moratility across turbine types (incl 95% CI)
##########################################
# Level ordering
mort_analysis_df$turbine_type <- factor(mort_analysis_df$turbine_type,
                                        levels=c("Kaplan",
                                                 "Francis",
                                                 "VLH",
                                                 "Screw",
                                                 "Cross-flow",
                                                 "Water wheel"))
bc_df_length$turbine_type <- factor(bc_df_length$turbine_type,
                                    levels=c("Kaplan",
                                             "Francis",
                                             "VLH",
                                             "Screw",
                                             "Cross-flow",
                                             "Water wheel"))

purple_shades <- sample(tinter(material_colors("deep purple"), steps = 6, direction="both", crop = 3),6)
# Facet label order must much level ordering
facet_label_names <- list(
  "Kaplan" = "Kaplan" ,
  "Francis" = "Francis",
  "VLH" = "VLH",
  "Screw" = "Archimedes' screw",
  "Cross-flow" = "Cross-flow",
  "Water wheel" = "Water wheel")
facet_labeller <- function(variable,value){
  return(facet_label_names[value])
}

tt="Cross-flow"
plot(mort_analysis_df$avg_length[mort_analysis_df$turbine_type==tt],
     mort_analysis_df$mort_rate[mort_analysis_df$turbine_type==tt])

cairo_pdf("./Results_Figs_Tabs/Length_mort_relationship.pdf",width=6,height=6)
ggplot(data=bc_df_length,
       aes(x=avg_length,
           y=statistic))+
  geom_point(data=mort_analysis_df,aes(x=avg_length,y=mort_rate,group=turbine_type),
             shape=21,alpha=40,colour="grey",size=1.1)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),
             fill=material_colors("deep purple"),alpha=0.3)+
  scale_fill_manual(values=purple_shades)+
  geom_line(color=material_colors("deep purple"))+
  scale_color_manual(values=purple_shades)+
  facet_wrap(~turbine_type,labeller=facet_labeller,ncol=2)+
  scale_x_continuous(name="Fish length [m]",breaks=c(0,25,50,75,100),labels=c(0,25,50,75,100)/100, limits=c(0,105))+
  scale_y_continuous(name="Mortality rate",limits=c(0,1))+
  theme_bw()+
  theme(legend.title = element_blank(),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

