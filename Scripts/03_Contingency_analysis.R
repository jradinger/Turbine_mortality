#########################################################
## Contingency tabes, associatiation analysis
#########################################################
library(rcompanion)
library(MASS)
library(xtable)

###############################################
# Load data
turb_mort_df <- read.csv("./Data/turb_mort_df.csv",stringsAsFactors=T)
source("./Scripts/00_ACT.r")

# Set parameters for ACT analysis
alpha <- 0.05
nrep <- 50000

# Do not consider the rare turbine types (otherwise to sparse crosstable for analysis)
table(turb_mort_df$turbine_type)
turb_mort_df_ct <- turb_mort_df[turb_mort_df$turbine_type!="Other",]
turb_mort_df_ct <- droplevels(turb_mort_df_ct)

# reorder levels
turb_mort_df_ct$turbine_type <- factor(turb_mort_df_ct$turbine_type,
                                       levels = c("Kaplan", "Francis", "VLH",
                                                  "Screw","Cross-flow","Water wheel"))
####################################
# Capactiy class vs. turbine type
####################################
cramerV(x=turb_mort_df_ct$capacity_class,
        y=turb_mort_df_ct$turbine_type,
        bias.correct = T)
ct_capacity_class_turbine_type <- table(turb_mort_df_ct$capacity_class,turb_mort_df_ct$turbine_type)

# ACT
report.ACT_capacity_class_turbine_type <- ACT_I(ct_capacity_class_turbine_type, alpha=alpha, Rtype='ADJ', nrep=nrep)
report.ACT_capacity_class_turbine_type
report.ACT_capacity_class_turbine_type$Famwise_Significant

# Prepare output data or residual analysis for MS
ACT_capacity_class_turbine_type_cell_sign <- as.matrix(round(report.ACT_capacity_class_turbine_type$Residuals,2))
ACT_capacity_class_turbine_type_cell_sign_df_ms <- as.data.frame.matrix(ifelse(report.ACT_capacity_class_turbine_type$Famwise_Significant, 
                                                                               paste(ACT_capacity_class_turbine_type_cell_sign, "*",sep = ""),
                                                                               paste(ACT_capacity_class_turbine_type_cell_sign, "",sep = "")))
ACT_capacity_class_turbine_type_expected_freq_df_ms <- as.data.frame.matrix(as.matrix(round(report.ACT_capacity_class_turbine_type$ExpectedFrequencies,2)))

# Output dataframe
capacity_class_turbine_type_df_ms <- rbind(as.data.frame.matrix(ct_capacity_class_turbine_type),
                                           ACT_capacity_class_turbine_type_expected_freq_df_ms,
                                           ACT_capacity_class_turbine_type_cell_sign_df_ms)
print(xtable(capacity_class_turbine_type_df_ms, type = "latex"), file = "./Results_Figs_Tabs/ACT_capacity_class_turbine_type.tex")


####################################
# Turbine type vs. family
####################################
#Only select families with >=5 cases
fam_excl <- levels(turb_mort_df_ct$family)[table(turb_mort_df_ct$family)<5]
turb_mort_df_ct2 <- turb_mort_df_ct[!(turb_mort_df_ct$family %in% fam_excl),]
turb_mort_df_ct2$family <- droplevels(turb_mort_df_ct2$family)

cramerV(x=turb_mort_df_ct2$turbine_type,
        y=turb_mort_df_ct2$family,
        bias.correct = T)
ct2_turbine_type_family <- table(turb_mort_df_ct2$turbine_type,turb_mort_df_ct2$family)

# ACT
report.ACT_turbine_type_family <- ACT_I(ct2_turbine_type_family, alpha=alpha, Rtype='ADJ', nrep=nrep)
report.ACT_turbine_type_family
## Omnibus H of independence rejected --> Alternative H (dependence) accepted; i.e. there is some relatoinship between turbine type and family
report.ACT_turbine_type_family$Famwise_Significant

# Prepare output data or residual analysis for MS
ACT_turbine_type_family_cell_sign <- as.matrix(round(report.ACT_turbine_type_family$Residuals,2))
ACT_turbine_type_family_cell_sign_df_ms <- as.data.frame.matrix(ifelse(report.ACT_turbine_type_family$Famwise_Significant, 
                                                                       paste(ACT_turbine_type_family_cell_sign, "*",sep = ""),
                                                                       paste(ACT_turbine_type_family_cell_sign, "",sep = "")))
ACT_turbine_type_family_expected_freq_df_ms <- as.data.frame.matrix(as.matrix(round(report.ACT_turbine_type_family$ExpectedFrequencies,2)))

# Output dataframe
turbine_type_family_df_ms <- rbind(as.data.frame.matrix(ct2_turbine_type_family),
                                   ACT_turbine_type_family_expected_freq_df_ms,
                                   ACT_turbine_type_family_cell_sign_df_ms)
#turbine_type_family_df_ms <- turbine_type_family_df_ms[,c(colnames(turbine_type_family_df_ms)[colnames(turbine_type_family_df_ms)!='Spp.'],'Spp.')]
print(xtable(turbine_type_family_df_ms, type = "latex"), file = "./Results_Figs_Tabs/ACT_turbine_type_family.tex")


# largest disrepancies within a family
colSums(abs(ACT_turbine_type_family_cell_sign),na.rm=T)
rowSums(abs(ACT_turbine_type_family_cell_sign),na.rm=T)



####################################
# Capacity class vs. family
####################################
cramerV(x=turb_mort_df_ct2$capacity_class,
        y=turb_mort_df_ct2$family,
        bias.correct = T)
ct2_capacity_class_family <- table(turb_mort_df_ct2$capacity_class,turb_mort_df_ct2$family)

# ACT
report.ACT_capacity_class_family <- ACT_I(ct2_capacity_class_family, alpha=alpha, Rtype='ADJ', nrep=nrep)
report.ACT_capacity_class_family
## Omnibus H of independence rejected --> Alternative H (dependence) accepted; i.e. there is some relatoinship between capacity class and family
report.ACT_capacity_class_family$Famwise_Significant

# Prepare output data or residual analysis for MS
ACT_capacity_class_family_cell_sign <- as.matrix(round(report.ACT_capacity_class_family$Residuals,2))
ACT_capacity_class_family_cell_sign_df_ms <- as.data.frame.matrix(ifelse(report.ACT_capacity_class_family$Famwise_Significant, 
                                                                         paste(ACT_capacity_class_family_cell_sign, "*",sep = ""),
                                                                         paste(ACT_capacity_class_family_cell_sign, "",sep = "")))
ACT_capacity_class_family_expected_freq_df_ms <- as.data.frame.matrix(as.matrix(round(report.ACT_capacity_class_family$ExpectedFrequencies,2)))

# Output dataframe
capacity_class_family_df_ms <- rbind(as.data.frame.matrix(ct2_capacity_class_family),
                                     ACT_capacity_class_family_expected_freq_df_ms,
                                     ACT_capacity_class_family_cell_sign_df_ms)

print(xtable(capacity_class_family_df_ms, type = "latex"), file = "./Results_Figs_Tabs/ACT_capacity_class_family.tex")






