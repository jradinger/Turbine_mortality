##############################################
### Mortality of fish passing turbines
### Descriptive statistics of turbine mortality dataset
### Taxa violin plots
### Author: Johannes Radinger
### October 2021
###############################################

## Load libraries
library(tinter)
library(see)
library(vioplot)
library(dplyr)
library(forcats)

###############################################
# Load data
source("./Scripts/00_Helper_functions.R")
turb_mort_df <- read.csv("./Data/turb_mort_df.csv",stringsAsFactors=T)

# handling corrected or all cases
handling_label <- ""
handling <- F
if(handling){
  print("handling corrected data used instead of complete data")
  handling_label <- "_handling_corrected"
  turb_mort_df <- turb_mort_df[turb_mort_df$catch_related_uncertainty==0 & turb_mort_df$post_monitoring_uncertainty==0,]
}

###############################################
# triangle probability
turb_mort_df$triangle_weight_mort <- 1

# Aggregate individuals of less studied orders into order "others"
turb_mort_df$taxorder_agg <- turb_mort_df$taxorder
sum_taxorder <- aggregate(turb_mort_df$N,by=list(taxorder = turb_mort_df$taxorder),FUN=sum)
turb_mort_df$taxorder_agg <- fct_collapse(turb_mort_df$taxorder_agg,
                                          "Others" = as.character(sum_taxorder$taxorder[sum_taxorder$x<5000]))

#############################################
# Generate bootstrapped replicate dataframes
#############################################
B_violin <- 100 # number of replicates
rep.data.violin <- list()
i <- 1
while(i <= B_violin){
  print(i)
  rep.data.violin_i <- try(.cases.within.cases.resamp(dat=turb_mort_df,
                                                      cluster = c("experiment",".id"),
                                                      resample = c(TRUE,TRUE),
                                                      resample.mort = list(n_mort = "n_mort",
                                                                           n_not_mort = "n_not_mort",
                                                                           uncertainty_lwr = "catch_related_uncertainty",
                                                                           uncertainty_upr = "post_monitoring_uncertainty",
                                                                           triangle_weight_mort = "triangle_weight_mort"),
                                                      resample.length = FALSE),
                           silent=TRUE)
  if(isTRUE(class(rep.data.violin_i)=="try-error")) {
    next}else{
      rep.data.violin_i$mort_rate_resampled <- (rep.data.violin_i$n_mort/rep.data.violin_i$N)
      rep.data.violin[[i]] <- rep.data.violin_i
      i <- i+1}
}


########################################################
# Select taxa for violin plots
#########################################################
# Get list of top most studied species/taxonomic orders
tax_rank <- "taxorder"
#tax_rank <- "species"
if(tax_rank=="taxorder"){
  no_parts <- 1
  selected_taxa <- levels(turb_mort_df$taxorder_agg)
  tax_col <- "taxorder_agg"
}
if(tax_rank=="species"){
  no_parts <- 2
  selected_sp <- levels(turb_mort_df$species_code)[table(turb_mort_df$species_code)>10]
  species_lookup <- unique(turb_mort_df[c("species_code", "species")])
  selected_taxa <- as.vector(species_lookup$species[match(selected_sp,species_lookup$species_code)])
  tax_col <- "species"
  
}

# Get sortings of taxa (descending by mortality)
taxa_sort <- turb_mort_df %>%
  group_by(taxon=get(tax_col)) %>%
  summarise_at("mort_rate", funs(median, mean)) %>%
  arrange(desc(mean)) %>%
  arrange(desc(median)) %>%
  filter(taxon %in% selected_taxa)
taxa_sort_df <- as.data.frame(taxa_sort)
selected_taxa_sort_complete <- taxa_sort_df$taxon
selected_taxa_sort_complete <- factor(selected_taxa_sort_complete,
                                      levels=as.character(selected_taxa_sort_complete))

# Split data into parts /  Make subsets if needed
selected_taxa_sort_parts <- split(selected_taxa_sort_complete,
      sort(rep_len(c(1:no_parts),length.out=length(selected_taxa_sort_complete))))

for(j in 1:no_parts){
  #j=2
  selected_taxa_sort <- selected_taxa_sort_parts[[j]]
  
  # Select subset of bootstrapped dataframe with only selected taxa
  rep.data.violin_sub <- lapply(rep.data.violin,
                                FUN=function(x){
                                  x <- as.data.frame(x)
                                  x <- droplevels(x[x[,tax_col] %in% selected_taxa_sort,])
                                  x[,tax_col] <- factor(x[,tax_col],
                                                        levels=levels(selected_taxa_sort))
                                  x <- droplevels(x)
                                  x
                                }
  )
  
  # Select subset of original dataframe with only selected taxa
  orig.data.violin_sub  <- droplevels(turb_mort_df[turb_mort_df[,tax_col] %in% selected_taxa_sort,])
  orig.data.violin_sub[,tax_col] <- factor(orig.data.violin_sub[,tax_col],
                                           levels=levels(selected_taxa_sort))
  orig.data.violin_sub <- droplevels(orig.data.violin_sub)
  
  ########################################################
  # Build violin plots
  ########################################################
  # Get redrawn mortality estimates (redrawn from triangular uncertainty distribution around measure value)
  # and display as overlaid violin plots
  cairo_pdf(paste("./Results_Figs_Tabs/Violins_",tax_rank,handling_label,"_part",j,".pdf",sep=""),
            width=8,height=5)
  par(mar=c(9,4,0.5,0.5))
  B_violin <- 100
  for(i in seq(1,B_violin+1)){
    print(i)
    if(i==1){
      add=F
      plot.new()}else{add=T}
    
    # Overlay of original data after all bootstrapped violins are printed out
    if(i==B_violin+1){
      # Split across species
      Mort_list <- split(orig.data.violin_sub$mort_rate,orig.data.violin_sub[,tax_col])
      # Plot average values
      vioplot2(Mort_list,
               names=NULL,
               col=NA,
               colPercline=NA,
               rectCol=material_colors("teal"),
               border="#00645a",
               lwd=2,
               colMed=NA,
               wex=0.8,
               ylim=c(0,1),
               colMedline=material_colors("teal"),
               colqline=NA,
               #wex=0.8,
               add=add
      )
      mtext("Mortality rate",side=2,line=2)
      
    }else{
      
      
      #Select bootstrap sample i
      subset_turb_mort_df_resampled_i <- rep.data.violin_sub[[i]]
      
      # Split across species
      Mort_list_resampled <- split(subset_turb_mort_df_resampled_i$mort_rate_resampled,
                                   subset_turb_mort_df_resampled_i[,tax_col])
      
      # Plot resampled values
      if(tax_rank=="taxorder"){
        las.x.axis=1
      }
      if(tax_rank=="species"){
        las.x.axis=2
      }
      vioplot2(Mort_list_resampled, 
               col=NA,
               border=adjustcolor("grey80", alpha.f = 0.55),
               colPercline=NA,
               colMed=NA,
               wex=0.8,
               colMedline=adjustcolor("grey80", alpha.f = 0.55),
               colqline=NA,
               rectCol=NA,
               add=add,
               ylim=c(0,1),
               names=names(Mort_list_resampled),
               cex.axis=0.7,
               las.x.axis=las.x.axis)
    }
  }
  dev.off()
}
