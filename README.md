# Turbine mortality #
## Meta-analysis of fish mortality at hydroelectric turbines :fish::+1: R-analysis framework, data and scripts
<img src="https://user-images.githubusercontent.com/2982536/136748303-65b46c91-ffeb-42b7-9d48-1b5da56d0eb5.png" width="650">



## Content - Scripts
### 01 Descriptive statistics
* Alluvial /Sankey diagrams to visually describe observed relationships between taxonomic order, hydropower scale, turbine types and lethal vs. sub-lethal effects on fish
* Summary (e.g. mean mortality) by species, taxonomic order

### 02 Violin plot 
* Violin plots showing distributions of mortality estimates across species/taxonomic orders

### 03 Analysis of contingency tables
* Analysis based on ACT method (Analysis of Contingency Tables, García-Pérez et al., 2014)
* Capactiy class vs. turbine type
* Turbine type vs. family
* Capacity class vs. family

### 04 Geographical distribution
* Map plot of investigated turbine mortality studies across space

### 05 Case bootstrapping of mean mortatilty
* Calculate bootstrapped confidence intervall of turbine-related mortality based on individual counts following case bootstrapping
  
### 06 Generalized linear mixed model
* Bootstrapping of GLMMs investigating the relationship between morality rate - avgerage length, turbine_type etc.


---
## Content - Data
### Dataset: turb_df.csv

| Column  | Description |
| ---------- | ---------- |
| study_ID  | ID of each study  |
| study  | specific study |
| experiment  | specific experiment within study |
| first_author  | first author of study |
| year  | year of study |
| n_turbines  | number of turbines at HPP location |
| n_turbines_studied  | number of turbines studied at a given location |
| total_capacity_study_turbines_kW  | total generating capacity of all studied turbines at HPP location |
| total_capacity_kW  | total generating HPP capacity at location |
| capacity_class  | HPP generating capacity class |
| turbine_type  | type of turbine |
| max_diameter  | maximum diameter of turbine |
| n_blades  | number of turbine blades |
| rpm  | revolutions per minute |
| head  | appr. (hydraulic) head HPP |
| max_flow_rate_HPP  | total max. flow rate at HPP |
| sampling_method  | sampling method |
| country  | country |
| global_region  | continent |
| location  | specific location |
| waterbody  | name of waterbody/river |
| MQ_waterbody  | mean discharge study river |
| Lat  | Latitude in decimal degrees |
| Long  | Longitude in decimal degrees |

### Dataset: turb_mort_df.csv

| Column  | Description |
| ---------- | ---------- |
| species_code  | species code |
| study_ID  | ID of each study  |
| study  | specific study |
| experiment  | specific experiment within study |
| first_author  | first author of study |
| year  | year of study |
| taxorder  | taxonomic order |
| family  | taxonomic family |
| genus  | taxonomic genus |
| species  | scientific species name |
| common_name  | common species name |
| min_length  | min. of studied range of body lengths |
| max_length  | max. of studied range of body lengths |
| avg_length  | average of studied range of body lengths |
| N  | total number of individuals studied |
| n_mort  | number of  dead and medium to severely injured fish |
| n_not_mort  | number of alive or sublethally injured fish |
| mort_rate  | mortality rate |
| mort_rate_lwr  | lower estimate of mortality rate |
| mort_rate_upr  | upper estimate of mortality rate |
| post_monitoring_uncertainty  | assigned uncertainty related to post-monitoring |
| catch_related_uncertainty  | assigned uncertainty related to catch/handling |
| bshape  | body shape |
| a3.0  | body shape parameter *) |
| swim_bladder  | type of swim bladder |
| rel_height  | relative body height |
| rel_width  | relative body width |
| n_turbines  | number of turbines at HPP location |
| n_turbines_studied  | number of turbines studied at a given location |
| total_capacity_study_turbines_kW  | total generating capacity of all studied turbines at HPP location |
| total_capacity_kW  | total generating HPP capacity at location |
| capacity_class  | HPP generating capacity class |
| turbine_type  | type of turbine |
| max_diameter  | maximum diameter of turbine |
| n_blades  | number of turbine blades |
| rpm  | revolutions per minute |
| head  | appr. (hydraulic) head HPP |
| max_flow_rate_HPP  | total max. flow rate at HPP |
| sampling_method  | sampling method |
| global_region  | continent |
| country  | country |
| location  | specific location |
| waterbody  | name of waterbody/river |
| MQ_waterbody  | mean discharge study river |
| Lat  | Latitude in decimal degrees |
| Long  | Longitude in decimal degrees |

*) body shape parameter calculated from: 
```R
a3.0 <- apply(turb_mort_df, MARGIN=1, FUN=function(x){a3.0(a=as.numeric(x["a_bayes"]),
                                                b=as.numeric(x["b_bayes"]),
                                                S=-1.358)})
```
