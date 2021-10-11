# Turbine mortality #
## Meta-analysis of fish mortality at hydroelectric turbines :fish::+1: R-analysis framework, data and scripts
<img src="https://user-images.githubusercontent.com/2982536/136748303-65b46c91-ffeb-42b7-9d48-1b5da56d0eb5.png" width="650">

=======
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
* Bootstrapping of GLMMs investigating the relationship between morality rate ~ avgerage length, turbine_type etc.
