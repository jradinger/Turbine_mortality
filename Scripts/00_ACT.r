####################################################
### Mortality of fish passing turbines
### Analysis of residuals in contingency tables
### Author: (see referenced code)
### October 2021
####################################################


### This script is from the ACT software for Analysis of Contingency Tables presented in García-
### Pérez, Núñez-Antón, and Alcalá-Quintana (“Analysis of residuals in contingency tables: Another nail in
### the coffin of conditional approaches to significance testing,” Behavior Research Methods, 2014,
### http://dx.doi.org/10.3758/s13428-014-0472-0).

####################
### Independence ###
####################
ACT_I <- function(observed, alpha, Rtype, nrep){
  ifelse (is.character(Rtype), Rtype <- toupper(Rtype), stop('Error: Invalid Rtype (must be a string)', call.=FALSE)) 
  observed <- data.matrix(observed)
  if (length(dim(observed))!=2 || length(observed) < 4) {
    stop('Error: Not a two-way table; please fix and rerun)', call.=FALSE)
    } else if (!is.numeric(observed)){
        stop('Error: Non-numeric table; please fix and rerun)', call.=FALSE)
    } else if (!is.numeric(alpha) || alpha<=0 || alpha >=0.5 || length(alpha)!=1) {
        stop('Error: Invalid alpha (0 < alpha < 0.5)', call.=FALSE)
    } else if(!any (identical(Rtype, 'ADJ'), identical(Rtype, 'MC'))){
        stop('Error: Invalid Rtype (must be "MC" for "moment corrected residuals"
               or "ADJ" for "adjusted residuals"; case insensitive)', call.=FALSE)
    } else if (!is.numeric(nrep) || !is.wholenumber(nrep) || nrep<=0 || length(nrep)!=1) {
        stop('Error: Invalid nrep (must be a positve scalar)', call.=FALSE)
    } else if (any(is.na(observed))) {
        stop ('Error: Table cannot contain missing values; please fix and rerun', call.=FALSE)
    }
  ### Computes expected frequencies, probabilities, and residuals for the data entered
  nfil <- nrow(observed); ncolumn <- ncol(observed)
  margin.col <- t(colSums(observed))
  margin.row <- t(t(rowSums(observed)))
  if (any(margin.row==0) || any (margin.col==0)) {
   stop('Error: table cannot have empty rows or columns; please fix and rerun', call.=FALSE)
  } 
  n <- sum(observed)
  expected <- margin.row%*%margin.col/n
  probabilities <- expected/n
  ifelse (Rtype=='ADJ',
          variances <- (1-margin.row/n)%*%(1-margin.col/n),
          variances <- matrix((ncolumn-1)*(nfil-1)/(ncolumn*nfil),nfil,ncolumn))
  residuals <- (observed-expected)/sqrt(expected*variances) 
  ### Simulates
  simulated <- array(rmultinom(nrep,n,probabilities),c(nfil,ncolumn,nrep))
  sim.expected <- array (NA,c(nfil,ncolumn,nrep))
  sim.variances <- array (NA,c(nfil,ncolumn,nrep))
  for (i in 1:nrep){
    sim.expected[,,i] <- t(t(rowSums(simulated[,,i])))%*%t(colSums(simulated[,,i]))/n
  }
  if (Rtype=='ADJ'){
    for (i in 1:nrep){      
      sim.variances[,,i] <- (1-t(t(rowSums(simulated[,,i])))/n)%*%(1-t(colSums(simulated[,,i]))/n)
    }      
  } else {sim.variances <- array((ncolumn-1)*(nfil-1)/(ncolumn*nfil),c(nfil,ncolumn,nrep))}
  sim.residuals <- (simulated-sim.expected)/sqrt(sim.expected*sim.variances)
  toKeep <- which(is.finite(apply(sim.residuals,3,sum)))
  valid <- length(toKeep)
  if (valid==0) {
    stop('Table is too sparse to produce valid replicates; consider merging rows or columns', call.=FALSE)
  } else if (valid <= nrep/2) {
    warning ('Table seems to be too sparse; consider merging rows or columns', call.=FALSE)
  }
  sim.residuals <- sim.residuals[,,toKeep]
  total <- length(sim.residuals)
  zmin <- qnorm(1-alpha); zmax <- 10
  for (i in 1:25){
    z.omnibus <- (zmax+zmin)/2
    type1 <- mean(apply(abs(sim.residuals)>z.omnibus,3,any))
    ifelse (type1>alpha, zmin <- z.omnibus, zmax <- z.omnibus)
  }
  z.residuals <- qnorm(1-alpha/2)
  type1.cell <- mean(abs(sim.residuals)>z.residuals)
  alpha.star <- 2*pnorm(-z.omnibus)
  signif.omnibus <- (abs(residuals)>z.omnibus)
  signif.residual <- (abs(residuals)>z.residuals)
  ifelse (any(signif.omnibus),omnibus.test <- 'Rejected', omnibus.test <- 'Not rejected')
  ### Output
  stri <- switch (Rtype,
                 'MC'='moment-corrected',
                 'ADJ'='adjusted')
  if (total > 1000){
    report <- list(Problem = paste('Omnibus test of independence and'
                                   ,stri,'residual analysis'),
                   InputTable = observed,
                   NominalTestSize = alpha,
                   NumReplicates = nrep,
                   ValidReplicates = valid,
                   ExpectedFrequencies = expected,
                   Residuals = residuals,
                   Cellwise_CriticalValue = z.residuals,
                   Cellwise_Significant = signif.residual,
                   Cellwise_ExactTestSize = type1.cell,                
                   Famwise_AlphaStar = alpha.star,
                   Famwise_CriticalValue = z.omnibus,
                   Famwise_Significant = signif.omnibus,
                   Famwise_ExactTestSize = type1,
                   OmnibusHypothesis = omnibus.test)
  } else {
      report <- list(Problem = paste('Omnibus test of independence and'
                                     ,stri,'residual analysis'),
                   InputTable = observed,
                   NominalTestSize = alpha,
                   NumReplicates = nrep,
                   ValidReplicates = valid,
                   ExpectedFrequencies = expected,
                   Residuals = residuals,
                   Cellwise_CriticalValue = z.residuals,
                   Cellwise_Significant = signif.residual,
                   Cellwise_ExactTestSize ='Not computed; insufficent valid replicates',                
                   Famwise_AlphaStar = 'Not computed; insufficent valid replicates',
                   Famwise_CriticalValue = 'Not computed; insufficent valid replicates',
                   Famwise_Significant = 'Not computed; insufficent valid replicates',
                   Famwise_ExactTestSize = 'Not computed; insufficent valid replicates',
                   OmnibusHypothesis = 'Not conducted; insufficent valid replicates')
  }
  class(report) <- 'ACT'
  return(report)}

###################
### Homogeneity ###
###################
ACT_H <- function(observed, alpha, Rtype, nrep){
  ###Check conditions
  ifelse (is.character(Rtype), Rtype <- toupper(Rtype), stop('Error: Invalid Rtype (must be a string)', call.=FALSE)) 
  observed <- data.matrix(observed)
  if (length(dim(observed))!=2 || length(observed) < 4) {
    stop('Error: Not a two-way table; please fix and rerun)', call.=FALSE)
    } else if (!is.numeric(observed)){
  stop('Error: Non-numeric table; please fix and rerun)', call.=FALSE)
  } else if (!is.numeric(alpha) || alpha<=0 || alpha >=0.5 || length(alpha)!=1) {
      stop('Error: Invalid alpha (0 < alpha < 0.5)', call.=FALSE)
  } else if(!any (identical(Rtype, 'ADJ'), identical(Rtype, 'MC'))){
      stop('Error: Invalid Rtype (must be "MC" for "moment corrected residuals" or
           "ADJ" for "adjusted residuals"; case insensitive)', call.=FALSE)
  } else if (!is.numeric(nrep) || !is.wholenumber(nrep) || nrep<=0 || length(nrep)!=1) {
      stop('Error: Invalid nrep (must be a positve scalar)', call.=FALSE)
  } else if (any(is.na(observed))) {
      stop ('Error; Table cannot contain missing values; please fix and rerun', call.=FALSE)
  }
  ### Computes expected frequencies, probabilities, and residuals for the data entered
  nfil <- nrow(observed); ncolumn <- ncol(observed)
  margin.col <- t(colSums(observed))
  margin.row <- t(t(rowSums(observed)))
  if (any(margin.row==0) || any (margin.col==0)) {
    stop('Error: table cannot have empty rows or columns; please fix and rerun', call.=FALSE)
  } 
  n <- sum(observed)
  expected <- margin.row%*%margin.col/n
  probabilities <- margin.col/n
  ifelse (Rtype=='ADJ',
          variances <- (1-margin.row/n)%*%(1-margin.col/n),
          variances <- matrix((ncolumn-1)*(nfil-1)/(ncolumn*nfil),nfil,ncolumn))
  residuals <- (observed-expected)/sqrt(expected*variances) 
  ### Simulates
  simulated <- array (NA,c(nfil,ncolumn,nrep))
  for (i in 1:nfil){
    simulated[i,,] <- rmultinom(nrep,margin.row[i,1],probabilities)
  }
  sim.expected <- array (NA,c(nfil,ncolumn,nrep))
  sim.variances <- array (NA,c(nfil,ncolumn,nrep))
  for (i in 1:nrep){
    sim.expected[,,i] <- t(t(rowSums(simulated[,,i])))%*%t(colSums(simulated[,,i]))/n
  }
  if (Rtype=='ADJ'){
    for (i in 1:nrep){      
      sim.variances[,,i] <- (1-t(t(rowSums(simulated[,,i])))/n)%*%(1-t(colSums(simulated[,,i]))/n)
    }      
  } else {sim.variances <- array((ncolumn-1)*(nfil-1)/(ncolumn*nfil),c(nfil,ncolumn,nrep))}
  sim.residuals <- (simulated-sim.expected)/sqrt(sim.expected*sim.variances)
  toKeep <- which(is.finite(apply(sim.residuals,3,sum)))
  valid <- length(toKeep)
  if (valid==0) {
    stop('Table is too sparse to produce valid replicates; consider merging rows or columns', call.=FALSE)
  } else if (valid <= nrep/2) {
    warning ('Table seems to be too sparse; consider merging rows or columns', call.=FALSE)
  }
  sim.residuals <- sim.residuals[,,toKeep]
  total <- length(sim.residuals)
  zmin <- qnorm(1-alpha); zmax <- 10
  for (i in 1:25){
    z.omnibus <- (zmax+zmin)/2
    type1 <- mean(apply(abs(sim.residuals)>z.omnibus,3,any))
    ifelse (type1>alpha, zmin <- z.omnibus, zmax <- z.omnibus)
  }
  z.residuals <- qnorm(1-alpha/2)
  type1.cell <- mean(abs(sim.residuals)>z.residuals)
  alpha.star <- 2*pnorm(-z.omnibus)
  signif.omnibus <- (abs(residuals)>z.omnibus)
  signif.residual <- (abs(residuals)>z.residuals)
  ifelse (any(signif.omnibus),omnibus.test <- 'Rejected', omnibus.test <- 'Not rejected')
  ### Output
  stri <- switch (Rtype,
                  'MC'='moment-corrected',
                  'ADJ'='adjusted')
if (total > 1000){
  report <- list(Problem = paste('Omnibus test of homogeneity and',stri,'residual analysis'),
                 InputTable = observed,
                 NominalTestSize = alpha,
                 NumReplicates = nrep,
                 ValidReplicates = valid,
                 ExpectedFrequencies = expected,
                 Residuals = residuals,
                 Cellwise_CriticalValue = z.residuals,
                 Cellwise_Significant = signif.residual,
                 Cellwise_ExactTestSize = type1.cell,                
                 Famwise_AlphaStar = alpha.star,
                 Famwise_CriticalValue = z.omnibus,
                 Famwise_Significant = signif.omnibus,
                 Famwise_ExactTestSize = type1,
                 OmnibusHypothesis = omnibus.test)
  } else {
    report <- list(Problem = paste('Omnibus test of homogeneity and',stri,'residual analysis'),
                   InputTable = observed,
                   NominalTestSize = alpha,
                   NumReplicates = nrep,
                   ValidReplicates = valid,
                   ExpectedFrequencies = expected,
                   Residuals = residuals,
                   Cellwise_CriticalValue = z.residuals,
                   Cellwise_Significant = signif.residual,
                   Cellwise_ExactTestSize ='Not computed; insufficent valid replicates',                
                   Famwise_AlphaStar = 'Not computed; insufficent valid replicates',
                   Famwise_CriticalValue = 'Not computed; insufficent valid replicates',
                   Famwise_Significant = 'Not computed; insufficent valid replicates',
                   Famwise_ExactTestSize = 'Not computed; insufficent valid replicates',
                   OmnibusHypothesis = 'Not conducted; insufficent valid replicates')
                 }
class(report) <- 'ACT'
return(report)}

###########
### Fit ###
###########
ACT_F <- function(observed, model, alpha, Rtype, nrep){
  ### Check conditions
  ifelse (is.character(Rtype), Rtype <- toupper(Rtype), stop('Error: Invalid Rtype (must be a string)', call.=FALSE)) 
  observed <- data.matrix(observed)
  if (length(dim(observed))!=2 || length(observed) < 4) {
    stop('Error: Not a two-way table; please fix and rerun)', call.=FALSE)
    } else if (!is.numeric(observed)){
      stop('Error: Non-numeric table; please fix and rerun)', call.=FALSE)
    } else if (!identical(dim(observed), dim(model))){
      stop('Error: Model does not much table; please fix and rerun)', call.=FALSE)
    } else if (!is.numeric(alpha) ||alpha<=0 || alpha >=0.5|| length(alpha)!=1) {
      stop('Error: Invalid alpha (0 < alpha < 0.5)', call.=FALSE)
    } else if(!any (identical(Rtype, 'ADJ'), identical(Rtype, 'MC'))){
      stop('Error: Invalid Rtype (must be "MC" for "moment corrected residuals" or
                "ADJ" for "adjusted residuals"; case insensitive)', call.=FALSE)
    } else if (!is.numeric(nrep) || !is.wholenumber(nrep) || nrep<=0 || length(nrep)!=1) {
      stop('Error: Invalid nrep (must be a positve scalar)', call.=FALSE)
    } else if (any(is.na(observed))) {
      stop ('Error; Table cannot contain missing values; please fix and rerun', call.=FALSE)
    }
  ### Computes expected frequencies, probabilities, and residuals for the data entered
  nfil <- nrow(observed); ncolumn <- ncol(observed)
  margin.col <- t(colSums(observed))
  model <- model/matrix(colSums(model), nfil, ncolumn, byrow=T)
  if (any(model==0)) {
    warning('Model assigns zero probability to one or more cells. One or more residuals will be Inf or NaN, and will not be considered', call.=FALSE)
  }
  validrows <- colSums(model>0)
  expected <- model*matrix(margin.col, nrow=nfil, ncol=ncolumn, byrow=T)
  ifelse (Rtype=='ADJ',
          variances <- 1-model,
          variances <- matrix((validrows-1)/validrows, nfil, ncolumn, byrow=T))
  residuals <- (observed-expected)/sqrt(expected*variances)
  simulated <- array (NA,c(nfil,ncolumn,nrep))
  ### Simulates
  for (i in 1:ncolumn){
    simulated[,i,] <- rmultinom(nrep,margin.col[1,i],model[,i])
  }
  sim.expected <- array (expected,c(nfil,ncolumn,nrep))
  sim.residuals <- (simulated-sim.expected)/sqrt(sim.expected*array(variances,c(nfil,ncolumn,nrep)))
  toKeep <- apply(sim.residuals,c(1,2,3), is.finite)
  sim.residuals[!toKeep] <- NA
  zmin <- qnorm(1-alpha); zmax <- 10
  for (i in 1:25){
    z.omnibus <- (zmax+zmin)/2
    type1 <- mean(apply(abs(sim.residuals)>z.omnibus,3,any, na.rm=TRUE), na.rm=TRUE)
    ifelse (type1>alpha, zmin <- z.omnibus, zmax <- z.omnibus)
  }
  z.residuals <- qnorm(1-alpha/2)
  type1.cell <- mean(abs(sim.residuals)>z.residuals, na.rm=TRUE)
  alpha.star <- 2*pnorm(-z.omnibus)
  signif.omnibus <- (abs(residuals)>z.omnibus)
  signif.omnibus[!is.finite(residuals)] <- FALSE
  signif.residual <- (abs(residuals)>z.residuals)
  signif.residual[!is.finite(residuals)] <- FALSE
  ifelse (any(signif.omnibus),omnibus.test <- 'Rejected', omnibus.test <- 'Not rejected')
  ### Output
  stri <- switch (Rtype,
                  'MC'='moment-corrected',
                  'ADJ'='adjusted')
  report <- list(Problem = paste('Omnibus test of fit and',stri,'residual analysis'),
                 InputTable = observed,
                 Model = model,
                 NominalTestSize = alpha,
                 NumReplicates = nrep,
                 ExpectedFrequencies = expected,
                 Residuals = residuals,
                 Cellwise_CriticalValue = z.residuals,
                 Cellwise_Significant = signif.residual,
                 Cellwise_ExactTestSize = type1.cell,
                 Famwise_AlphaStar = alpha.star,
                 Famwise_CriticalValue = z.omnibus,
                 Famwise_Significant = signif.omnibus,
                 Famwise_ExactTestSize = type1,
                 OmnibusHypothesis = omnibus.test)
  class(report) <- 'ACT'
  return(report)}

###################
### wholenumber ###
###################
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

###############################
### Format & print matrices ###
###############################
print.ACT.matrix <- function(m){
  write.table(format(m, justify="right"),
              row.names=T, col.names=F, sep="  ", quote=F)
}

############################
### Print method for ACT ###
############################
print.ACT <- function (x){
  nrows <- nrow(x$InputTable)
  cat(rep ('-', times=nchar(x$Problem)+4), sep="", '\n') 
  cat(' ',x$Problem, '\n')
  cat(rep ('-', times=nchar(x$Problem)+4), sep="", '\n') 
  cat('             Input table:',format(x$InputTable[1,], trim=FALSE, big.interval="  "),'\n')
  for (i in 2:nrows){
    cat('                         ',format(x$InputTable[i,], trim=FALSE),'\n')
  } 
  cat('\n')
  if(any(grep('fit', x$Problem))){
    row.names(x$Model) <- 
      c('                   Model:',rep('                         ', times=nrows-1))
    print.ACT.matrix (x$Model)
    cat('\n')
  }
  cat('       Nominal test size:', x$NominalTestSize, '\n')
  cat('    Number of replicates:', x$NumReplicates, '\n')
  if(!any(grep('fit', x$Problem))) cat('        Valid replicates:', x$NumReplicates, '\n')
  cat('\n')
  row.names(x$ExpectedFrequencies) <- 
    c('    Expected frequencies:',rep('                         ', times=nrows-1))
  print.ACT.matrix (x$ExpectedFrequencies)
  cat('\n')
  row.names(x$Residuals) <- 
    c('               Residuals:',rep('                         ', times=nrows-1))
  print.ACT.matrix (x$Residuals)
  cat('\n')
  cat(' Cellwise critical value:', x$Cellwise_CriticalValue, '\n')
  cat('\n')
  cat('    Cellwise significant:', format (x$Cellwise_Significant[1,]), '\n')
  for (i in 2:nrows){
    cat('                         ',format (x$Cellwise_Significant[i,]),'\n')
  } 
  cat('\n')
  cat('Cellwise exact test size:', x$Cellwise_ExactTestSize, '\n')
  cat('   Familywise alpha star:', x$Famwise_AlphaStar, '\n')
  cat('  Famwise critical value:', x$Famwise_CriticalValue, '\n') 
  if (length(dim(x$Famwise_Significant)) == 2){
    cat('\n')
    cat('  Familywise significant:', format(x$Famwise_Significant[1,]), '\n')
    for (i in 2:nrows){
      cat('                         ',format(x$Famwise_Significant[i,]),'\n')
    }
    cat('\n')
  } else {cat('  Familywise significant:', x$Famwise_Significant, '\n')
  }
  cat(' Famwise exact test size:', x$Famwise_ExactTestSize, '\n')
  cat('      Omnibus hypothesis:', x$OmnibusHypothesis, '\n')
  invisible(x)
}

