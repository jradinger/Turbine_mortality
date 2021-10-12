####################################################
### Mortality of fish passing turbines
### Helper functions
### Author: Johannes Radinger et al. (see referenced code)
### October 2021
####################################################


#####################################################################
### Helper function to draw values from a distribution describing the mortality range
### Combination of uniform and triangle probability
### mortc =  c(100,250) # first value =  not mort (alive), second = mort
### uncertainty_bounds = c(0.15,0) # first value is potential maximum lower mortality, second is higher mortality
#####################################################################
# Function for a mixed triangular-uniform distribtution
library(triangle)
citation("triangle")
rmix_triangular_uniform <- function (n = 1, a = 0, b = 1, c = (a + b)/2, triangle_weight=0.5){
    #set.seed(iseed)
    p1 <- rtriangle(n=n,a=a,b=b,c=c)
    p2 <- runif(n=n,min=a,max=b)
    #
    flag <- rbinom(n=n,size=1,prob=triangle_weight)
    p <- p1*(flag) + p2*(1-flag)
    #
    return(p)
}

# Calculated redrawn mortality values from a (mixed) triangular-uniform distribtution
mort_uncertainty <- function(n_mort,
                             n_not_mort,
                             uncertainty_lwr,
                             uncertainty_upr,
                             mort_rate=F,
                             triangle_weight=0.5){
  
  n_total = sum(n_mort,n_not_mort)
  modal_mort <- n_mort/n_total
  
  mort_prob <- rmix_triangular_uniform(n = 1, 
                     a = max(0,modal_mort-uncertainty_lwr), 
                     b = min(1,modal_mort+uncertainty_upr), 
                     c = modal_mort,
                     triangle_weight = triangle_weight)
  
  mort_result <- c(round(n_total*mort_prob),round(n_total*(1-mort_prob)))
  
  if(mort_rate){
    return(mort_result[1]/sum(mort_result))
  }else{
  return(mort_result)}
}

# Calculated redrawn length values from a mixed triangular-uniform distribtution
length_uncertainty <- function(avg_length,
                               min_length,
                               max_length,
                               triangle_weight=0.5){
  
  resampled_length <- rmix_triangular_uniform(n = 1, 
                             a = min_length, 
                             b = max_length,
                             c=avg_length,
                             triangle_weight = triangle_weight)
  return(resampled_length)
}


#######################################################
### Helper function for non-parameteric bootstrapping
### Code based on code snippets from the package lmersampler
######################################################
### Function to resample cases within cases of structured data-analysis for bootstrap analysis
### Code based on lmersampler for case.bootstrapping
.cases.within.cases.resamp <- function(dat, 
                                       cluster,
                                       resample, 
                                       resample.mort=FALSE, resample.length=FALSE) {
  
  # exit early for trivial data
  if(nrow(dat) == 1 || all(resample==FALSE))
    return(dat)
  
  # ver <- as.numeric_version(packageVersion("dplyr"))
  res <- dat
  
  # within case-resampling of mortality data
  if(is.list(resample.mort)){
    #check if variables are provided
    if(!all(names(resample.mort) %in% c("n_mort","n_not_mort","uncertainty_lwr","uncertainty_upr","triangle_weight_mort"))){
      stop("Incorrect variables specified in resamp.mort()")
    }
    res[,c(resample.mort$n_mort,resample.mort$n_not_mort)] <- 
      t(apply(res, MARGIN=1, 
              FUN=function(x){mort_uncertainty(n_mort = as.numeric(x[resample.mort$n_mort]),
                                               n_not_mort = as.numeric(x[resample.mort$n_not_mort]),
                                               uncertainty_lwr = as.numeric(x[resample.mort$uncertainty_lwr]),
                                               uncertainty_upr = as.numeric(x[resample.mort$uncertainty_upr]),
                                               triangle_weight = as.numeric(x[resample.mort$triangle_weight_mort]),
                                               mort_rate=F)}))
  }
  # within case-resampling of length data
  if(is.list(resample.length)){
    #check if variables are provided
    if(!all(names(resample.length) %in% c("avg_length","min_length","max_length","triangle_weight_length"))){
      stop("Incorrect variables specified in resamp.length()")
    }
    res[,resample.length$avg_length] <- 
      as.vector(apply(res, MARGIN=1, 
                      FUN=function(x){length_uncertainty(avg_length=as.numeric(x[resample.length$avg_length]),
                                                         min_length=as.numeric(x[resample.length$min_length]),
                                                         max_length=as.numeric(x[resample.length$max_length]),
                                                         triangle_weight=as.numeric(x[resample.length$triangle_weight_length]))}))
  }
  
  for(i in 1:length(cluster)) {
    if(i==1 & resample[i]) {
      dots <- cluster[1]
      grouped <- dplyr::group_by_(res, dots)
      g_rows <- dplyr::group_rows(grouped)
      cls <- sample(seq_along(g_rows), replace = resample[i])
      idx <- unlist(g_rows[cls], recursive = FALSE)
      res <- res[idx, ]
    } else{
      if(i == length(cluster) & resample[i]) {
        dots <- cluster[-i]
        grouped <- dplyr::group_by_(res, .dots = dots)
        res <- dplyr::sample_frac(grouped, size = 1, replace = TRUE)
      } else{
        if(resample[i]) {
          dots <- cluster[i]
          res <- split(res, res[, cluster[1:(i-1)]], drop = TRUE)
          res <- plyr::ldply(res, function(df) {
            grouped <- dplyr::group_by_(df, .dots = dots)
            g_rows <- dplyr::group_rows(grouped)
            # g_rows <- ifelse(ver >= "0.8.0", dplyr::group_rows(grouped), attributes(grouped)$indices)
            cls <- sample(seq_along(g_rows), replace = resample[i])
            idx <- unlist(g_rows[cls], recursive = FALSE)
            grouped[idx, ]
          }, .id = NULL)
        }
      }
    }
    
  }
  return(res)
}

### Function to rerun models based on bootstrapped dataframes
### Code based on lmersampler for case.bootstrapping
.cases.completion <- function(model, rep.data, fn){
  t0 <- fn(model) #what results to extract and use for bootstrap analysis
  
  # Refit the model and apply 'fn' to it using lapply
  form <- model@call$formula
  control_param <- eval(model@call$control)
  family_param <- as.character(model@call$family)
  tstar <- list()
  for(i in (1:length(rep.data))){
    print(i)
    tstar_i <- try(fn(lme4::glmer(formula = form, 
                                  data = rep.data[[i]], 
                                  control = control_param,
                                  family = family_param)),
                   silent=TRUE)
    if(isTRUE(class(tstar_i)=="try-error")) {
      tstar[[i]] <- NULL}
    else{
      tstar[[i]] <- tstar_i}
  }
  tstar <- tstar[lapply(tstar,length)==length(t0)]
  tstar <- do.call("cbind", tstar) # Can these be nested?
  rownames(tstar) <- names(t0)
  
  RES <- structure(list(t0 = t0, 
                        t = t(tstar),
                        R = nrow(t(tstar)), 
                        data = model@frame,
                        seed = .Random.seed, statistic = fn,
                        sim = "case.within.case", call = match.call()),
                   class = "boot")
  
  return(RES)
}

#####################################################################
###  Functions to calculate mean and its CI for bootstrap samples (groups) of a single column
#####################################################################
fn_avg_boot_mort_weighted <- function(x,orig.data,group_var){
  if(!is.na(group_var)){
    t <- sapply(split(x,x[,group_var]),
                FUN=function(x){weighted.mean(x$mort_rate,w=x[,"N"],na.rm=T)})
  }else{
    t <-  weighted.mean(x$mort_rate,w=x[,"N"],na.rm=T)}
  return(t)
  }

bootstrap_avg_mort_weighted <- function(orig.data, rep.data, fn=fn_avg_boot_mort_weighted, group_var=NA){
  t0 <- fn(x=orig.data,orig.data=orig.data,group_var=group_var)
  # Refit stat and apply 'fn'
  tstar <- array(numeric(),c(length(t0),length(rep.data)))
  rownames(tstar) <- names(t0)
  if(is.na(group_var)){rownames(tstar) <- "mean"}
  for(i in (1:length(rep.data))){
    print(i)
    #i=1
    df_i <- as.data.frame(rep.data[[i]])
    tstar_i <- try(
      fn(x=df_i,orig.data=orig.data,group_var=group_var),
      silent=TRUE)
    
    if(isTRUE(class(tstar_i)=="try-error")) {
      tstar[,i]<-NA}
    else{
      if(is.na(group_var)){names(tstar_i) <- "mean"}
      tstar[names(tstar_i),i]<-tstar_i}
  }

  RES <- structure(list(t0 = cbind(t0), 
                        t = t(tstar),
                        R = nrow(t(tstar)), 
                        data = orig.data,
                        seed = .Random.seed, 
                        statistic = "w mean",
                        sim = "case.within.case", call = match.call()),
                   class = "boot")
  return(RES)
}


## Extract CI of bootstrapped dataframe for multiple variables/columns (based on boot.ci)
multi.boot.ci <- function(x, conf=0.95){
  #x=bootstrap_avg_mort_species
  r_ls <- list()
  for(i in c(1:nrow(x$t0))){
    #print(i)
    if(sum(!is.na(x$t[,i]))==0){
      r_ls[[row.names(x$t0)[i]]] <- c(NA,NA,NA)
    }
     else if(var(x$t[,i],na.rm=T)>0){
      r_i <- boot.ci(x, conf=conf, type="perc", index=i, na.rm=TRUE)
      r_ls[[row.names(x$t0)[i]]] <- c(r_i$t0,r_i$percent[c(4,5)])
    }else{
      r_ls[[row.names(x$t0)[i]]] <- c(x$t0[i],x$t0[i],x$t0[i])
    }
  }
  o <- do.call("rbind",r_ls)
  colnames(o) <- c("statistic","conf.low","conf.high")
  as.data.frame(o)
}


#####################################################################
### Dataframe with only complete cases for defined columns
#####################################################################
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec,desiredCols ])
}

#####################################################################
# species-specific body form factor (a 3.0; Froese 2006)
#####################################################################
a3.0 <- function(a,b,S=-1.358){
  10^(log10(a)-(S*(b-3)))}


#####################################################################
### cnt2bin (count to binary) takes a data frame with 2-column ######
### "count" response variable of successes and failures and    ######
### converts it to long format, with one column showing        ######
### 0s and 1s for failures and successes.                      ######
### data is data frame with 2-column response variable         ######
### suc and fail are character expressions for columns         ######
### containing counts of successes and failures respectively   ######
#####################################################################
cnt2bin <- function(data, suc, fail) {
  xvars <- names(data)[names(data)!=suc & names(data)!=fail]
  list <- lapply(xvars, function(z) with(data, rep(get(z), get(suc)+get(fail))))
  names(list) <- xvars
  df <- as.data.frame(list)
  with(data,data.frame(bin=rep(rep(c(1,0),nrow(data)),c(rbind(get(suc),get(fail)))),
                       df))
}

#######################################################################
##### Violin Plots, adapted version with horizontal lines for med, q1 and q3
#####################################################################
vioplot2<-function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
                    horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
                    lwd = 1, ltyq=2, ltymed=1, colMedline="black", colqline="grey",
                    rectCol = "black", colMed = "white", pchMed = 19, 
                    at, add = FALSE, wex = 1, drawRect = TRUE, colPercline="black",
                    cex.axis = 1, las.x.axis = 1)
{
  
  require("sm")
  if(!is.list(x)){
    datas <- list(x, ...)
  } else{
    datas<-x
  }
  #datas=list(rnorm(100))
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  q1_density <- vector(mode = "numeric", length = n)
  q3_density <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  med_density <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if(length(col)!=n){
    col <- rep(col,length.out=n)
  }
  if(length(names)!=n){
    names <- rep(names,length.out=n)
  }
  if(length(border)!=n){
    border <- rep(border,length.out=n)
  }
  if(length(ltyq)!=n){
    ltyq <- rep(ltyq,length.out=n)
  }
  if(length(ltymed)!=n){
    ltymed <- rep(ltymed,length.out=n)
  }
  if(length(colMedline)!=n){
    colMedline <- rep(colMedline,length.out=n)
  }
  if(length(colqline)!=n){
    colqline <- rep(colqline,length.out=n)
  }
  if(length(colPercline)!=n){
    colPercline <- rep(colPercline,length.out=n)
  }
  if(length(rectCol)!=n){
    rectCol <- rep(rectCol,length.out=n)
  }
  if (!(is.null(h))) 
    args <- c(args, h = h)
  for (i in 1:n) {
    #i=1
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    med_density[[i]] <- (do.call("sm.density", c(list(data, xlim = est.xlim,
                                                      eval.points=med[i]),args))$estimate)*hscale
    q1_density[[i]] <- (do.call("sm.density", c(list(data, xlim = est.xlim), 
                                              eval.points=q1[i],args))$estimate)*hscale
    q3_density[[i]] <- (do.call("sm.density", c(list(data, xlim = est.xlim), 
                                              eval.points=q3[i],args))$estimate)*hscale
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if(!add){
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  
  boxwidth <- 0.05 * wex
  
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label, cex.axis = cex.axis, las = las.x.axis)
    }
    box()
    for (i in 1:n) {
      segments(x0 = at[i]-med_density[[i]], x1 = at[i]+med_density[[i]],
               y0 = med[i], y1=med[i], lty=ltymed[i],col=colMedline[i])
      segments(x0 = at[i]-q1_density[[i]], x1 = at[i]+q1_density[[i]],
               y0 = q1[i], y1=q1[i], lty=ltyq[i],col=colqline[i])
      segments(x0 = at[i]-q3_density[[i]], x1 = at[i]+q3_density[[i]],
               y0 = q3[i], y1=q3[i], lty=ltyq[i],col=colqline[i])
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col[i], border = border[i], 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty,col=colPercline[i])
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
             q3[i], col = rectCol[i],border=rectCol[i])
        points(at[i], med[i], pch = pchMed, col = colMed)
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), col = col[i], border = border[i], 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty,col=colPercline[i])
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
               boxwidth/2, col = rectCol[i],border=rectCol[i])
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}

#####################################################################
### overdispersion test http://glmm.wikidot.com/faq
#####################################################################
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}


