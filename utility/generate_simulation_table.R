#######################################    Comparison of point prediction      ###########################
#Predict conditional mean, SD, quantiles
#s:  number of points for calculating MSPE
#fit.gcde: the gcde fitted model for evaluating
#x.test: Test data. If NULL, generate some new data.
#y.test: Test data. If NULL, generate some new daat.
#fit.gcde.extra: If non null, more gcde fitted model can be compared together.
#extra.names: The name of extra gcde fitted model.
#otherparamter: The parameter for KDE, NNKCDE, FLEXCODE.   If NULL, these methods will NOT be compared.
#clnum: How many nodes are used for parallelization. If 1, no parallelization.
#m.of.g: number of observations calculated by GCDE to get the point estimate(expectation).

generate.table <- function(s = 1000, fit.gcde, x.test = NULL, y.test = NULL, 
                           fit.gcde.extra = NULL, 
                           extra.names = NULL,
                           otherparameter = NULL, clnum = NULL, m.of.g = 10000){
    #If data for generating data is given, then just use it.
    if(!is.null(x.test)) {
        s <- nrow(x.test)
        random.x.all <- x.test
    } else {
        #if x.test is not provided, create them.
        if(fit.gcde$diagnostic.info$x.dist != "normal") 
            stop("If the distribution of X is not normal, please provide it to function 'generate.table'.")
        set.seed(155)
        random.x.all <- matrix(rnorm(fit.gcde$dim.x * s), nrow = s, ncol = fit.gcde$dim.x)}
    
    true.point <- matrix(0, nrow = s, ncol = 7)
    dmr.point <- matrix(0, nrow = s, ncol = 7)
    
    extra.gcde.length <- length(fit.gcde.extra)
    if(!is.null(fit.gcde.extra)){
        if(is.null(extra.names)) stop("please provide the names for all GCDE fit.")
        extra.dmr.point.list <- list()
        extra.t.list         <- list()
        for (i in 1:extra.gcde.length) {
            extra.dmr.point.list[[i]] <- dmr.point
        }
    }
    
    nncde.point <- matrix(0, nrow = s, ncol = 7)
    ckde.point <- matrix(0, nrow = s, ncol = 7)
    flex.point <- matrix(0, nrow = s, ncol = 7)
    quantile.proba <- c(0.05, 0.25, 0.5, 0.75, 0.95)
    
    
    y.true.range <- matrix(0, nrow = s, ncol = 2)
    message("generate table:GCDE...")
    #Evaluate performace for GCDE fit
    for (i in 1:s) {
        random.x <- random.x.all[i, ]
        #generate true data
        if(!is.null(fit.gcde$diagnostic.info)){
                diag.data <- data.generate(n = 5000, x = random.x,
                                           dim.x = fit.gcde$dim.x, dim.y = fit.gcde$dim.y, 
                                           dim.epsilon = fit.gcde$diagnostic.info$dim.epsilon,
                                           epsilon.dist = fit.gcde$diagnostic.info$epsilon.dist,
                                           model = fit.gcde$diagnostic.info$model,
                                           seed = NULL,
                                           testdata = FALSE,
                                           extrapara = fit.gcde$diagnostic.info$extrapara)
            y.true <- diag.data$y
        } else y.true <- NULL
        #extrat x matrix from true data
        x.matrix <- matrix(random.x, nrow = m.of.g, ncol = fit.gcde$dim.x, byrow = TRUE)
        
        if(fit.gcde$eta.dist == "normal"){
            eta <- matrix(rnorm(fit.gcde$dim.eta * m.of.g), nrow = m.of.g, ncol = fit.gcde$dim.eta)
        }
        if(fit.gcde$eta.dist == "uniform"){
            eta <- matrix(runif(fit.gcde$dim.eta * m.of.g,-1,1), nrow = m.of.g, ncol = fit.gcde$dim.eta)
        }
        
        t <- as.array(predict(fit.gcde$generator.t,list(x.matrix, eta)))
        if(!is.null(fit.gcde.extra)){
            for (j in 1:extra.gcde.length) {
                if(fit.gcde.extra[[j]]$dim.eta != fit.gcde$dim.eta){
                    if(fit.gcde.extra[[j]]$eta.dist == "normal"){
                        eta <- matrix(rnorm(fit.gcde.extra[[j]]$dim.eta * m.of.g), nrow = m.of.g, 
                                      ncol = fit.gcde.extra[[j]]$dim.eta)
                    }
                    if(fit.gcde.extra[[j]]$eta.dist == "uniform"){
                        eta <- matrix(runif(fit.gcde.extra[[j]]$dim.eta * m.of.g,-1,1), nrow = m.of.g, 
                                      ncol = fit.gcde.extra[[j]]$dim.eta)
                    }
                }
                extra.t.list[[j]] <- as.array(predict(fit.gcde.extra[[j]]$generator.t,list(x.matrix, eta)))
            }
        }
        #Store the range of true y value for use in next part: other methods
        y.true.range[i, ] <- range(y.true)
        #mean
        true.point[i,1] <- diag.data$y.mean[1,]
        #sd
        true.point[i,2] <- diag.data$y.sd[1]
        #quantiles
        true.point[i,3:7] <- quantile(y.true, quantile.proba)
        #single new y value, use just the first one if the test data is not provided.
        #if(is.null(y.test)) {true.point[i,8] <- diag.data$y[1,]
        #} else true.point[i,8] <- y.test[i,]
        
        dmr.point[i,1] <- mean(t)
        dmr.point[i,2] <- sd(t)
        dmr.point[i,3:7] <- quantile(t, quantile.proba)
        #dmr.point[i,8] <- dmr.point[i,1]
        
        if(!is.null(fit.gcde.extra)){
            for (j in 1:extra.gcde.length) {
                extra.dmr.point.list[[j]][i,1] <- mean(extra.t.list[[j]])
                extra.dmr.point.list[[j]][i,2] <- sd(extra.t.list[[j]])
                extra.dmr.point.list[[j]][i,3:7] <- quantile(extra.t.list[[j]], quantile.proba)
                #extra.dmr.point.list[[j]][i,8] <- extra.dmr.point.list[[j]][i,1]
            }
        }
        
        if(i%%20 == 0) invisible(gc())
    }
    
    
    ###  !!! Note that this part has not been updated to include MSPE of new Y
    #
    #Evaluate performance for non GCDE fit, parallelization is used to speed up.
    if(!is.null(otherparameter)){
        message("generate table:other methods...")
        #For each point, we evaluate the performance by parfun function.
        parfun <- function(i){
            random.x <- random.x.all[i, ]
            x.matrix <- matrix(random.x, nrow = 10000, ncol = fit.gcde$dim.x, byrow = TRUE)
            #other methods
            grid.length <- 1000
            y.grid <- seq(min(y.true.range[i,]) - max(1.5, diff(y.true.range[i,])*0.3), 
                          max(y.true.range[i,]) + max(1.5, diff(y.true.range[i,])*0.3), len = grid.length)
            grid.bandwidth <- diff(y.grid)[1]
            x.grid <- matrix(random.x, nrow = length(y.grid), ncol = fit.gcde$dim.x, byrow = TRUE)
            
            #NNKCDE
            nncde <- otherparameter$nncde.fit$predict(matrix(random.x, nrow = 1), y.grid)
            if ( abs(sum(nncde[1, ])*grid.bandwidth - 1) >0.001) warning("Numerical integration of PDF not 1.")
            nncde.mean <- sum(y.grid * nncde[1, ])*grid.bandwidth
            nncde.point.temp <- c(nncde.mean,
                                  sqrt(sum((y.grid - nncde.mean)^2 * nncde[1, ])*grid.bandwidth),
                                  sapply(quantile.proba, function(x, vec) {y.grid[which.min(abs(vec-x))]},
                                         cumsum(nncde[1, ])*grid.bandwidth) )
            #CKDE
            ckde <- fitted(npcdens(exdat=x.grid, 
                                   eydat=y.grid, 
                                   bws=otherparameter$ckde.bw))
            if ( abs(sum(ckde)*grid.bandwidth - 1) >0.001) warning("Numerical integration of PDF not 1.")
            ckde.mean <- sum(y.grid * ckde)*grid.bandwidth
            ckde.point.temp <- c(ckde.mean,
                                 sqrt(sum((y.grid - ckde.mean)^2 * ckde)*grid.bandwidth),
                                 sapply(quantile.proba, function(x, vec) {y.grid[which.min(abs(vec-x))]},
                                        cumsum(ckde)*grid.bandwidth) )
            #Flexcode
            flex.pre=predict(otherparameter$flexcode.fit, matrix(random.x, nrow = 1), B=5000)
            if ( abs( sum(flex.pre$CDE[1,])*diff(flex.pre$z)[1] - 1) >0.001) 
                warning("Numerical integration of PDF not 1.")
            flex.mean <- sum(flex.pre$z * flex.pre$CDE[1,]) * diff(flex.pre$z)[1]
            flex.point.temp <- c(flex.mean,
                                 sqrt(sum((flex.pre$z - flex.mean)^2 * flex.pre$CDE[1,]) * diff(flex.pre$z)[1] ),
                                 sapply(quantile.proba, function(x, vec) {flex.pre$z[which.min(abs(vec-x))]},
                                        cumsum(flex.pre$CDE[1,]) * diff(flex.pre$z)[1]) )
            list(nncde.temp = nncde.point.temp, ckde.temp = ckde.point.temp, flex.temp = flex.point.temp)
        }
        
        #ifclnum is NULL. DO not use parallelization.
        if(is.null(clnum)){
            all.temp <- sapply(1:s, parfun)
        }
        else{
            cl <- makeCluster(clnum)
            clusterEvalQ(cl, library(FlexCoDE))
            clusterEvalQ(cl, library(NNKCDE))
            clusterEvalQ(cl, library(np))
            clusterExport(cl, c("random.x.all", "y.true.range", "parfun", "otherparameter"),
                          envir=environment())
            all.temp <- parSapply(cl, 1:s, parfun)
            stopCluster(cl)
        }
        nncde.point <- matrix(unlist(all.temp[1,]), nrow = s, ncol = 7, byrow = TRUE)
        ckde.point <- matrix(unlist(all.temp[2,]), nrow = s, ncol = 7, byrow = TRUE)
        flex.point <- matrix(unlist(all.temp[3,]), nrow = s, ncol = 7, byrow = TRUE)
    }
    
    
    
    table.mspe <- data.frame(GCDE = colMeans((dmr.point - true.point)^2))
    if(!is.null(fit.gcde.extra)){
        for (j in 1:extra.gcde.length) {
            table.mspe <- cbind(table.mspe, colMeans((extra.dmr.point.list[[j]] - true.point)^2))
        }
        colnames(table.mspe) <- extra.names
    }
    
    if(!is.null(otherparameter)){                        
        table.mspe$NNKCDE <- colMeans((nncde.point - true.point)^2)
        table.mspe$CKDE <- colMeans((ckde.point - true.point)^2)
        table.mspe$FlexCode <- colMeans((flex.point - true.point)^2)
    }
    
    row.names(table.mspe) <- c("E(Y|X)", "SD(Y|X)", paste0("$\\tau=",as.character(quantile.proba), "$"))
                               #"Y|X")
    return(table.mspe)
}
