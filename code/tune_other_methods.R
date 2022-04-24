########################################   Tune other methods  ################################################
#library(kedd)

simulation.tune <- function(fit.gcde, nnkcde.h.grid, nnkcde.k.grid, nimax){
    set.seed(153)
    plot.train.index <- 1:round(0.7*nrow(fit.gcde$train.x))
    #knn-cde
    nncde.fit <- NNKCDE$new(fit.gcde$train.x[plot.train.index, ], fit.gcde$train.y[plot.train.index, ])
    if(fit.gcde$dim.y == 1)
        nncde.fit$tune(fit.gcde$train.x[-plot.train.index, ], fit.gcde$train.y[-plot.train.index, ],
                       k_grid = nnkcde.k.grid, h_grid = nnkcde.h.grid)
    if(fit.gcde$dim.y == 2){
        nncde.fit$tune(fit.gcde$train.x[-plot.train.index, ], fit.gcde$train.y[-plot.train.index, ],
                       k_grid = nnkcde.k.grid)
        nncde.fit$h <- ks::Hpi(fit.gcde$train.y[plot.train.index, ])
        #nncde.fit$h <- h.amise.default(fit.gcde$train.y[plot.train.index, ])$h
    }
    #C-KDE
    bw <- npcdensbw(xdat = fit.gcde$train.x, ydat = fit.gcde$train.y, bwmethod = "normal-reference")
    #Flexcode
    if(fit.gcde$dim.y == 1)
        flexcode.fit <- fitFlexCoDE(fit.gcde$train.x[plot.train.index, ], fit.gcde$train.y[plot.train.index, ],
                                    fit.gcde$train.x[-plot.train.index, ],
                                    fit.gcde$train.y[-plot.train.index, ],
                                    nIMax = nimax, regressionFunction = regressionFunction.NN)
    else flexcode.fit <- NULL
    
    return(list(nncde.fit = nncde.fit, ckde.bw = bw, flexcode.fit = flexcode.fit))
}
