########################################   CDE  Graph  ##################################################
plot.cde.comparison <- function(x = NULL, n.pic, modelname, fit.gcde, fit.extra = NULL, 
                                otherparameter, true.model = NULL, 
                                ylabs2d = NULL, 
                                xlims2d = NULL, ylims2d = NULL, diag = FALSE, ...){
    if(fit.gcde$dim.y == 1){
        if(n.pic ==4){
            #x11(width = 20,height = 14)
            par(mfrow = c(2,2))}
        else {
            #x11(width = 6, height = 4)
            par(mar = c(3.6,3.6,2.5,1))
        }
        
        for(i in 1:n.pic){
            # density of DMR and truth
            if(is.null(x)) random.x <- rnorm(fit.gcde$dim.x)
            else random.x <- x
            x.matrix <- matrix(random.x, nrow = 50000, ncol = fit.gcde$dim.x, byrow = TRUE)
            if(fit.gcde$eta.dist == "normal"){
                eta <- matrix(rnorm(fit.gcde$dim.eta * 50000), nrow = 50000, ncol = fit.gcde$dim.eta)
            }
            if(fit.gcde$eta.dist == "uniform"){
                eta <- matrix(runif(fit.gcde$dim.eta * 50000,-1,1), nrow = 50000, ncol = fit.gcde$dim.eta)
            }
            t <- as.array(predict(fit.gcde$generator.t,list(x.matrix, eta)))
            
            #Compare primal and dual form.
            if(!is.null(fit.extra)){
                t.fen.extra <- as.array(fit.extra[[1]]$generator.t(list(x.matrix,eta)))
                t.don.extra <- as.array(fit.extra[[2]]$generator.t(list(x.matrix,eta)))
            }
            
            
            #generate true conditional y samples
            if(!is.null(fit.gcde$diagnostic.info)){
                diag.data <- data.generate(n = 50000, x = random.x,
                                           dim.x = fit.gcde$dim.x, dim.y = fit.gcde$dim.y, 
                                           dim.epsilon = fit.gcde$diagnostic.info$dim.epsilon,
                                           epsilon.dist = fit.gcde$diagnostic.info$epsilon.dist,
                                           model = fit.gcde$diagnostic.info$model,
                                           seed = NULL,
                                           testdata = FALSE)
                y.true <- diag.data$y
            } else y.true <- NULL
            if(modelname == "M1") t.bw <- 0.15
            else if(modelname == "M2") t.bw <- 0.45
            else if(modelname == "M3") t.bw <- 1
            else if(modelname == "M4") t.bw <- 0.3
            density.t <- density(t, bw = t.bw)
            density.y.true <- density(y.true)
            
            #set up grid
            y.grid <- seq(min(t, y.true)-0.5, max(t,y.true)+0.5, len = 1000)
            x.grid <- matrix(random.x, nrow = 1000, ncol = fit.gcde$dim.x, byrow = TRUE)
            
            #nn-CDE
            nncde <- otherparameter$nncde.fit$predict(matrix(random.x, nrow = 1), y.grid)
            #C-KDE
            ckde <- fitted(npcdens(exdat=x.grid, 
                                   eydat=y.grid, 
                                   bws=otherparameter$ckde.bw))
            #Flexcode
            flex.pre=predict(otherparameter$flexcode.fit, matrix(random.x, nrow = 1), B=500)
            
            #Plot comparison of primal and dual
            if(!is.null(fit.extra)){
                plot(density.t, lwd = 2, xlim = c(min(t, y.true), max(t,y.true)), 
                     ylim = c(0,max(density.t$y,density.y.true$y)),ann =  FALSE)
                title(xlab = "Y", ylab = "Density", cex.lab = 1.3, line = 2)
                title(main = paste("Conditional Density Estimation in Model", modelname))
                lines(density.y.true, col = "red", lwd = 2)
                lines(density(t.fen.extra), col = "blue", lwd = 2)
                lines(density(t.don.extra), col = "green", lwd = 2)
                legend("topright", lty = c(1,1,1), col = c("red", "black","blue", "green"),
                       legend = c("Truth", "GCDS(primal)", "GCDS(fenchel)","GCDS(donsker)"))
            }
            else{
                #plot comparison of primal and others.
                plot(density.t, lwd = 2, xlim = c(min(t, y.true), max(t,y.true)), 
                     ylim = c(0,max(density.t$y,density.y.true$y)),ann =  FALSE)
                title(xlab = "Y", ylab = "Density", cex.lab = 1.3, line = 2)
                title(main = paste("Conditional Density Estimation in Model", modelname))
                lines(density.y.true, col = "red", lwd = 2)
                lines(nncde[1,] ~y.grid, col = "green", lwd = 2)
                lines(ckde~y.grid, col = "blue", lwd = 2)
                lines(flex.pre$CDE[1,]~ flex.pre$z, col = "pink", lwd = 2)
                legend("topright", lty = c(1,1,1,1,1), col = c("red","black", "green", "blue","pink"), 
                       legend = c("Truth","GCDS", "NNCDE", "CKDE", "Flexcode"))
            }
        }
    }
    
    if(fit.gcde$dim.y == 2){
        x11(width = 8,height = 2*n.pic)
        layout(matrix(1:(5*n.pic), nrow = n.pic, ncol = 5, byrow = TRUE),
               widths = c(4, 1 ,3.5,3.5,3.5), heights = c(3.8,rep(3.5, n.pic-1)))
        #mgp controls the distance of axis title, axis tick label, axis to the edge of plot area.
        par(mgp = c(3, 0.8, 0))
        
        
        #generate data
        plot.size <- 10000
        for (nplot in 1:n.pic) {
            if(is.null(x)) random.x <- rnorm(fit.gcde$dim.x)
            else random.x <- x[nplot]
            x.matrix <- matrix(random.x, nrow = plot.size, ncol = fit.gcde$dim.x, byrow = TRUE)
            
            dim.eta  <- fit.gcde$dim.eta
            true.y <- data.generate(n=plot.size, x = random.x, dim.x = fit.gcde$dim.x ,dim.y = fit.gcde$dim.y
                                    , dim.epsilon = 2, epsilon.dist = "normal", model = true.model)$y
            t <- as.array(predict(fit.gcde$generator.t,
                                  list(x.matrix, 
                                                    matrix(rnorm(plot.size*dim.eta), nrow = plot.size, ncol = dim.eta))))
            
            #Grid for 2D plot
            if(!is.null(xlims2d)) y1.grid <- seq(xlims2d[[nplot]][1], xlims2d[[nplot]][2], len = 200)
            else y1.grid <- seq(min(true.y[,1]), max(true.y[,1]), len = 200)
            if(!is.null(ylims2d)) y2.grid <- seq(ylims2d[[nplot]][1], ylims2d[[nplot]][2], len = 200)
            else y2.grid <- seq(min(true.y[,2]), max(true.y[,2]), len = 200)
            y.grid <- expand.grid(y1.grid, y2.grid)
            #nnkcde
            nncde <- otherparameter$nncde.fit$predict(matrix(random.x, nrow = 1), y.grid)
            #C-KDE
            ckde <- fitted(npcdens(exdat=matrix(random.x, nrow = 40000, ncol = fit.gcde$dim.x, byrow = TRUE), 
                                   eydat=y.grid, 
                                   bws=otherparameter$ckde.bw))
            
            
            mar11 <- c(2.1, 4.1, 2.1, 0.2)
            mar12 <- c(2.1, 2.1, 2.1, 0.2)
            mar21 <- c(2.1, 4.1, 0.1, 0.2)
            mar22 <- c(2.1, 2.1, 0.1, 0.2)
            #color density legend
            library(fields)
            colorlegend <- function(maxd, nplot){
                topmargin <- ifelse(nplot == 1, 3, 2)
                savepar <- par(cex=0.7, lwd=0.5, mar = c(2.5, 0.5, topmargin, 2.4),
                               xaxs="i", yaxs="i")
                plot.new()
                plot.window(xlim=c(0, 0.4), ylim=c(0, maxd))
                rect(0, seq(0, maxd, length=65)[-65],
                     0.4, seq(0, maxd, length=65)[-1],
                     col=two.colors(n = 64, start="white", end=blues9,
                                    alpha=1.0), border=NA)
                box()
                axis(4, at=seq(0,maxd, by = 0.1), las=1)
                par(savepar)
            }
            getd <- function(){ 
                z <- get('dens', envir = parent.frame(1))
                dd <<- z}
            if(is.null(ylabs2d)) ylabs2d <- rep("", nplot)
            
            #Plot
            if(nplot == 1) {par(mar = mar11)
                titles <- c("True conditional samples", "GCDE", "NNKCDE", "CKDE")}
            else {par(mar = mar21)
                titles <- NULL}
            dd <- 0
            smoothScatter(true.y, nrpoints = 0, nbin = 200, bandwidth = 0.15,
                          main = titles[1], xlab ="", ylab = ylabs2d[nplot], cex.lab = 1.2,
                          postPlotHook = getd, xlim = xlims2d[[nplot]], ylim = ylims2d[[nplot]], ...)
            #add color density legend
            colorlegend(max(dd), nplot)
            
            if(nplot == 1) par(mar = mar12)
            else par(mar = mar22)
            smoothScatter(t, nrpoints = 0, nbin = 200, bandwidth = 0.15,
                          main = titles[2], xlab = "Y1", ylab = "", 
                          xlim = xlims2d[[nplot]], ylim = ylims2d[[nplot]],...)
            image(x = y1.grid, y = y2.grid, z = matrix(nncde, nrow = 200, ncol = 200), 
                  col = hcl.colors(30,"blues", rev = TRUE),
                  main = titles[3], xlab = "Y1", ylab = "", xaxs = "r", yaxs = "r",
                  xlim = xlims2d[[nplot]], ylim = ylims2d[[nplot]], ...)
            image(x = y1.grid, y = y2.grid, z = matrix(ckde, nrow = 200, ncol = 200), 
                  col = hcl.colors(30,"blues", rev = TRUE),
                  main = titles[4], xlab = "Y1", ylab = "", xaxs = "r", yaxs = "r",
                  xlim = xlims2d[[nplot]], ylim = ylims2d[[nplot]], ...)
        }
    }
    
    if(diag == TRUE) list(t = t, y.true = y.true)
    else return(invisible())
}