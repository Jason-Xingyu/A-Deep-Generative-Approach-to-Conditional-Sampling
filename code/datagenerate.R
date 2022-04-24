###################################   data generating function  ##########################################
#Function   data.generate
#Description of arguments:
#n:           Size of the dataset.
#x = NULL:    If not null, use given x to generator model.  If null, generate x from normal.
#dim.x:       Dimension of X.
#dim.y:       Dimension of Y.
#dim.epsilon: Dimension of epsilon. (error term)
#epsilon.dist:The distribution of epsilon.
#model:       Specify which model to use.
#seed = NULL: Seed for random number generator.
#testdata = FALSE:  If TRUE, will also generate a test dataset with the same setting as before.

data.generate <- function(n, x = NULL, x.dist = "normal",
                          dim.x, dim.y, dim.epsilon,
                          epsilon.dist,
                          model,
                          seed = NULL,
                          testdata = FALSE,
                          extrapara = NULL){
    if(!is.null(seed)) set.seed(seed)
    generate <- function(){
        #Generate epsilon
        if(epsilon.dist == "normal"){
            epsilon <- matrix(rnorm(n*dim.epsilon, mean = 0, sd = 1), nrow = n, ncol = dim.epsilon)
        }
        if(epsilon.dist == "uniform"){
            epsilon <- matrix(runif(n*dim.epsilon, -1,1), nrow = n, ncol = dim.epsilon)
        }
        if(epsilon.dist == "mixture"){
            epsilon <- matrix(ifelse(runif(n*dim.epsilon)<0.5, rnorm(n*dim.epsilon, mean = 2, sd = 1), 
                                     rnorm(n*dim.epsilon, mean = -2, sd = 1)), nrow = n, ncol = dim.epsilon)
        }
        
        #Generate x: IF x is nonnull, then just repeat x n times.
        if(!is.null(x)) {
            if(length(x) != dim.x) warning("The length of input x does not match dim.x. 
                                     There might be mistake If you are not prividing the whole x matrix.")
            x <- matrix(x, nrow = n, ncol = dim.x, byrow = TRUE)
            xgiven <- TRUE
            
        } else {
            xgiven <- FALSE
            if(x.dist == "space.corr"){
                message("Cov of X is using spacial correlation.")
                x.cov <- matrix(0, dim.x,dim.x)
                for (q in 1:dim.x) {
                    x.cov[q,] <- 0.5^abs(q - 1:dim.x)
                }
                x <- mvrnorm(n, mu = rep(0,dim.x), Sigma = x.cov)
            }else if(x.dist == "mixture"){
                x <- matrix(ifelse(runif(n*dim.x)>0.5, rnorm(n*dim.x, -2,1), rnorm(n*dim.x, 2, 1)), 
                                    nrow = n, ncol = dim.x)
            }else if(x.dist == "normal"){
                x <- matrix(rnorm(n*dim.x), nrow = n, ncol = dim.x)
            }else if(x.dist == "uniform"){
                x <- matrix(runif(n*dim.x, 0,1), nrow = n, ncol = dim.x)
            }else if(x.dist == "uniform2"){
                x <- matrix(runif(n*dim.x, -2,2), nrow = n, ncol = dim.x)
            }else stop("The specified distribution of X is not supported yet.")
        }
        
        #Generate y 
        y <- matrix(0, nrow = n, ncol = dim.y)
        y.mean <- y
        if(dim.y != 1) warning("The dimension of y is not 1. Cannot calculate SD.")
        y.sd <- rep(0,n)
        
        if(model == "add-toy")
            for(i in 1:n){
                u <- runif(1)
                if(u<1/2) y.mean[i, ] <- x[i,1]
                if(u>1/2) y.mean[i, ] <- -x[i,1]
                y[i,] <- y.mean[i,] + 0.25 * epsilon[i,1]
                y.mean[i, ] <- 0
                y.sd[i] <- sqrt(0.25^2 + x[i,1]^2)}
        
        else if (model == "add-mix-wgangp")
            for(i in 1:n){
                u <- runif(1)
                if(u<1/2) y.mean[i, ] <- x[i,1]
                if(u>1/2) y.mean[i, ] <- -x[i,1]
                y[i,] <- y.mean[i,] + 0.5 * epsilon[i,1]
                y.mean[i, ] <- 0
                y.sd[i] <- sqrt(0.5^2 + x[i,1]^2)}
        
        else if (model == "add-mix-unsymmetric")
            for(i in 1:n){
                u <- runif(1)
                if(u<2/3) y[i, ] <-  1 + x[i,1] + 0.5*x[i,2] + 1 * epsilon[i,1] 
                if(u>2/3) y[i, ] <- -1 - x[i,1] - 0.5*x[i,2] + 0.5 * epsilon[i,1] 
                #For mean and sd, see wiki:mixture distribution.
                y.mean[i, ] <- 1/3 * (1 + x[i,1] + 0.5*x[i,2])
                y.sd[i] <- sqrt(8/9 * (1 + x[i,1] + 0.5*x[i,2])^2 + 1/3/4 + 2/3)
                }
        
        
        else if(model == "add1") 
            for(i in 1:n){
                y.mean[i,] <- x[i, 1]^2 + x[i,2] + x[i,3] + x[i,4] + x[i,5]
                y[i,] <- y.mean[i,] + 1 * epsilon[i,1]
                y.sd[i] <- 1}
        
        else if(model == "add2")
            for(i in 1:n){
                y.mean[i,] <- x[i, 1]^2 + exp((x[i,2] + x[i,3])/3) + sin(x[i,4]*x[i,5])
                y[i,] <- y.mean[i,] + 1 * epsilon[i,1]
                y.sd[i] <- 1}
        
        else if(model == "add-hete")
            for(i in 1:n){
                y.mean[i,] <- x[i,1]^2 + exp((x[i,2] +x[i,3]/3)) + x[i,4] - x[i,5]
                y[i,] <- y.mean[i,] + (0.5+ x[i,2]^2/2 +x[i,5]^2/2)*epsilon[i,1]
                y.sd[i] <- (0.5+ x[i,2]^2/2 +x[i,5]^2/2)}
        
        else if(model == "mul-simple")
            for(i in 1:n){
                temp <- (x[i, 1]^2/3 + x[i,2]/2 + x[i,3]/2 + x[i,4] + x[i,5])
                y.mean[i,] <- temp * exp(1/8)
                y[i,] <- temp * exp(0.5 * epsilon[i,1])
                y.sd[i] <- temp * sqrt((exp(0.25)-1) * exp(0.25))}
        
        
        else if(model == "mul-hete")
            for(i in 1:n){
                y[i,] <- (4 + x[i,1]+x[i,2]^2+x[i,3]) * exp((0.1+abs(x[i,2]/3 + x[i,3])/3) * epsilon[i,1])
                y.mean[i,] <- (4 + x[i,1]+x[i,2]^2+x[i,3]) * exp(((0.1+abs(x[i,2]/3 + x[i,3])/3)^2) /2) 
                y.sd[i] <- (4 + x[i,1]+x[i,2]^2+x[i,3]) * sqrt((exp((0.1+abs(x[i,2]/3 + x[i,3])/3)^2) - 1) * 
                                                                   exp(((0.1+abs(x[i,2]/3 + x[i,3])/3))^2))}
        
        else if(model == "mul-simple-moved")
            for(i in 1:n){
                temp <- (5 + x[i, 1]^2/3 + x[i,2]/2 + x[i,3]/2 + x[i,4] + x[i,5])
                y[i,] <- temp * exp(0.5 * epsilon[i,1])
                if(epsilon.dist == "normal"){
                    y.mean[i,] <- temp * exp(1/8)
                    y.sd[i] <- temp * sqrt((exp(0.25)-1) * exp(0.25))}
                if(epsilon.dist == "mixture"){
                    mu1 <- temp * exp(1+1/8)
                    mu2 <- temp * exp(-1 + 1/8)
                    sigma1.square <- temp^2 * (exp(0.25)-1) * exp(2 + 0.25)
                    sigma2.square <- temp^2 * (exp(0.25)-1) * exp(-2 + 0.25)
                    y.mean[i, ] <- (mu1 + mu2)/2
                    y.sd[i] <-  sqrt(0.5 * (mu1^2 + mu2^2 + sigma1.square + sigma2.square - 2 * y.mean[i,]^2))  }}
        
        else if(model == "mul-simple-moved-x1new")
            for(i in 1:n){
                temp <- (5 + x[i, 1]^2 + x[i,2]/2 + x[i,3]/2 + x[i,4] + x[i,5])
                y[i,] <- temp * exp(0.5 * epsilon[i,1])
                if(epsilon.dist == "normal"){
                    y.mean[i,] <- temp * exp(1/8)
                    y.sd[i] <- temp * sqrt((exp(0.25)-1) * exp(0.25))}
                if(epsilon.dist == "mixture"){
                    mu1 <- temp * exp(1+1/8)
                    mu2 <- temp * exp(-1 + 1/8)
                    sigma1.square <- temp^2 * (exp(0.25)-1) * exp(2 + 0.25)
                    sigma2.square <- temp^2 * (exp(0.25)-1) * exp(-2 + 0.25)
                    y.mean[i, ] <- (mu1 + mu2)/2
                    y.sd[i] <-  sqrt(0.5 * (mu1^2 + mu2^2 + sigma1.square + sigma2.square - 2 * y.mean[i,]^2))  }}
        
        
        else if(model == "mul-simple-moved-log")
            for(i in 1:n){
                temp <- (5 + x[i, 1]^2/3 + x[i,2]/2 + x[i,3]/2 + x[i,4] + x[i,5])
                y[i,] <- log(temp) + 0.5 * epsilon[i,1]
                if(epsilon.dist == "normal"){
                    y.mean[i,] <- log(temp) * exp(1/8)
                    y.sd[i] <- temp * sqrt((exp(0.25)-1) * exp(0.25))}
                if(epsilon.dist == "mixture"){
                    mu1 <- temp * exp(1+1/8)
                    mu2 <- temp * exp(-1 + 1/8)
                    sigma1.square <- temp^2 * (exp(0.25)-1) * exp(2 + 0.25)
                    sigma2.square <- temp^2 * (exp(0.25)-1) * exp(-2 + 0.25)
                    y.mean[i, ] <- (mu1 + mu2)/2
                    y.sd[i] <-  sqrt(0.5 * (mu1^2 + mu2^2 + sigma1.square + sigma2.square - 2 * y.mean[i,]^2))  }}
        
        
        else if(model == "mixture")
            for(i in 1:n){
                u <- runif(1)
                if(u<1/3) y.mean[i, ] <- 2 * x[i,1]+x[i,2]^2
                if(1/3<u & u<2/3) y.mean[i, ] <- x[i,2] - 2*sin(x[i,3])
                if(u>2/3) y.mean[i, ] <- x[i,4] + x[i,5] + 2*x[i,1]^2
                y[i,] <- y.mean[i, ] + 1 * epsilon[i,1]
                y.sd[i] <- sqrt(1/3 * (3 + (2 * x[i,1]+x[i,2]^2)^2 + (x[i,2] - 2*sin(x[i,3]))^2 + 
                                           (x[i,4] + x[i,5] + 2*x[i,1]^2)^2 
                                       -3*((2 * x[i,1]+x[i,2]^2 +x[i,2] - 2*sin(x[i,3]) + x[i,4] + x[i,5] + 2*x[i,1]^2)/3)^2
                ) )}
        
        else if(model == "2dmixture"){
            for (i in 1:n) {
                u <- runif(1)
                if(u<1/4) y.mean[i, ] <- c(x[i,1], x[i,1])
                if(1/4<u & u<1/2) y.mean[i, ] <- c(-x[i,1], x[i,1])
                if(1/2<u & u<3/4) y.mean[i, ] <- c(-x[i,1], -x[i,1])
                if(3/4<u) y.mean[i, ] <- c(x[i,1], -x[i,1])
                y[i,] <- y.mean[i, ] + 0.15 * epsilon[i,]
                y.mean[i,] <- c(0,0)
            }
            y.sd <- NULL}
        
        else if(model == "2dmixtureofring"){
            for (i in 1:n) {
                u <- runif(1)
                y.mean[i, ] <- c(x[i,1]*cos(2*pi*u), x[i,1]*sin(2*pi*u))
                y[i,] <- y.mean[i, ] + 0.05 * epsilon[i,]
                y.mean[i,] <- c(0,0)
            }
            y.sd <- NULL}
        
        else if(model == "2dhyperbolic"){
            for (i in 1:n) {
                u <- runif(1, min = 0, max = 2*pi)
                y.mean[i, ] <- c(x[i,1]/u*cos(u), x[i,1]/u*sin(u))
                y[i,] <- y.mean[i, ] + 0.03 * epsilon[i,]
                y.mean[i,] <- c(0,0)
            }
            y.sd <- NULL}
        
        else if(model == "2dinvolute"){
            for (i in 1:n) {
                y.mean[i, ] <- c(x[i,1] * sin(2*x[i,1]), x[i,1] * cos(2*x[i,1]))
                y[i,] <- y.mean[i, ] + 0.9 * epsilon[i,]
            }
            y.sd <- NULL}
        
        else if(model == "2dcondi.involute"){
            for(i in 1:n) {
                u <- runif(1, 0, 2*pi)
                y[i, ] <- c(2*x[i,1] + u*sin(2*u), 2*x[i,1] + u*cos(2*u)) + extrapara * epsilon[i, ]
            }
            y.mean <- NULL
            y.sd <- NULL
        }
        
        else if(model == "ddr1"){
            core <- cbind(rowSums(x[,1:3]), rowSums(x[,3:6]), rowSums(x[,7:9]))
            if(dim.y !=1) stop("ddr model support only dim.y= 1 currently")
            y.mean <- 2*core[,1] + core[,2]^2 + 10 * sin(core[,3])
            y <- matrix(y.mean + extrapara * as.vector(epsilon), nrow = n, ncol = 1)
            y.mean <- matrix(y.mean, nrow = n, ncol = 1)
            y.sd <- rep(1, n)
        }
        
        else if(model == "ddr2"){
            #model A in KDR
            if(dim.y !=1) stop("ddr model support only dim.y= 1 currently")
            y.mean <- x[,1]/(0.5 + (x[,2]+1.5)^2) + (1 + x[,2])^2
            y <- matrix(y.mean + extrapara * as.vector(epsilon), nrow = n, ncol = 1)
            y.mean <- matrix(y.mean, nrow = n, ncol = 1)
            y.sd <- rep(1, n)
        }
        
        else if(model == "ddr3"){
            #model2 in ma and zhu 2013
            if(dim.y !=1) stop("ddr model support only dim.y= 1 currently")
            #use the 1,2 columns of x as x1, x2.   5, 6 columns of x as e1, e2
            if(xgiven) {
                #If X given, we calculate one row from true model and repeat n times.
                xx <- x[1,]
                xx[3] <- 0.2 * xx[1] + 0.2*(xx[2]+2)^2 + 0.2*xx[5]
                xx[4] <- 0.1 + 0.1*(xx[1] + xx[2]) + 0.3 * (xx[1] + 1.5)^2 + 0.2 *xx[6]
                xx[5] <- rbinom(1, 1, exp(xx[1])/(1+exp(xx[1])) )
                xx[6] <- rbinom(1, 1, exp(xx[2])/(1+exp(xx[2])) )
                x <- matrix(xx, nrow = n, ncol = dim.x, byrow = TRUE)
            } else {
            x[,3] <- 0.2 * x[,1] + 0.2*(x[,2]+2)^2 + 0.2*x[,5]
            x[,4] <- 0.1 + 0.1*(x[,1] + x[,2]) + 0.3 * (x[,1] + 1.5)^2 + 0.2 *x[,6]
            x[,5] <- rbinom(n, 1, exp(x[,1])/(1+exp(x[,1])) )
            x[,6] <- rbinom(n, 1, exp(x[,2])/(1+exp(x[,2])) )
            }
            beta <- c(1.3, -1.3, 1.0, -0.5, 0.5, -0.5)
            xbeta <- x[,1:6] %*% beta
            
            y.mean <- as.vector(sin(2 * xbeta) + 2*exp(2 + xbeta)) 
            y <- matrix(y.mean + as.vector(sqrt(log(2+xbeta^2))) * as.vector(epsilon),
                        nrow = n, ncol = 1)
            y.mean <- matrix(y.mean, nrow = n, ncol = 1)
            y.sd <- as.vector(sqrt(log(2+xbeta^2)))
        }
        
        else if(model == "ddr4"){
            #for high dimensional model testing
            if(dim.y !=1) stop("ddr model support only dim.y= 1 currently")
            y.mean <- (x[,1]+x[,2]+x[,3])/(0.5 + (x[,4]+1.5)^2) + (1 + x[,5])^2 +
                x[,6]^2 + (0.5 * x[,7] + 0.5 * x[,8])^2 + 10*sin(x[,9] + x[,10]+x[,11])
        
            y <- matrix(y.mean + extrapara * as.vector(epsilon), nrow = n, ncol = 1)
            y.mean <- matrix(y.mean, nrow = n, ncol = 1)
            y.sd <- rep(1, n)
        }
        
        else if(model == "ddr5"){
            #Efficient Sparse Estimate of Sufficient Dimension Reduction in High Dimension
            #[Xin Chen, Wenhui Sheng, Xiangrong Yin, 2017 Technometrics]
            if(dim.y !=1) stop("ddr model support only dim.y= 1 currently")
            y.mean <- 0.5*(x[,1] + x[,2]+x[,3]+x[,4])^2 + 3*sin(0.25*(x[,1]-x[,2]+x[,3]-x[,4]))
            y <- matrix(y.mean + 0.2*as.vector(epsilon), nrow = n, ncol = 1)
            y.mean <- matrix(y.mean, nrow = n, ncol = 1)
            y.sd <- rep(0.2, n)
        }
        
        else if(model == "ddr6"){
            #[Qian Lin, Zhigen Zhao, Jun S. Liu, 2018 JASA]
            #Sparse Sliced Inverse Regression Via Lasso
            #model VI
            if(dim.y !=1) stop("ddr model support only dim.y= 1 currently")
            z1 <- x[,1]+x[,2] + x[,3]+x[,4]
            z2 <- x[,5] + x[,6]+x[,7]
            y.mean <- abs((z2/4 + 2)^3)* z1/abs(z1)
            y <- matrix(y.mean + 1*as.vector(epsilon), nrow = n, ncol = 1)
            y.mean <- matrix(y.mean, nrow = n, ncol = 1)
            y.sd <- rep(1, n)
        }
        else if(model == "ddr7"){
            #[Qian Lin, Zhigen Zhao, Jun S. Liu, 2018 JASA]
            #Sparse Sliced Inverse Regression Via Lasso
            #Model VII
            if(dim.y !=1) stop("ddr model support only dim.y= 1 currently")
            z1 <- x[,1]+x[,2] + x[,3]+x[,4] + x[,5] + x[,6] + x[,7]
            z2 <- x[,8] + x[,9]+x[,10] + x[,11] + x[,12]
            y.mean <- z1 * exp(z2/8)
            y <- matrix(y.mean + 1*as.vector(epsilon), nrow = n, ncol = 1)
            y.mean <- matrix(y.mean, nrow = n, ncol = 1)
            y.sd <- rep(1, n)
        }
        else if(model == "ddr8"){
            #[Qian Lin, Zhigen Zhao, Jun S. Liu, 2018 JASA]
            #Sparse Sliced Inverse Regression Via Lasso
            #Model IX
            if(dim.y !=1) stop("ddr model support only dim.y= 1 currently")
            z1 <- x[,1]+x[,2] + x[,3]+x[,4] + x[,5] + x[,6] + x[,7] + x[,8]
            z2 <- x[,9]+x[,10] + x[,11] + x[,12]
            y.mean <- z1 * (2 + z2/3)^2
            y <- matrix(y.mean + 1*as.vector(epsilon), nrow = n, ncol = 1)
            y.mean <- matrix(y.mean, nrow = n, ncol = 1)
            y.sd <- rep(1, n)
        }
        
        else stop("The model specified is not suppored.")
        
        #Error check
        if((model != "mul-simple-moved") & (epsilon.dist != "normal")
           & (model != "mul-simple-moved-x1new")) 
            stop("the sd calculation is not correct when the 
                                    distribution of epsilon is not normal")
        return(list(x = x, y = y, y.mean = y.mean, y.sd = y.sd, epsilon = epsilon))
    }
    
    data1 <- generate()
    
    
    ## Generate test data if needed
    
    if(testdata == TRUE){
        datatest <- generate()
    } else {datatest <- list(x = NULL,
                             y = NULL,
                             y.mean = NULL,
                             y.sd = NULL,
                             epsilon = NULL)}
    
    
    diagnostic.info <- list(dim.epsilon = dim.epsilon, epsilon.dist = epsilon.dist, model = model,
                            test.y.mean = datatest$y.mean, extrapara = extrapara,
                            x.dist = x.dist)
    
    return(list(x = data1$x, y = data1$y, y.mean = data1$y.mean, y.sd = data1$y.sd, epsilon = data1$epsilon,
                test.x = datatest$x, test.y = datatest$y, test.y.mean = datatest$y.mean,
                test.y.sd = datatest$y.sd, test.epsilon = datatest$epsilon,
                n = n, dim.x = dim.x, dim.y = dim.y, dim.epsilon = dim.epsilon, model = model, seed = seed,
                diagnostic.info = diagnostic.info, extrapara = extrapara))
}
#**********************************   END OF DATA Generating function   ***********************************