########################
#This is a copy of  fitgcde.R
#This file contains moditiocation of the network structure to fit the CGAN MNIST example.
########################
#1/4 of image is used as X .


##################################   SET UP  ##########################################
library(keras)
library(tensorflow)
library(xtable)
library(progress)#progress bar
library(parallel)#Used in generate table
library(pracma)#gramschimidt for orthonormal matrix
#For simulation comparison
require(FlexCoDE)
require(NNKCDE)
require(np)

#Set to use double-precision floating point
is_keras_available()
tf$keras$backend$set_floatx('float32')
#set to only allocate GPU memory when needed.
if(length(tf$config$experimental$list_physical_devices("GPU")))
    tf$config$experimental$set_memory_growth(tf$config$experimental$list_physical_devices("GPU")[[1]], TRUE)

#######################################################################################################
############################     Main function for fitting GCDE  ######################################
#######################################################################################################
#Description of the arguments:
## train.x, train.y, test.x, test.y:  training data and testing data
#load.gene = NULL, load.disc = NULL:  If not null, directly load weights from the file path.
#diagnostic.info:    If null, there is no in-training diagnostic.
#                    This object is created by function 'data.gene'.
#v.label:       "new" or "old.  If "old, use original formulation by Yuling, 1 means true data.
# If "new, use new formulation.  -1 means true data.
#                   Note, this will only affect the primal form.
#kl.dual:         TRUE or FALSE.   If TURE, use dual form of KL divergence.
#dual.form:      "fenchel" OR "donsker". Fenchel is the default.
#iterations: Number of max iterations of updating generator.
#dim.eta;    The dimension of eta which is used in the generator.
#eta.dist:   The distribution of eta.
#useddr:     TRUE or FALSE. If true, ddr is used to reduced dimension of x before gene and disc.
#lambda.ortho.loss:     If 0,  no orthonormal constraints loss added. If not 0, then this is the 
#                        lambda value used to orthonormal constraints.
#dr.nodes.num:            The nodes number of dimension reduction. Currently only one number.
#dr.regularizer：         regularizer for dr network.
#gene.nodes.num:   vector.   The number of neurons in each layer of generator.
#dis.nodes.num:    vector.   The number of neurons in each layer of discriminator
#disc.train.type: "batch" or "fit".   Train discriminator by batch or by fit.
#                  "fit" is used when you want to update discriminator more often. This is faster.
#step.disc:     If dis.train.type is "batch", this term determines how many times it is updated per one
#               update iteration of generator.
#step.gene;     This terms gives how many times generator is updated per one iteration.
#display:       per how many iterations, you want to display the diagnostic information.

#batch_size:    batch size in training.
#gene.learning.rate:   Learning rate of generator.
#disc.learning.rate:   Learning rate of discriminator.
#gene.lr.decay:        Decay of generator learning rate. If 0, no decay.  lr = lr * 1/(1+ decay * iterations)
#disc.lr.decay:        Decay of discriminator learning rate. If 0, no decay.
#gene.regularizer.1:   Regularizer of first layer(input layer) of generator.
#                      regularizer_l1(l = 0.1)  or regularizer_l2(l = 0.001)
#gene.regularizer.o:   Regularizer of the rest layers(other than first) of generator.
#disc.regularizer.1:   Regularizer of the first layer.
#disc.regularizer:     Regularizer of discriminator.
#activation.disc:      Activation function of discriminator. Default "relu".
#activation.gene:      Activation function of generator. Default "relu".
#activation.gene.out:  Activation funcion for the output layer of generator. Default is "NULL",i.e. linear.
#dropout.rate:         Rate for dropout layer. Default is 0, i.e. no dropout.
#lambda.kl:            The lambda weight of KL divergence in loss funtion.
#lambda.mse:           The lambda weight of MSE in loss function.
#if this term is 0, no mse penalty term is used.
#m:                    The number of observations used to calculate conditional mean( in generator.g).
#extra.function:       Extra function ran in diagnostic step. Only for simulation. Default NULL.
#vanillaNN:            If true, the vanilla NN will be trained. Then the diagnostic table can direct tuning.
fitgcde.mnist.xhalfimage.cnn <- function(train.x, train.y, test.x, test.y,
                                     load.gene = NULL,
                                     load.disc = NULL,
                                     diagnostic.info,
                                     show.diagnostic.pic = TRUE,
                                     v.label,
                                     kl.dual,
                                     dual.form = "fenchel",
                                     iterations,
                                     dim.eta,
                                     eta.dist,
                                     useddr = FALSE,
                                     lambda.ortho.loss = 0,
                                     dr.nodes.num = NULL,
                                     dr.regularizer = NULL,
                                     gene.nodes.num,
                                     dis.nodes.num,
                                     disc.train.type,
                                     step.disc,
                                     step.gene,
                                     display,
                                     batch_size,
                                     gene.learning.rate,
                                     disc.learning.rate,
                                     gene.lr.decay,
                                     disc.lr.decay,
                                     gene.regularizer.1,
                                     gene.regularizer.o,
                                     disc.regularizer.1 = NULL,
                                     disc.regularizer,
                                     activation.disc,
                                     activation.gene,
                                     activation.gene.out = NULL,
                                     dropout.rate = 0,
                                     lambda.kl,
                                     lambda.mse,         
                                     m,
                                     extra.function = NULL,
                                     vanillaNN = TRUE,
                                     savepath = NULL
)
{
    #set up basic information:
    #if(class(train.x) != "matrix" | class(train.y) != "matrix" | 
    #   class(test.x) != "matrix"  | class(test.y) != "matrix" ) stop("The dataset are not matrix!")
    dim.x <- NULL
    dim.y <- NULL
    trainningsetsize <- nrow(train.x)
    
    
    
    
    ###################################   Define discriminator   ##########################################
    #Whether to freeze the weight of reducer when updating discriminator
    #freeze_weights(reducer)
    
    discriminator_input_x <- layer_input(shape = c(28, 14,1))
    discriminator_input_y <- layer_input(shape = c(28, 28, 1))
    
    disc.x.cnn <- layer_conv_2d(discriminator_input_x, 32, kernel_size = 5, strides = 2,
                                padding = "same")%>%
        layer_batch_normalization(trainable = TRUE)%>%
        layer_activation_leaky_relu(alpha = 0.2)%>%
        layer_conv_2d(64, kernel_size = 4, strides = 1, padding = "same")%>%
        layer_batch_normalization(trainable = TRUE)%>%
        layer_activation_leaky_relu(alpha = 0.2)%>%
        layer_flatten()
    disc.y.cnn <- layer_conv_2d(discriminator_input_y, 64, kernel_size = 5, strides = 2,
                                padding = "same")%>%
        layer_batch_normalization(trainable = TRUE)%>%
        layer_activation_leaky_relu(alpha = 0.2)%>%
        layer_conv_2d(128, kernel_size = 5, strides = 2, padding = "same")%>%
        layer_batch_normalization(trainable = TRUE)%>%
        layer_activation_leaky_relu(alpha = 0.2)%>%
        layer_flatten()
    disc.xy <- layer_concatenate(list(disc.x.cnn, disc.y.cnn), axis = -1L)%>%
        #layer_conv_2d(128, kernel_size = 7, strides = 1, padding = "valid")%>%
        #layer_batch_normalization(trainable = TRUE)%>%
        #layer_activation_leaky_relu(alpha = 0.2)%>%
        #layer_flatten()%>%
        # 64 -64
        layer_dense(256)%>%
        layer_batch_normalization(trainable = TRUE)%>%
        layer_activation_leaky_relu(alpha = 0.2)%>%
        layer_dense(128)%>%
        layer_batch_normalization(trainable = TRUE)%>%
        layer_activation_leaky_relu(alpha = 0.2)
    
    disc.output <- layer_dense(disc.xy, 1, activation = NULL)
    
    
    
    discriminator <- keras_model(inputs = c(discriminator_input_x, discriminator_input_y), 
                                 outputs = disc.output)
    
    discriminator.structure <- capture.output(summary(discriminator))
    #summary(discriminator)
    
    # To stabilize training, we use learning rate decay
    # and gradient clipping (by value) in the optimizer.
    #discriminator.optimizer <- optimizer_adam( 
    #    lr = disc.learning.rate, 
    #    decay = disc.lr.decay)
    discriminator.optimizer <- tf$keras$optimizers$Adam(learning_rate = disc.learning.rate, beta_1 = 0.5)
    
    
    #Note that y_pred is actually V. V is the label of {-1,1} indicating True or generated Y.
    if(kl.dual == FALSE)
        discri.loss <- function(y_true, y_pred){
            #Add orthonormal constraints here
            if(useddr & (lambda.ortho.loss != 0)){
                weightmatrix <- reducer$weights[[1]]
                btb <- tf$matmul(tf$transpose(weightmatrix), weightmatrix)
                loss.orthonormal <- 
                    tf$norm(btb - tf$linalg$diag(k_constant(1, dtype = "float32", shape = dr.nodes.num)))
            } else loss.orthonormal <- tf$constant(0, dtype = "float32")
            
            k_mean(tf$math$log(k_constant(1, dtype = "float32") + tf$math$exp(-y_true*y_pred) ))+ 
                k_constant(lambda.ortho.loss , dtype = "float32") * loss.orthonormal
        }
    
    #Note that y_true is v labels. V is 1 if it is from true data. V is -1 if it is from generated(fake).
    if(kl.dual == TRUE){
        if(v.label == "old")
            stop("Use new labels in dual form, please.")
        
        if(v.label == "new")
            if(dual.form == "fenchel"){
                message("Fenchel dual used.")
                discri.loss <- function(y_true, y_pred){
                    #Fenchel dual
                    - k_mean(y_pred * (y_true + 1) * (1/2)  - tf$math$exp(y_pred) * (y_true - 1) * (-1/2) )
                }
            } else if(dual.form == "donsker"){
                message("Donsker dual used.")
                discri.loss <- function(y_true, y_pred){
                    #Donsker dual
                    - k_mean(y_pred * (y_true + 1) * (1/2) * 2) + 
                        tf$math$log(k_mean(tf$math$exp(y_pred) * (y_true - 1) * (-1/2)) * 2  )
                }
            } else stop("Either fenchel or donsker dual may be used in dual form.")
    }
    
    
    discriminator %>% compile(
        optimizer = discriminator.optimizer,
        loss = discri.loss
    )
    #*********************************   END of Define discriminator  ********************************
    
    
    
    ###################################   Define GENERATOR   ##########################################
    # Set discriminator weights to non-trainable
    # (will only apply to the `generator` model)
    freeze_weights(discriminator)
    #freeze_weights(discriminator)
    
    
    gene.input.x              <- layer_input(shape = c(28, 14, 1))
    gene.input.eta            <- layer_input(shape = c(dim.eta))
    
    #gene.conv <- layer_conv_2d(gene.input.x, 32, kernel_size = 5, strides = 2,
    #                           padding = "same")%>%
    #    layer_batch_normalization(trainable = TRUE)%>%
    #    layer_activation_leaky_relu(alpha = 0.2)%>%
    #    layer_conv_2d(64, kernel_size = 4, strides = 2, padding = "valid")%>%
    #    layer_batch_normalization(trainable = TRUE)%>%
    #    layer_activation_leaky_relu(alpha = 0.2)%>%
    #    layer_flatten()
    gene.conv <- layer_flatten(gene.input.x)%>%
        #64
        layer_dense(128)%>%
        layer_batch_normalization(trainable = TRUE)%>%
        layer_activation_leaky_relu(alpha = 0.2)
    gene.combined <- layer_concatenate(list(gene.conv, gene.input.eta))%>%
        layer_reshape(target_shape = c(1,1,228))
    gene.output.t <- layer_conv_2d_transpose(gene.combined, 256, kernel_size = c(7,7), strides = c(1,1), 
                                             padding = "valid")%>%
        layer_batch_normalization(trainable = TRUE)%>%
        layer_activation_leaky_relu(alpha = 0.2)%>%
        layer_conv_2d_transpose(128, kernel_size = c(5,5), strides = c(2,2), padding = "same")%>%
        layer_batch_normalization(trainable = TRUE)%>%
        layer_activation_leaky_relu(alpha = 0.2)%>%
        layer_conv_2d_transpose(1, kernel_size = c(5,5), strides = c(2,2), padding = "same")%>%
        layer_activation(activation = "tanh")
    
    
    
    #calculate the log density ratio
    gene.output.d <- discriminator(list(gene.input.x, gene.output.t))
    
    #estimate the conditional expectation
    generator.t <- keras_model(inputs = c(gene.input.x, gene.input.eta), outputs = c(gene.output.t))
    
    
    #t is T(X,eta). d is D(T(X,eta)). g is g(X).
    
    generator <- keras_model(inputs = c(gene.input.x, gene.input.eta), 
                             outputs = c(gene.output.t, gene.output.d))
    
    generator.structure <- capture.output(summary(generator))
    #summary(generator)
    
    
    #generator.optimizer <- optimizer_adam( 
    #    lr = gene.learning.rate, 
    #    decay = gene.lr.decay)
    
    generator.optimizer <- tf$keras$optimizers$Adam(learning_rate = gene.learning.rate, beta_1 = 0.5)
    
    loss.zero <- function(y_true, y_pred){
        #Add orthonormal loss if use ddr
        if(useddr & (lambda.ortho.loss != 0)){
            weightmatrix <- reducer$weights[[1]]
            btb <- tf$matmul(tf$transpose(weightmatrix), weightmatrix)
            loss.orthonormal <- 
                tf$norm(btb - tf$linalg$diag(k_constant(1, dtype = "float32", shape = dr.nodes.num)))
        } else loss.orthonormal <- tf$constant(0, dtype = "float32")
        k_constant(lambda.ortho.loss , dtype = "float32") * loss.orthonormal
        #k_constant(0, dtype = "double")
    }
    
    #The kl loss part for both primal and dual form are the same
    if(v.label == "new"){
        loss.kl <- function(y_true, y_pred){
            k_mean(y_pred)}
    }
    
    if(v.label == "old"){
        if(kl.dual == TRUE) stop("If you want to use dual form, please use new label.")
        loss.kl <- function(y_true, y_pred){
            k_mean(-y_pred)
        }
    }
    
    if(lambda.mse == 0){
        loss.mse <- function(y_true, y_pred)
            k_constant(0, dtype = "float32")
    } else {
        loss.mse <- function(y_true, y_pred){
            k_mean(tf$math$squared_difference(y_true, y_pred)) * dim.y
        }
    }
    
    generator %>% compile(
        optimizer = generator.optimizer, loss = list(loss.zero, loss.kl),
        loss_weights = list(1, lambda.kl)
    )
    
    #*********************************   end of Define GENERATOR  ************************************
    
    
    
    
    ###################################    Define diagnostic picture  #####################################
    diagnostic.pic <- function(){
        x11()
        par(mfrow = c(3,2))
        #Predict conditional expectation E(Y|X) = g(X) on test set
        if(!is.null(diagnostic.info)){
            expect.y <- NULL
            for (j in 1:100) {
                diagonostic.xs <- array(0 , dim = c(trainningsetsize/100, m, dim.x))
                diagonostic.etas <- array(0 , dim = c(trainningsetsize/100, m, dim.eta))
                for(i in 1:(trainningsetsize/100)){
                    diagonostic.xs[i, ,1:dim.x] <- 
                        matrix(test.x[(i + trainningsetsize/100*(j-1)), ], byrow = TRUE, nrow = m, ncol = dim.x)
                    if(eta.dist == "normal"){
                        diagonostic.etas[i, ,] <- rnorm(m*dim.eta)
                    }
                    
                    if(eta.dist == "uniform"){
                        diagonostic.etas[i, ,] <- runif(m*dim.eta, 0,1)
                    }
                }
                expect.y <- c(expect.y, as.array(generator.g(list(diagonostic.xs, diagonostic.etas))))
            }
            
            plot(expect.y ~ diagnostic.info$test.y.mean,ylab = "Estimated Conditional expectation of Y", 
                 xlab = "True Conditional mean of y", main = "g_theta(X) VS E(Y|X) on Testset")
            abline(0,1, col = "red", lwd = 3)
            
            
        }
        
        #If true model not known ,show g vs true y.
        if(is.null(diagnostic.info)){
            expect.y <- NULL
            for (j in 1:100) {
                diagonostic.xs <- array(0 , dim = c(trainningsetsize/100, m, dim.x))
                diagonostic.etas <- array(0 , dim = c(trainningsetsize/100, m, dim.eta))
                for(i in 1:(trainningsetsize/100)){
                    diagonostic.xs[i, ,1:dim.x] <- 
                        matrix(test.x[(i + trainningsetsize/100*(j-1)), ], byrow = TRUE, nrow = m, ncol = dim.x)
                    if(eta.dist == "normal"){
                        diagonostic.etas[i, ,] <- rnorm(m*dim.eta)
                    }
                    
                    if(eta.dist == "uniform"){
                        diagonostic.etas[i, ,] <- runif(m*dim.eta,0,1)
                    }
                }
                expect.y <- c(expect.y, as.array(generator.g(list(diagonostic.xs, diagonostic.etas))))
            }
            plot(expect.y~test.y,ylab = "Estimated Conditional expectation of Y", 
                 xlab = "True single yi", main = "g_theta(X) VS y on Testset(True model unknown)")
            abline(0,1, col = "red", lwd = 3)
        }
        
        #plot the loss path
        plot(gene.loss[1:step,1]~c(1:step), type = "l", main = "Gene Losses", ylim = c(min(gene.loss),max(gene.loss)),
             xlab = "iterations", ylab = "Losses")
        lines((lambda.kl*gene.loss[1:step,2])~c(1:step), type = "l", col = "red")
        lines((lambda.mse*gene.loss[1:step,3]/2)~c(1:step), type = "l", lty = 2,col ="green")
        legend("topright",lty = c(1,1,2), col = c("black","red","green"), legend = c("Totol","KL","MSE"))
        
        #plot the disc loss
        plot(disc.loss, type = "l", main = "Discriminator loss", ylab = "loss", xlab = "epoch")
        
        #plot((lambda.kl*gene.loss[1:step,3])~c(1:step), type = "l", col = "red", 
        #     ylab = "KL loss", xlab = "iterations", main = "KL loss")
        
        
        # Marginal Density comparison
        if(!is.null(diagnostic.info)){
            density.estimated <- density(expect.y)
            plot(density(test.y), col = "red", ylim = c(0, max(density(test.y)$y,
                                                               density.estimated$y)),
                 lwd = 2,main = "Density comparison of y and g(X)")
            lines(density.estimated,
                  lwd = 2)
            legend("topright", lty = c(1,1), col = c("red","black"), legend = c("True y", "Estimated g(X)"))
        }
        
        
        
        #Pictures of conditional distribution
        x11(width = 30,height = 20, title = paste0("Iterations:",step))
        par(mfrow = c(5,4))
        for (i in 1:10) {
            random.x <- rnorm(dim.x)
            x.matrix <- matrix(random.x, nrow = 10000, ncol = dim.x, byrow = TRUE)
            if(eta.dist == "normal"){
                eta <- matrix(rnorm(dim.eta * 10000), nrow = 10000, ncol = dim.eta)
            }
            if(eta.dist == "uniform"){
                eta <- matrix(runif(dim.eta * 10000,0,1), nrow = 10000, ncol = dim.eta)
            }
            t <- as.array(generator.t(list(x.matrix, eta)))
            
            if(!is.null(diagnostic.info)){
                diag.data <- data.generate(n = 10000, x = random.x,
                                           dim.x = dim.x, dim.y = dim.y, dim.epsilon = diagnostic.info$dim.epsilon,
                                           epsilon.dist = diagnostic.info$epsilon.dist,
                                           model = diagnostic.info$model,
                                           seed = NULL,
                                           testdata = FALSE)
                y.true <- diag.data$y
            } else y.true <- NULL
            
            density.t <- density(t)
            density.y.true <- density(y.true)
            
            plot(density.t, lwd = 2, xlim = c(min(t, y.true), max(t,y.true)), 
                 ylim = c(0,max(density.t$y,density.y.true$y)), main = "Conditional Density")
            if(!is.null(diagnostic.info)) lines(density.y.true, col = "red", lwd = 2)
            y.mean <- signif(mean(y.true), digits = 3)
            t.mean <- signif(mean(t), digits = 3)
            y.sd   <- signif(sd(y.true), digits = 3)
            t.sd   <- signif(sd(t), digits = 3)
            legend("topright", lty = c(1,1), col = c("red","black"), 
                   legend = c(paste("y|x M:",y.mean, "SD:",y.sd), 
                              paste("T|x M:",t.mean, "SD:", t.sd)))
            
            
            #estimated density ratio
            conditional.r <- discriminator(list(x.matrix[1:5000, , drop=FALSE], 
                                                matrix(seq(min(t,y.true),max(t,y.true), length.out = 5000), nrow = 5000)))
            
            plot(exp(ifelse(v.label == "old", -1, 1)* as.array(conditional.r)) ~ 
                     seq(min(t, y.true),max(t,y.true), length.out = 5000), type = "l",
                 main = "Estimated conditional density ratio", ylab = "Estimated r", xlab = "y")
            
            
        }
        
    }
    ###############################   End of diagnostic picture   #########################
    
    
    
    
    ##############################   Diagnostic table  ####################################
    #Construct table
    diagnostic.table <- function(s, t.size = 5000, y.size = 5000){
        ks <- rep(0,s)
        p <- rep(0,s)
        true.point <- matrix(0, nrow = s, ncol = 7)
        dmr.point <- matrix(0, nrow = s, ncol = 7)
        quantile.proba <- c(0.05, 0.25, 0.5, 0.75, 0.95)
        set.seed(155)
        random.x.all <- matrix(rnorm(dim.x * s), nrow = s, ncol = dim.x)
        y.true.range <- matrix(0, nrow = s, ncol = 2)
        
        for(i in 1:s){
            random.x <- random.x.all[i, ]
            if(!is.null(diagnostic.info)){
                diag.data <- data.generate(n = y.size, x = random.x,
                                           dim.x = dim.x, dim.y = dim.y, 
                                           dim.epsilon = diagnostic.info$dim.epsilon,
                                           epsilon.dist = diagnostic.info$epsilon.dist,
                                           model = diagnostic.info$model,
                                           seed = NULL,
                                           testdata = FALSE,
                                           extrapara = diagnostic.info$extrapara)
                y.true <- diag.data$y
            } else y.true <- NULL
            x.matrix <- diag.data$x
            if(eta.dist == "normal"){
                eta <- matrix(rnorm(dim.eta * t.size), nrow = t.size, ncol = dim.eta)
            }
            if(eta.dist == "uniform"){
                eta <- matrix(runif(dim.eta * t.size,0,1), nrow = t.size, ncol = dim.eta)
            }
            
            t <- as.array(generator.t(list(x.matrix, eta)))
            
            
            #Store the range of true y value for use in next part: other methods
            y.true.range[i, ] <- range(y.true)
            #mean
            true.point[i,1] <- diag.data$y.mean[1,]
            #sd
            true.point[i,2] <- diag.data$y.sd[1]
            #quantiles
            true.point[i,3:7] <- quantile(y.true, quantile.proba)
            
            dmr.point[i,1] <- mean(t)
            dmr.point[i,2] <- sd(t)
            dmr.point[i,3:7] <- quantile(t, quantile.proba)
            
            #KS test:
            dmr.ks.test <- ks.test(y.true, t)
            ks[i] <- dmr.ks.test$statistic
            p[i] <- dmr.ks.test$p.value
            
        }
        table.mspe <- data.frame(GCDE = c(colMeans((dmr.point - true.point)^2), mean(ks),mean(p)))
        row.names(table.mspe) <- c("Mean", "SD", paste0("$\\tau=",as.character(quantile.proba), "$"), "ks","ks.p")
        table.mspe
    }
    
    
    ##############################  End of diagnostic table  ##############################
    
    
    
    
    
    
    
    
    
    
    ###################################################################################
    #**********************************************************************************
    ###############################   TRAINNING PROCESS   #######################
    #**********************************************************************************
    ###################################################################################
    #if continue is true, load the weight before training #
    continue <- FALSE
    if(!is.null(load.disc)){
        load_model_weights_hdf5(generator, load.gene)
        load_model_weights_hdf5(discriminator, load.disc)
    }
    # Log
    #dir.create("result", showWarnings = FALSE)
    cat("This is the log of trainning. ", as.character(Sys.time()), "\n", 
        file = paste0(ifelse(is.null(savepath), "./", savepath),"traininglosslog.txt"), 
        append = TRUE)
    
    stepfrom <- 1
    if(continue){
        stepfrom <- max(which(kl.loss[,1]!=0))
        if(is.infinite(stepfrom)){
            disc.loss <- NULL
            gene.loss <- matrix(0, nrow = iterations * step.gene, ncol = 3)
            stepfrom <- 1
        } else {
            gene.loss <- rbind(gene.loss[1:stepfrom, ],matrix(0, nrow = iterations * step.gene, ncol = 3))
            stepfrom <- stepfrom + 1
        }
    } else {
        gene.loss <- matrix(0, nrow = iterations * step.gene, ncol = 3)
        disc.loss <- NULL
    }
    
    time.start <- Sys.time()
    
    #create labels V for using in training discriminator:
    #New formulation: -1 for true, 1 for generated
    v.whole.labels <- rbind(matrix(1, nrow = trainningsetsize, ncol = 1),
                            matrix(-1, nrow = trainningsetsize, ncol = 1))
    v.labels <- rbind(matrix(1, nrow = batch_size, ncol = 1),
                      matrix(-1, nrow = batch_size, ncol = 1))
    
    if(v.label == "new"){
        v.whole.labels <- -v.whole.labels
        v.labels <- -v.labels
    }
    
    #create dummy labels for using in training generator:
    t.labels <- matrix(0, nrow = batch_size, ncol = 1)
    d.labels <- matrix(0, nrow = batch_size, ncol = 1)
    
    
    #Set sgd parameter
    sgd.n <- floor(trainningsetsize/batch_size)
    sgd.count.disc <- 1
    sgd.count.gene <- 1
    sgd.sample.disc <- sample(trainningsetsize, trainningsetsize)
    sgd.sample.gene <- sample(trainningsetsize, trainningsetsize)
    
    table.count <- 1
    table.list <- list()
    
    #*********** Main loop for training **********************************************************
    
    #Set progress bar
    pb <- progress_bar$new(format = " training [:bar] :percent in :elapsed  ETA: :eta",
                           total = (iterations-stepfrom +1), width = 70)
    for (step in stepfrom:iterations) {
        train.disc.count <- 1
        train.gene.count <- 1
        
        #Update discriminator by train on batch
        if(disc.train.type == "batch"){
            while (train.disc.count <= step.disc) {
                #get batch index
                sgdsample <- sgd.sample.disc[(1 + sgd.count.disc*batch_size - batch_size):(sgd.count.disc*batch_size)]
                sgd.count.disc <- sgd.count.disc + 1
                if(sgd.count.disc == sgd.n) 
                {sgd.count.disc <- 1
                sgd.sample.disc <- sample(trainningsetsize, trainningsetsize)}
                # generate random sample from y and y.tilde
                sgd.y <- train.y[sgdsample, , , ,drop = FALSE]
                sgd.x.disc <- train.x[sgdsample, , , ,drop = FALSE]
                #use the same batch of sgdsample as y?
                #Generate eta and y.tilde
                if(eta.dist == "normal"){
                    eta.tilde <- matrix(rnorm(batch_size*dim.eta), nrow = batch_size, ncol = dim.eta)}
                if(eta.dist == "uniform"){
                    eta.tilde <- matrix(runif(batch_size*dim.eta,0,1), nrow = batch_size, 
                                        ncol = dim.eta)}
                sgd.y.tilde <- as.array(predict(generator.t,
                                                list(sgd.x.disc, eta.tilde), training = FALSE))
                combined <- abind(sgd.y, sgd.y.tilde, along = 1)
                sgd.x.disc <- abind(sgd.x.disc, sgd.x.disc, along = 1)
                
                # Trains the discriminator
                d_loss <- discriminator %>% train_on_batch(list(sgd.x.disc, combined), v.labels)       #time comsuming
                train.disc.count <- train.disc.count + 1
                disc.loss <- c(disc.loss, d_loss)
            }
        }
        
        
        #Update discriminator by fit
        if(disc.train.type == "fit"){
            if(eta.dist == "normal"){
                eta.tilde <- matrix(rnorm(trainningsetsize*dim.eta), nrow = trainningsetsize, ncol = dim.eta)}
            if(eta.dist == "uniform"){
                eta.tilde <- matrix(runif(trainningsetsize*dim.eta,0,1), nrow = trainningsetsize, 
                                    ncol = dim.eta)}
            y.tilde <- as.array(generator.t(list(train.x, eta.tilde)))
            combined.y <- rbind(train.y, y.tilde)
            combined.x <- rbind(train.x, train.x)
            history <- fit(discriminator, list(combined.x, combined.y), v.whole.labels, epoch = 1, verbose = 0)
            disc.loss <- c(disc.loss, history$metrics$loss)
        }
        
        
        #Next, we train the generator
        while(train.gene.count <= step.gene){
            #get batch index
            sgdsample <- sgd.sample.gene[(1 + sgd.count.gene*batch_size - batch_size):(sgd.count.gene*batch_size)]
            sgd.count.gene <- sgd.count.gene + 1
            if(sgd.count.gene == sgd.n) 
            {sgd.count.gene <- 1
            sgd.sample.gene <- sample(trainningsetsize, trainningsetsize)}
            
            #Create sgd batch data
            sgd.x <- train.x[sgdsample, , , ,drop = FALSE]
            if(eta.dist == "normal"){
                sgd.eta <- matrix(rnorm(batch_size*dim.eta), nrow = batch_size, ncol = dim.eta)
            }
            if(eta.dist == "uniform"){
                sgd.eta <- matrix(runif(batch_size*dim.eta,0,1), nrow = batch_size, ncol = dim.eta)
            }
            #sgd.g.xs <- array(0 , dim = c(batch_size, m, dim.x))
            #sgd.g.etas <- array(0 , dim = c(batch_size, m, dim.eta))
            #for(i in 1:batch_size){
            #    sgd.g.xs[i, ,] <- matrix(sgd.x[i,], byrow = TRUE, nrow = m, ncol = dim.x)
            #    if(eta.dist == "normal"){
            #        sgd.g.etas[i, ,] <- rnorm(m*dim.eta)
            #    }
            #    if(eta.dist == "uniform"){
            #        sgd.g.etas[i, ,] <- runif(m*dim.eta,0,1)
            #    }
            #}
            sgd.y <- train.y[sgdsample, , , ,drop = FALSE]
            
            # Trains the generator (where the discriminator weights are frozen)
            g_loss <- generator %>% train_on_batch(list(sgd.x, sgd.eta), 
                                                   list(t.labels, d.labels))
            
            #save loss
            g.loss.temp <- unlist(g_loss)
            gene.loss[step,] <- c(g.loss.temp[1], g.loss.temp[3], 0)
            train.gene.count <- train.gene.count+1
        }
        #Free memory
        if(step %% 10 == 0) invisible(gc())
        
        
        # Occasionally adjust step.disc and saves models and print results
        if (step %% display == 0) {
            # Saves model weights
            #save_model_weights_hdf5(generator, "generator.h5")
            #save_model_weights_hdf5(discriminator, "discriminator.h5")
            # Prints metrics
            cat("Step:", step, " ")
            cat("Generator loss -", "Total:", gene.loss[step,1],
                " KL:", lambda.kl*gene.loss[step,2], 
                " LS:", lambda.mse/2 * gene.loss[step,3], " ")
            cat("Est:",
                format(round(difftime(Sys.time(), time.start, units = "mins")/step*(iterations-step))))
            print(Sys.time())
            cat("\n")
            #write log
            
            cat("Step:", step, " ",
                file = paste0(ifelse(is.null(savepath), "./", savepath),"traininglosslog.txt"), append = TRUE)
            cat("Generator loss -", "Total:", gene.loss[step,1],
                " KL:", lambda.kl* gene.loss[step,2], 
                " LS:", lambda.mse/2 * gene.loss[step,3], "  ",
                file = paste0(ifelse(is.null(savepath), "./", savepath),"traininglosslog.txt"), append = TRUE)
            cat("Est time:", 
                format(round(difftime(Sys.time(), time.start, units = "mins")/step*(iterations-step))), "  ",
                file = paste0(ifelse(is.null(savepath), "./", savepath),"traininglosslog.txt"), append = TRUE)
            cat(as.character(Sys.time()), "\n", 
                file = paste0(ifelse(is.null(savepath), "./", savepath),"traininglosslog.txt"), append = TRUE)
            
            #plot diagnostic image
            if(!is.null(diagnostic.info)){ 
                if(show.diagnostic.pic) diagnostic.pic()
                table.list[[table.count]] <- 
                    diagnostic.table(s=300, t.size = 2000, y.size = 2000)
                print(table.list[[table.count]])
                table.count <- table.count+1
            }
            #extra function to diagnostic while training
            if(!is.null(extra.function)) 
                extra.function(generator.t = generator.t, ite = step, dual = ifelse(kl.dual,"fenchel","primal"))
            
            #free memory
            invisible(gc())
        }
        pb$tick()
    }
    
    
    
    #********************************|    END OF  TRAINNING PROCESS   *****************************
    return(list(generator = generator, generator.t = generator.t,
                discriminator = discriminator,
                train.x            = train.x,
                train.y            = train.y,
                test.x             = test.x,
                test.y             = test.y,
                diagnostic.info    = diagnostic.info,
                v.label            = v.label,
                kl.dual            = kl.dual,
                iterations         = iterations,
                dim.x              = dim.x,
                dim.y              = dim.y,
                dim.eta            = dim.eta,
                eta.dist           = eta.dist,
                gene.nodes.num     = gene.nodes.num,
                dis.nodes.num      = dis.nodes.num,
                disc.train.type    = disc.train.type,           
                step.disc          = step.disc,              
                step.gene          = step.gene,                  
                display            = display,             
                batch_size         = batch_size,
                gene.learning.rate = gene.learning.rate,
                disc.learning.rate = disc.learning.rate,
                gene.lr.decay      = gene.lr.decay,            
                disc.lr.decay      = disc.lr.decay, 
                gene.regularizer.1 = gene.regularizer.1,   
                gene.regularizer.o = gene.regularizer.o,
                disc.regularizer   = disc.regularizer,   
                activation.disc    = activation.disc,
                activation.gene    = activation.gene,
                lambda.kl          = lambda.kl,
                lambda.mse         = lambda.mse,
                m                  = m,
                table.list         = table.list))
}

















#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
my_code_path <- "../../GCDE/code"
sapply(list.files(my_code_path, full.names = TRUE), source)

library(abind)
mnist <- dataset_mnist(path = "mnist.npz")
image(mnist$train$x[4,,])
mnist$train$y[1:4]

#x is the leftbottom 1/4 image. Y is the rest of the image.
x_train <- array(mnist$train$x[, 1:28, 1:14], dim = c(60000, 28,14,1))
x_train <- (x_train/255-0.5)/0.5
y_train <- array((mnist$train$x/255 - 0.5)/0.5, dim = c(60000, 28, 28, 1))


x_test <- array(mnist$test$x[, 1:28, 1:14], dim = c(10000, 28,14,1))
x_test <- (x_test/255 - 0.5)/0.5
y_test <- array((mnist$test$x/255 - 0.5)/0.5, dim = c(10000, 28, 28, 1))


#diagnos
#Extra diagnostic function 
plot.whiletraining.mnist.xhalfimage.cnn <- function(generator.t, ite, train.x = x_test, 
                                                    orig.image = mnist$test$x[1:62,,],
                                                labels = mnist$train$y, padding = "black", dual = "fenchel",
                                                pic.num = 19){
    #x11(width = 12, height = 8, title = paste0("iteration = ", ite))
    test.id <- c(4, 3, 2, 19, 5, 16, 12, 1, 62, 8)
    real.plot <- function(){
        par(mfrow = c(10,pic.num + 3), mar = c(0,0,0,0))
        for (i in test.id) {
            plot.new()
            plot.image <- matrix(ifelse(padding == "white", 255, -255), 32, 32)
            plot.image[3:30,3:16] <- train.x[i, , ,1] * 255
            
            n <- 1
            dim.eta  <- 100
            
            ori.im <- matrix(ifelse(padding == "white", 255, 0), 32, 32)
            ori.im[3:30,3:30] <- orig.image[i,,]
            image(t(ori.im)[, ncol(ori.im):1], axes = FALSE, 
                  col = gray.colors(255, start = 0, end = 1))
            rect(0,0,1,1, border = "white")
            #text(0,0.95, labels[i], cex=1, pos=4, col="red", offset=0)
            for(j in 1:pic.num){
                t <- as.array(predict(generator.t,
                                      list(array(train.x[i, , , ], dim = c(1,28, 14,1)), 
                                               matrix(rnorm(n*dim.eta), nrow = n, ncol = dim.eta)), 
                                          training = FALSE))
                plot.image[3:30, 17:30] <- t[1, 1:28, 15:28,1]*255
                if(j == 1){
                    plot.image[3:30, 17:30] <- -255
                }
                image(t(plot.image)[, ncol(plot.image):1], axes = FALSE, col = gray.colors(255, start = 0, end = 1))
                rect(0,0,0.5,1, border = "gray")
                #text(0,0.95, labels[i], cex=1, pos=4, col="red", offset=0)
            }
            plot.new()
        }
    }
    #dir.create("training_plot/xhalfimage")
    pdf(file = paste0("training_plot/xhalfimage/mnist_xhalfimage_cnn_ite=",ite, dual, ".pdf"), 
        width = pic.num+3,height = 10)
    real.plot()
    dev.off()
    
    real.plot()
}




#20000ite on lr 0.0002



#Start of the training
gcde.mnist.xhalfimage.cnn <- fitgcde.mnist.xhalfimage.cnn(x_train, y_train, x_test, y_test,
                                                  load.gene = "trained/mnist_xhalfimage_cnn_fenchel_generator.h5",
                                                  load.disc = "trained/mnist_xhalfimage_cnn_fenchel_discriminator.h5",
                                                  diagnostic.info    = NULL,
                                                  v.label            = "new",
                                                  kl.dual            = TRUE,
                                                  dual.form          = "fenchel",
                                                  iterations         = 1,
                                                  dim.eta            = 100,
                                                  eta.dist           = "normal",
                                                  gene.nodes.num     = c(200,200),
                                                  dis.nodes.num      = c(200,200),
                                                  disc.train.type    = "batch",            
                                                  step.disc          = 1,                  
                                                  step.gene          = 1,                   
                                                  display            = 4000,              
                                                  batch_size         = 128,
                                                  gene.learning.rate = 0.00004,
                                                  disc.learning.rate = 0.00004,
                                                  gene.lr.decay      = 0.0,              
                                                  disc.lr.decay      = 0.0, #1/(1+ decay*iteration)
                                                  gene.regularizer.1 = regularizer_l2(l = 0.000000),   
                                                  gene.regularizer.o = regularizer_l2(l = 0.0000),
                                                  disc.regularizer   = regularizer_l2(l = 0.00000),   
                                                  activation.disc    = "relu",
                                                  activation.gene    = "relu",
                                                  activation.gene.out= "sigmoid",
                                                  lambda.kl          = 20,
                                                  lambda.mse         = 0,       
                                                  m                  = 1,
                                                  extra.function     = plot.whiletraining.mnist.xhalfimage.cnn,
                                                  vanillaNN          = FALSE)
save_model_weights_hdf5(gcde.mnist.xhalfimage.cnn$generator, "trained/mnist_xhalfimage_cnn_fenchel_generator.h5")
save_model_weights_hdf5(gcde.mnist.xhalfimage.cnn$discriminator, "trained/mnist_xhalfimage_cnn_fenchel_discriminator.h5")

plot.whiletraining.mnist.xhalfimage.cnn(gcde.mnist.xhalfimage.cnn$generator.t, "pic19", dual = "fenchel",pic.num = 19)


sapply(5:10, function(i) plot.whiletraining.mnist.xhalfimage.cnn(
    gcde.mnist.xhalfimage.cnn$generator.t, paste0("pic",i), dual = "fenchel", pic.num = i))




#################################################################################################
#Primal form
gcde.mnist.xhalfimage.cnn.primal <- fitgcde.mnist.xhalfimage.cnn(x_train, y_train, x_test, y_test,
                                                         #load.gene = "trained/mnist_xhalfimage_cnn_primal_generator.h5",
                                                         #load.disc = "trained/mnist_xhalfimage_cnn_primal_discriminator.h5",
                                                         diagnostic.info    = NULL,
                                                         v.label            = "new",
                                                         kl.dual            = FALSE,
                                                         dual.form          = NULL,
                                                         iterations         = 20000,
                                                         dim.eta            = 100,
                                                         eta.dist           = "normal",
                                                         gene.nodes.num     = c(200,200),
                                                         dis.nodes.num      = c(200,200),
                                                         disc.train.type    = "batch",            
                                                         step.disc          = 1,                  
                                                         step.gene          = 1,                   
                                                         display            = 2000,              
                                                         batch_size         = 128,
                                                         gene.learning.rate = 0.0001,
                                                         disc.learning.rate = 0.0001,
                                                         gene.lr.decay      = 0.0,              
                                                         disc.lr.decay      = 0.0, #1/(1+ decay*iteration)
                                                         gene.regularizer.1 = regularizer_l2(l = 0.000000),   
                                                         gene.regularizer.o = regularizer_l2(l = 0.0000),
                                                         disc.regularizer   = regularizer_l2(l = 0.00000),   
                                                         activation.disc    = "relu",
                                                         activation.gene    = "relu",
                                                         activation.gene.out= "sigmoid",
                                                         lambda.kl          = 20,
                                                         lambda.mse         = 0,       
                                                         m                  = 1,
                                                         extra.function     = plot.whiletraining.mnist.xhalfimage.cnn,
                                                         vanillaNN          = FALSE)
save_model_weights_hdf5(gcde.mnist.xhalfimage.cnn.primal$generator, "trained/mnist_xhalfimage_cnn_primal_generator.h5")
save_model_weights_hdf5(gcde.mnist.xhalfimage.cnn.primal$discriminator, "trained/mnist_xhalfimage_cnn_primal_discriminator.h5")

plot.whiletraining.mnist.xhalfimage.cnn(gcde.mnist.xhalfimage.cnn.primal$generator.t, "primal test" )

sapply(5:10, function(i) plot.whiletraining.mnist.xhalfimage.cnn(
    gcde.mnist.xhalfimage.cnn$generator.t, paste0("pic",i), dual = "primal", pic.num = i))












