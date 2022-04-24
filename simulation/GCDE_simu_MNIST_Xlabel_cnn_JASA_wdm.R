####################################################################################################
##############  2021.6.20 JASA Review Response
####################################################################################################
### Part 2: different width,depth, m
my_code_path <- "../../GCDE/code"
sapply(list.files(my_code_path, full.names = TRUE), source)

#diagnos
#Extra diagnostic function 
plot.whiletraining.mnist.xlabel.cnn <- function(fit, ite, padding = "white", dual = "fenchel"){
    #x11(width = 12, height = 6, title = paste0("iteration = ", ite))
    
    real.plot <- function(){
        par(mfrow = c(10,14), mar = c(0,0,0,0))
        
        for (i in 1:10) {
            x.label <- array(0, dim = c(12, 10))
            x.label[ , i] <- 1
            
            n <- 12
            dim.eta  <- fit$dim.eta
            set.seed(131*i)
            t <- as.array(fit$generator.t(list(x.label, 
                                           array(rnorm(n*dim.eta), dim = c(n, 1, 1, dim.eta))), training = FALSE))
            for(j in 1:14){
                if(j == 1 || j == 14) {plot.new(); next}
                plot.image <- array(t[j-1, , ,1] * 255, dim = c(28,28))
                im <- matrix(ifelse(padding == "white", 255, -255), 32, 32)
                im[3:30, 3:30] <- plot.image
                image(t(im)[, ncol(im):1], axes = FALSE, col = gray.colors(255, start = 0, end = 1))
            }
        }
    }
    
    pdf(file = paste0("training_plot/mnist_xlabel/mnist_xlabel_cnn_ite=",ite, dual, ".pdf"), width = 14,height = 10)
    real.plot()
    dev.off()
    
    real.plot()
}

mnist <- dataset_mnist(path = "mnist.npz")
image(mnist$train$x[4,,])
mnist$train$y[1:4]

#X is the one hot label. Y is the image.
x_train <- array(0, dim = c(60000,10))
for(i in 1:60000) {
    x_train[i, mnist$train$y[i] + 1] <- 1}

y_train <- array((mnist$train$x/255 - 0.5)/0.5, dim = c(60000, 28, 28, 1))



#Trained by full dataset.
gcde.mnist <- fitgcde.mnist.xlabel.cnn(x_train, y_train, x_train, y_train,
                                       load.gene          = "trained/mnist_xlabel_cnn_fenchel_generator.h5",
                                       load.disc          = "trained/mnist_xlabel_cnn_fenchel_discriminator.h5",
                                       diagnostic.info    = NULL,
                                       v.label            = "new",
                                       kl.dual            = TRUE,
                                       dual.form          = "fenchel",
                                       iterations         = 1,
                                       dim.eta            = 100,
                                       eta.dist           = "normal",
                                       gene.nodes.num     = c(500,500),
                                       dis.nodes.num      = c(500,500),
                                       disc.train.type    = "batch",            
                                       step.disc          = 1,                  
                                       step.gene          = 1,                   
                                       display            = 1000,              
                                       batch_size         = 128,
                                       gene.learning.rate = 0.0001,
                                       disc.learning.rate = 0.0001,
                                       gene.lr.decay      = 0.0,              
                                       disc.lr.decay      = 0.0, #1/(1+ decay*iteration)
                                       gene.regularizer.1 = regularizer_l2(l = 0.00000),   
                                       gene.regularizer.o = regularizer_l2(l = 0.0000),
                                       disc.regularizer   = regularizer_l2(l = 0.00000),   
                                       activation.disc    = "relu",
                                       activation.gene    = "relu",
                                       activation.gene.out= "sigmoid",
                                       dropout.rate       = 0.5,
                                       lambda.kl          = 20,
                                       lambda.mse         = 0,       
                                       m                  = 1,
                                       extra.function     = plot.whiletraining.mnist.xlabel.cnn,
                                       vanillaNN          = FALSE)
save_model_weights_hdf5(gcde.mnist$generator, "trained/mnist_xlabel_cnn_fenchel_generator.h5")
save_model_weights_hdf5(gcde.mnist$discriminator, "trained/mnist_xlabel_cnn_fenchel_discriminator.h5")





#width and dapth halved
gcde.mnist.wdhalf <- fitgcde.mnist.xlabel.cnn(x_train, y_train, NULL, NULL,
                                              widthanddepth = 0.5,
                                          load.gene          = "trained/mnist_xlabel_cnn_fenchel_jasawdhalf_generator.h5",
                                          load.disc          = "trained/mnist_xlabel_cnn_fenchel_jasawdhalf_discriminator.h5",
                                          diagnostic.info    = NULL,
                                          v.label            = "new",
                                          kl.dual            = TRUE,
                                          dual.form          = "fenchel",
                                          iterations         = 2,
                                          dim.eta            = 100,
                                          eta.dist           = "normal",
                                          gene.nodes.num     = c(500,500),
                                          dis.nodes.num      = c(500,500),
                                          disc.train.type    = "batch",            
                                          step.disc          = 1,                  
                                          step.gene          = 1,                   
                                          display            = 500,              
                                          batch_size         = 128,
                                          gene.learning.rate = 0.0001,
                                          disc.learning.rate = 0.0001,
                                          gene.lr.decay      = 0.0,              
                                          disc.lr.decay      = 0.0, #1/(1+ decay*iteration)
                                          gene.regularizer.1 = regularizer_l2(l = 0.00000),   
                                          gene.regularizer.o = regularizer_l2(l = 0.0000),
                                          disc.regularizer   = regularizer_l2(l = 0.00000),   
                                          activation.disc    = "relu",
                                          activation.gene    = "relu",
                                          activation.gene.out= "sigmoid",
                                          dropout.rate       = 0.5,
                                          lambda.kl          = 20,
                                          lambda.mse         = 0,       
                                          m                  = 1,
                                          extra.function     = plot.whiletraining.mnist.xlabel.cnn,
                                          vanillaNN          = FALSE)
save_model_weights_hdf5(gcde.mnist.wdhalf$generator, "trained/mnist_xlabel_cnn_fenchel_jasawdhalf_generator.h5")
save_model_weights_hdf5(gcde.mnist.wdhalf$discriminator, "trained/mnist_xlabel_cnn_fenchel_jasawdhalf_discriminator.h5")
#log: lr0.0001 for 30000 iterations

plot.whiletraining.mnist.xlabel.cnn(gcde.mnist.wdhalf$generator.t, 100)




#wddouble
gcde.mnist.wddouble <- fitgcde.mnist.xlabel.cnn(x_train, y_train, NULL, NULL,
                                                widthanddepth = 2,
                                          load.gene          = "trained/mnist_xlabel_cnn_fenchel_jasawddouble_generator.h5",
                                          load.disc          = "trained/mnist_xlabel_cnn_fenchel_jasawddouble_discriminator.h5",
                                          diagnostic.info    = NULL,
                                          v.label            = "new",
                                          kl.dual            = TRUE,
                                          dual.form          = "fenchel",
                                          iterations         = 1,
                                          dim.eta            = 100,
                                          eta.dist           = "normal",
                                          gene.nodes.num     = c(500,500),
                                          dis.nodes.num      = c(500,500),
                                          disc.train.type    = "batch",            
                                          step.disc          = 1,                  
                                          step.gene          = 1,                   
                                          display            = 500,              
                                          batch_size         = 128,
                                          gene.learning.rate = 0.0001,
                                          disc.learning.rate = 0.0001,
                                          gene.lr.decay      = 0.0,              
                                          disc.lr.decay      = 0.0, #1/(1+ decay*iteration)
                                          gene.regularizer.1 = regularizer_l2(l = 0.00000),   
                                          gene.regularizer.o = regularizer_l2(l = 0.0000),
                                          disc.regularizer   = regularizer_l2(l = 0.00000),   
                                          activation.disc    = "relu",
                                          activation.gene    = "relu",
                                          activation.gene.out= "sigmoid",
                                          dropout.rate       = 0.5,
                                          lambda.kl          = 20,
                                          lambda.mse         = 0,       
                                          m                  = 1,
                                          extra.function     = plot.whiletraining.mnist.xlabel.cnn,
                                          vanillaNN          = FALSE)
save_model_weights_hdf5(gcde.mnist.wddouble$generator, "trained/mnist_xlabel_cnn_fenchel_jasawddouble_generator.h5")
save_model_weights_hdf5(gcde.mnist.wddouble$discriminator, "trained/mnist_xlabel_cnn_fenchel_jasawddouble_discriminator.h5")
#log: lr0.0001 for 30000 iterations

plot.whiletraining.mnist.xlabel.cnn(gcde.mnist.wddouble, 100)



#m10
gcde.mnist.m10 <- fitgcde.mnist.xlabel.cnn(x_train, y_train, NULL, NULL,
                                          load.gene          = "trained/mnist_xlabel_cnn_fenchel_jasam10_generator.h5",
                                          load.disc          = "trained/mnist_xlabel_cnn_fenchel_jasam10_discriminator.h5",
                                          diagnostic.info    = NULL,
                                          v.label            = "new",
                                          kl.dual            = TRUE,
                                          dual.form          = "fenchel",
                                          iterations         = 1,
                                          dim.eta            = 10,
                                          eta.dist           = "normal",
                                          gene.nodes.num     = c(500,500),
                                          dis.nodes.num      = c(500,500),
                                          disc.train.type    = "batch",            
                                          step.disc          = 1,                  
                                          step.gene          = 1,                   
                                          display            = 500,              
                                          batch_size         = 128,
                                          gene.learning.rate = 0.0001,
                                          disc.learning.rate = 0.0001,
                                          gene.lr.decay      = 0.0,              
                                          disc.lr.decay      = 0.0, #1/(1+ decay*iteration)
                                          gene.regularizer.1 = regularizer_l2(l = 0.00000),   
                                          gene.regularizer.o = regularizer_l2(l = 0.0000),
                                          disc.regularizer   = regularizer_l2(l = 0.00000),   
                                          activation.disc    = "relu",
                                          activation.gene    = "relu",
                                          activation.gene.out= "sigmoid",
                                          dropout.rate       = 0.5,
                                          lambda.kl          = 20,
                                          lambda.mse         = 0,       
                                          m                  = 1,
                                          extra.function     = plot.whiletraining.mnist.xlabel.cnn,
                                          vanillaNN          = FALSE)
save_model_weights_hdf5(gcde.mnist.m10$generator, "trained/mnist_xlabel_cnn_fenchel_jasam10_generator.h5")
save_model_weights_hdf5(gcde.mnist.m10$discriminator, "trained/mnist_xlabel_cnn_fenchel_jasam10_discriminator.h5")
#log: lr0.0001 for 30000 iterations

plot.whiletraining.mnist.xlabel.cnn(gcde.mnist.m10$generator.t, 100)


#mdouble
gcde.mnist.mdouble <- fitgcde.mnist.xlabel.cnn(x_train, y_train, NULL, NULL,
                                             load.gene          = "trained/mnist_xlabel_cnn_fenchel_jasamdouble_generator.h5",
                                             load.disc          = "trained/mnist_xlabel_cnn_fenchel_jasamdouble_discriminator.h5",
                                             diagnostic.info    = NULL,
                                             v.label            = "new",
                                             kl.dual            = TRUE,
                                             dual.form          = "fenchel",
                                             iterations         = 2,
                                             dim.eta            = 200,
                                             eta.dist           = "normal",
                                             gene.nodes.num     = c(500,500),
                                             dis.nodes.num      = c(500,500),
                                             disc.train.type    = "batch",            
                                             step.disc          = 1,                  
                                             step.gene          = 1,                   
                                             display            = 500,              
                                             batch_size         = 128,
                                             gene.learning.rate = 0.0001,
                                             disc.learning.rate = 0.0001,
                                             gene.lr.decay      = 0.0,              
                                             disc.lr.decay      = 0.0, #1/(1+ decay*iteration)
                                             gene.regularizer.1 = regularizer_l2(l = 0.00000),   
                                             gene.regularizer.o = regularizer_l2(l = 0.0000),
                                             disc.regularizer   = regularizer_l2(l = 0.00000),   
                                             activation.disc    = "relu",
                                             activation.gene    = "relu",
                                             activation.gene.out= "sigmoid",
                                             dropout.rate       = 0.5,
                                             lambda.kl          = 20,
                                             lambda.mse         = 0,       
                                             m                  = 1,
                                             extra.function     = plot.whiletraining.mnist.xlabel.cnn,
                                             vanillaNN          = FALSE)
save_model_weights_hdf5(gcde.mnist.mdouble$generator, "trained/mnist_xlabel_cnn_fenchel_jasamdouble_generator.h5")
save_model_weights_hdf5(gcde.mnist.mdouble$discriminator, "trained/mnist_xlabel_cnn_fenchel_jasamdouble_discriminator.h5")
#log: lr0.0001 for 30000 iterations

plot.whiletraining.mnist.xlabel.cnn(gcde.mnist.mdouble$generator.t, 100)






plot.mnist.jasa <- function(gene1, gene2, gene3, padding = "white", title, seed = NULL){
    #x11(width = 12, height = 6, title = paste0("iteration = ", ite))
    
    real.plot <- function(){
        par(mfrow = c(10,20), mar = c(0,0,0,0))
        
        for (i in 1:10) {
            n <- 6
            x.label <- array(0, dim = c(n, 10))
            x.label[ , i] <- 1
        
            if(is.null(seed)) set.seed(131*i)
            else set.seed(seed)
            for (k in 1:3) {
                if(k==1) {
                    dim.eta <- gene1$dim.eta
                    t <- as.array(gene1$generator.t(list(x.label, 
                                                  array(rnorm(n*dim.eta), dim = c(n, 1, 1, dim.eta))), training = FALSE))
                    }
                if(k==2) {
                    dim.eta <- gene2$dim.eta
                    t <- as.array(gene2$generator.t(list(x.label, 
                                                  array(rnorm(n*dim.eta), dim = c(n, 1, 1, dim.eta))), training = FALSE))
                    }
                if(k==3) {
                    dim.eta <- gene3$dim.eta
                    t <- as.array(gene3$generator.t(list(x.label, 
                                                  array(rnorm(n*dim.eta), dim = c(n, 1, 1, dim.eta))), training = FALSE))
                }
                
                for(j in 1:6){
                    #plot.new()
                    plot.image <- array(t[j, , ,1] * 255, dim = c(28,28))
                    im <- matrix(ifelse(padding == "white", 255, -255), 32, 32)
                    im[3:30, 3:30] <- plot.image
                    image(t(im)[, ncol(im):1], axes = FALSE, col = gray.colors(255, start = 0, end = 1))
                }
                if(k!=3) plot.new()
            }
            
            
        }
    }
    
    pdf(file = paste0("training_plot/mnist_xlabel/mnist_jasa_",title,".pdf"),
        width = 20,height = 10)
    real.plot()
    dev.off()
    
    real.plot()
}




#Plot of comparing different width and depth.
plot.mnist.jasa(gcde.mnist.wdhalf, gcde.mnist, 
                   gcde.mnist.wddouble, title = "widthanddepth")


#Plot of comparing different m
plot.mnist.jasa(gcde.mnist.m10, gcde.mnist, 
                gcde.mnist.mdouble, title = "m", seed =  3)



#Calculate FID:
# get the # of digits in testing dataset.
digi <- c(table(mnist$test$y))

generate.jasa.set <- function(fit){
    res <- NULL
    
    for (i in 1:10) {
        x.label <- array(0, dim = c(digi[i], 10))
        x.label[ , i] <- 1
        
        n <- digi[i]
        dim.eta  <- fit$dim.eta
        set.seed(131*i)
        
        t <- as.array(fit$generator.t(list(x.label, 
                                           array(rnorm(n*dim.eta), dim = c(n, 1, 1, dim.eta))), training = FALSE))
        
        res <- abind(res, t, along = 1)
    }
    
    res <- (res/2 + 0.5)*255
    result <- array(0, dim = c(10000, 28, 28, 3))
    result[, , ,1] <- res
    result[, , ,2] <- res
    result[, , ,3] <- res
    result
}


set.wdh <- generate.jasa.set(gcde.mnist.wdhalf)
set.wdd <- generate.jasa.set(gcde.mnist.wddouble)
set.m10 <- generate.jasa.set(gcde.mnist.m10)
set.m200 <- generate.jasa.set(gcde.mnist.mdouble)
set.base <- generate.jasa.set(gcde.mnist)

true.set <- array(0, dim = c(10000, 28, 28, 3))
true.set[,,,1] <- mnist$test$x
true.set[,,,2] <- true.set[,,,1]
true.set[,,,3] <- true.set[,,,1]
#FID and scale_images function is in GCDE-GIT/code
fid.wdh <- fid(set.wdh, true.set)
fid.wdd <- fid(set.wdd, true.set)
fid.m10 <- fid(set.m10, true.set)
fid.m200 <- fid(set.m200, true.set)
fid.base <- fid(set.base, true.set)


base::save(fid.wdh,fid.wdd, fid.m10, fid.m200, fid.base, file = "trained/mnist_jasa_wdm.RData")

fid.wdh
fid.base
fid.wdd

fid.m10
fid.base
fid.m200
