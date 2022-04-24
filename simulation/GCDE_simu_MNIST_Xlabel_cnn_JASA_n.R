####################################################################################################
##############  2021.6.20 JASA Review Response
####################################################################################################
### Part 1: different training set size
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
            dim.eta  <- 100
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


x_train.n1 <- x_train[1:2000, ]
y_train.n1 <- y_train[1:2000, , , , drop=FALSE]
x_train.n2 <- x_train[1:5000, ]
y_train.n2 <- y_train[1:5000, , , , drop=FALSE]
x_train.n3 <- x_train[1:10000, ]
y_train.n3 <- y_train[1:10000, , , , drop=FALSE]

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





#n1, smallest training set.
gcde.mnist.n1 <- fitgcde.mnist.xlabel.cnn(x_train.n1, y_train.n1, NULL, NULL,
                                           load.gene          = "trained/mnist_xlabel_cnn_fenchel_jasan1_generator.h5",
                                           load.disc          = "trained/mnist_xlabel_cnn_fenchel_jasan1_discriminator.h5",
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
save_model_weights_hdf5(gcde.mnist.n1$generator, "trained/mnist_xlabel_cnn_fenchel_jasan1_generator.h5")
save_model_weights_hdf5(gcde.mnist.n1$discriminator, "trained/mnist_xlabel_cnn_fenchel_jasan1_discriminator.h5")
#log: lr0.0001 for 30000 iterations

plot.whiletraining.mnist.xlabel.cnn(gcde.mnist.n1$generator.t, 100)




#n2
gcde.mnist.n2 <- fitgcde.mnist.xlabel.cnn(x_train.n2, y_train.n2, NULL, NULL,
                                          load.gene          = "trained/mnist_xlabel_cnn_fenchel_jasan2_generator.h5",
                                          load.disc          = "trained/mnist_xlabel_cnn_fenchel_jasan2_discriminator.h5",
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
save_model_weights_hdf5(gcde.mnist.n2$generator, "trained/mnist_xlabel_cnn_fenchel_jasan2_generator.h5")
save_model_weights_hdf5(gcde.mnist.n2$discriminator, "trained/mnist_xlabel_cnn_fenchel_jasan2_discriminator.h5")
#log: lr0.0001 for 30000 iterations

plot.whiletraining.mnist.xlabel.cnn(gcde.mnist.n2$generator.t, 100)



#n3
gcde.mnist.n3 <- fitgcde.mnist.xlabel.cnn(x_train.n3, y_train.n3, NULL, NULL,
                                          load.gene          = "trained/mnist_xlabel_cnn_fenchel_jasan3_generator.h5",
                                          load.disc          = "trained/mnist_xlabel_cnn_fenchel_jasan3_discriminator.h5",
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
save_model_weights_hdf5(gcde.mnist.n3$generator, "trained/mnist_xlabel_cnn_fenchel_jasan3_generator.h5")
save_model_weights_hdf5(gcde.mnist.n3$discriminator, "trained/mnist_xlabel_cnn_fenchel_jasan3_discriminator.h5")
#log: lr0.0001 for 30000 iterations

plot.whiletraining.mnist.xlabel.cnn(gcde.mnist.n3$generator.t, 100)



plot.mnist.multi.n <- function(gene1, gene2, gene3, padding = "white"){
    #x11(width = 12, height = 6, title = paste0("iteration = ", ite))
    
    real.plot <- function(){
        par(mfrow = c(10,20), mar = c(0,0,0,0))
        
        for (i in 1:10) {
            x.label <- array(0, dim = c(6, 10))
            x.label[ , i] <- 1
            
            n <- 6
            dim.eta  <- 100
            set.seed(121*i)
            for (k in 1:3) {
                if(k==1) t <- as.array(gene1(list(x.label, 
                                        array(rnorm(n*dim.eta), dim = c(n, 1, 1, dim.eta))), training = FALSE))
                if(k==2) t <- as.array(gene2(list(x.label, 
                                        array(rnorm(n*dim.eta), dim = c(n, 1, 1, dim.eta))), training = FALSE))
                if(k==3) t <- as.array(gene3(list(x.label, 
                                        array(rnorm(n*dim.eta), dim = c(n, 1, 1, dim.eta))), training = FALSE))
                
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
    
    pdf(file = paste0("training_plot/mnist_xlabel/mnist_jasa_multi_n.pdf"),
        width = 23,height = 10)
    real.plot()
    dev.off()
    
    real.plot()
}



#plot to compare different models trained by 2000,5000, 10000,60000 images.
plot.mnist.multi.n(gcde.mnist.n2$generator.t, 
                   gcde.mnist.n3$generator.t, gcde.mnist$generator.t)



#Calculate FID:
# get the # of digits in testing dataset.
digi <- c(table(mnist$test$y))

generate.jasa.set <- function(fit){
    res <- NULL
    
    for (i in 1:10) {
        x.label <- array(0, dim = c(digi[i], 10))
        x.label[ , i] <- 1
        
        n <- digi[i]
        dim.eta  <- 100
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


set.n1 <- generate.jasa.set(gcde.mnist.n1)
set.n2 <- generate.jasa.set(gcde.mnist.n2)
set.n3 <- generate.jasa.set(gcde.mnist.n3)
set.n <- generate.jasa.set(gcde.mnist)

true.set <- array(0, dim = c(10000, 28, 28, 3))
true.set[,,,1] <- mnist$test$x
true.set[,,,2] <- true.set[,,,1]
true.set[,,,3] <- true.set[,,,1]
#FID and scale_images function is in GCDE-GIT/code
fid.n1 <- fid(set.n1, true.set)
fid.n2 <- fid(set.n2, true.set)
fid.n3 <- fid(set.n3, true.set)
fid.n <- fid(set.n, true.set)

base::save(fid.n1,fid.n2, fid.n3, fid.n, file = "trained/mnist_jasa_n.RData")

fid.n1
fid.n2
fid.n3
fid.n
