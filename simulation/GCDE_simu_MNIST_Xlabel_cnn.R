my_code_path <- "../../GCDE/code"
sapply(list.files(my_code_path, full.names = TRUE), source)

mnist <- dataset_mnist(path = "mnist.npz")
image(mnist$train$x[4,,])
mnist$train$y[1:4]

#X is the one hot label. Y is the image.
x_train <- array(0, dim = c(60000,10))
for(i in 1:60000) {
    x_train[i, mnist$train$y[i] + 1] <- 1}

y_train <- array((mnist$train$x/255 - 0.5)/0.5, dim = c(60000, 28, 28, 1))



#diagnos
#Extra diagnostic function 
plot.whiletraining.mnist.xlabel.cnn <- function(generator.t, ite, padding = "white", dual = "fenchel"){
    #x11(width = 12, height = 6, title = paste0("iteration = ", ite))
    
    real.plot <- function(){
        par(mfrow = c(10,14), mar = c(0,0,0,0))
    
        for (i in 1:10) {
            x.label <- array(0, dim = c(12, 10))
            x.label[ , i] <- 1
            
            n <- 12
            dim.eta  <- 100
            set.seed(131*i)
            t <- as.array(generator.t(list(x.label, 
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

plot.reals <- function(labels = mnist$train$y, ims = y_train, padding = "white"){
    pdf(file = paste0("training_plot/mnist_xlabel/mnist_xlabel_cnn_realimages.pdf"), width = 12,height = 10)
    par(mfrow = c(10,14), mar = c(0,0,0,0))
    
    for (i in 0:9) {
        index <- which(labels == i)
        
        for(j in 1:14){
            if(j == 1 || j == 14) {plot.new(); next}
            plot.image <- ims[index[j-1], , ,1] * 255
            im <- matrix(ifelse(padding == "white", 255, -255), 32, 32)
            im[3:30, 3:30] <- plot.image
            image(t(im)[, ncol(im):1], axes = FALSE, col = gray.colors(255, start = 0, end = 1))
        }
    }
    dev.off()
}


lr_schedule <- tf$keras$optimizers$schedules$ExponentialDecay(
    initial_learning_rate = 0.0002,
    decay_steps = 500L,
    decay_rate = 0.95,
    staircase = TRUE)




#Start of the training
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
#log: lr0.0001 for 30000 iterations

plot.whiletraining.mnist.xlabel.cnn(gcde.mnist$generator.t, 100)
#plot real images for comparison
plot.reals()

tt <- gcde.mnist$generator.t(list(matrix(x_train[1,], nrow = 2, ncol = 14^2), 
                                  matrix(rnorm(100), nrow = 2, ncol = 50)))





#################################################################################################
#Primal form
gcde.mnist.xlabel.primal <- fitgcde.mnist.xlabel.cnn(x_train, y_train, x_train, y_train,
                                       load.gene          = NULL,
                                       load.disc          = NULL,
                                       diagnostic.info    = NULL,
                                       v.label            = "new",
                                       kl.dual            = FALSE,
                                       dual.form          = NULL,
                                       iterations         = 20000,
                                       dim.eta            = 100,
                                       eta.dist           = "normal",
                                       gene.nodes.num     = c(500,500),
                                       dis.nodes.num      = c(500,500),
                                       disc.train.type    = "batch",            
                                       step.disc          = 1,                  
                                       step.gene          = 1,                   
                                       display            = 2000,              
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
save_model_weights_hdf5(gcde.mnist$generator, "trained/mnist_xlabel_cnn_primal_generator.h5")
save_model_weights_hdf5(gcde.mnist$discriminator, "trained/mnist_xlabel_cnn_primal_discriminator.h5")
#log: lr0.0001 for 30000 iterations

plot.whiletraining.mnist.xlabel.cnn(gcde.mnist$generator.t, 100)









