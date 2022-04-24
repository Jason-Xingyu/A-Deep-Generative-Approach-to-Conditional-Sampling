my_code_path <- "../../GCDE/code"
sapply(list.files(my_code_path, full.names = TRUE), source)

mnist <- dataset_mnist(path = "mnist.npz")
image(mnist$train$x[4,,])
mnist$train$y[1:4]

#x is the leftbottom 1/4 image. Y is the rest of the image.
x_train <- array(mnist$train$x[, 15:28, 1:14], dim = c(60000, 14^2))
x_train <- (x_train/255-0.5)/0.5
y_train <- cbind(array(mnist$train$x[, 1:14, 1:28], dim = c(60000, 14*28)), 
                 array(mnist$train$x[, 15:28, 15:28], dim = c(60000, 14^2)))
y_train <- (y_train/255 - 0.5)/0.5

x_test <- array(mnist$test$x[, 15:28, 1:14], dim = c(10000, 14^2))
x_test <- (x_test/255 - 0.5)/0.5
y_test <- cbind(array(mnist$test$x[, 1:14, 1:28], dim = c(10000, 14*28)), 
                array(mnist$test$x[, 15:28, 15:28], dim = c(10000, 14^2)))
y_test <- (y_test/255 - 0.5)/0.5

#diagnos
#Extra diagnostic function 
plot.whiletraining.mnist <- function(generator.t, ite, train.x = x_test, orig.image = mnist$test$x[1:40,,],
                                     labels = mnist$train$y, padding = "black"){
    x11(width = 12, height = 8, title = paste0("iteration = ", ite))
    par(mfrow = c(12,20), mar = c(0,0,0,0))
    
    for (i in 1:12) {
        plot.image <- matrix(ifelse(padding == "white", 255, -255), 32, 32)
        plot.image[17:30,3:16] <- train.x[i, ] * 255
        
        n <- 1
        dim.eta  <- 100
        
        ori.im <- matrix(ifelse(padding == "white", 255, 0), 32, 32)
        ori.im[3:30,3:30] <- orig.image[i,,]
        image(t(ori.im)[, ncol(ori.im):1], axes = FALSE, 
              col = gray.colors(255, start = 0, end = 1))
        rect(0,0,0.5,0.5, border = "gray")
        #text(0,0.95, labels[i], cex=1, pos=4, col="red", offset=0)
        for(j in 1:19){
            t <- as.array(generator.t(list(matrix(train.x[i,], nrow = n, ncol = 14^2), 
                                           matrix(rnorm(n*dim.eta), nrow = n, ncol = dim.eta))))
            plot.image[3:16, 3:30] <- t[1, 1:(14*28)]*255
            plot.image[17:30, 17:30] <- t[1, (14*28+1):(14*28+14^2)]*255
            image(t(plot.image)[, ncol(plot.image):1], axes = FALSE, col = gray.colors(255, start = 0, end = 1))
            rect(0,0,0.5,0.5, border = "gray")
            #text(0,0.95, labels[i], cex=1, pos=4, col="red", offset=0)
        }
    }
}


lr_schedule <- tf$keras$optimizers$schedules$ExponentialDecay(
    initial_learning_rate = 0.00005,
    decay_steps = 10000L,
    decay_rate = 0.9,
    staircase = FALSE)


#60000ite on lr 0.0002

#Start of the training
gcde.mnist.ximage <- fitgcde.mnist.ximage(x_train, y_train, x_test, y_test,
                                   load.gene = "trained/mnist_ximage_fenchel_generator.h5",
                                   load.disc = "trained/mnist_ximage_fenchel_discriminator.h5",
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
                      display            = 2500,              
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
                      extra.function     = plot.whiletraining.mnist,
                      vanillaNN          = FALSE)
save_model_weights_hdf5(gcde.mnist$generator, "trained/mnist_ximage_fenchel_generator.h5")
save_model_weights_hdf5(gcde.mnist$discriminator, "trained/mnist_ximage_fenchel_discriminator.h5")




plot.whiletraining.mnist(gcde.mnist.ximage$generator.t,100 )














