source("D:/Dropbox/UIowa\ Study\ Materials/Ph.D\ research/GCDE/GCDE.R")
#######################################################################################################
############################     Model T1           ###################################################
#######################################################################################################
source("D:/Dropbox/UIowa\ Study\ Materials/Ph.D\ research/GCDE/GCDE.R")
show.diagnostic.pic <- FALSE
modelname <- "t1"
savepath <- "D:/Dropbox/UIowa Study Materials/Ph.D research/GCDE/code_simulation_0/trained/"



t1fun <- function(r){
    seed <- 131 * r
    cat("Start the", r, "simulation replicates...")
    data.t1 <- data.generate(n = 5000,
                             x = NULL,
                             dim.x = 5,
                             dim.y = 1,
                             dim.epsilon = 1,
                             epsilon.dist = "normal",
                             model = "add2",
                             seed = seed,
                             testdata = TRUE)
    
    
    
    gcde.t1.primal <- fitgcde(data.t1$x, data.t1$y, data.t1$test.x, data.t1$test.y,
                              diagnostic.info    = data.t1$diagnostic.info,
                              show.diagnostic.pic= show.diagnostic.pic,
                              v.label            = "new",
                              kl.dual            = FALSE,
                              iterations         = 1500,
                              dim.eta            = 3,
                              eta.dist           = "normal",
                              gene.nodes.num     = c(30),
                              dis.nodes.num      = c(40,20),
                              disc.train.type    = "fit",            
                              step.disc          = 10,                  
                              step.gene          = 5,                   
                              display            = 300,               
                              batch_size         = 128,
                              gene.learning.rate = 0.001,
                              disc.learning.rate = 0.0001,
                              gene.lr.decay      = 0,              
                              disc.lr.decay      = 0, 
                              gene.regularizer.1 = regularizer_l2(l = 0.000005),   
                              gene.regularizer.o = regularizer_l2(l = 0.000005),
                              disc.regularizer   = regularizer_l2(l = 0.0000),   
                              activation.disc    = "relu",
                              activation.gene    = "relu",
                              lambda.kl          = 20,
                              lambda.mse         = 0,       
                              m                  = 1, vanillaNN = FALSE)
    
    gcde.t1.fenchel <- fitgcde(data.t1$x, data.t1$y, data.t1$test.x, data.t1$test.y,
                            diagnostic.info    = data.t1$diagnostic.info,
                            show.diagnostic.pic= show.diagnostic.pic,
                            v.label            = "new",
                            kl.dual            = TRUE,
                            iterations         = 7500,
                            dim.eta            = 3,
                            eta.dist           = "normal",
                            gene.nodes.num     = c(30),
                            dis.nodes.num      = c(40,20),
                            disc.train.type    = "batch",            
                            step.disc          = 5,                  
                            step.gene          = 1,                   
                            display            = 1500,               
                            batch_size         = 128,
                            gene.learning.rate = 0.001,
                            disc.learning.rate = 0.0001,
                            gene.lr.decay      = 0,              
                            disc.lr.decay      = 0, 
                            gene.regularizer.1 = regularizer_l2(l = 0.000005),   
                            gene.regularizer.o = regularizer_l2(l = 0.000005),
                            disc.regularizer   = regularizer_l2(l = 0.0000),   
                            activation.disc    = "relu",
                            activation.gene    = "relu",
                            lambda.kl          = 20,
                            lambda.mse         = 0,       
                            m                  = 1, vanillaNN = FALSE)
    
    gcde.t1.donsker <- fitgcde(data.t1$x, data.t1$y, data.t1$test.x, data.t1$test.y,
                            diagnostic.info    = data.t1$diagnostic.info,
                            show.diagnostic.pic= show.diagnostic.pic,
                            v.label            = "new",
                            kl.dual            = TRUE,
                            dual.form          = "donsker",
                            iterations         = 7500,
                            dim.eta            = 3,
                            eta.dist           = "normal",
                            gene.nodes.num     = c(30),
                            dis.nodes.num      = c(40,20),
                            disc.train.type    = "batch",            
                            step.disc          = 5,                  
                            step.gene          = 1,                   
                            display            = 1500,               
                            batch_size         = 128,
                            gene.learning.rate = 0.001,
                            disc.learning.rate = 0.0001,
                            gene.lr.decay      = 0,              
                            disc.lr.decay      = 0, 
                            gene.regularizer.1 = regularizer_l2(l = 0.000005),   
                            gene.regularizer.o = regularizer_l2(l = 0.000005),
                            disc.regularizer   = regularizer_l2(l = 0.0000),   
                            activation.disc    = "relu",
                            activation.gene    = "relu",
                            lambda.kl          = 20,
                            lambda.mse         = 0,       
                            m                  = 1, vanillaNN = FALSE)
    gcde.save.model.weights(list(gcde.t1.primal, gcde.t1.fenchel, gcde.t1.donsker), modelname,
                            savepath, r)
    
    
    # generate full comparison table
    parameter.t1 <- simulation.tune(fit.gcde = gcde.t1.primal, nnkcde.h.grid = c(0.5, 0.8, 1.5, 3),
                                    nnkcde.k.grid = c(200, 500, 1000),
                                    nimax = 40)
    
    
    table.t1 <- generate.table(s = 2000, fit.gcde = gcde.t1.primal, 
                               fit.gcde.extra = list(gcde.t1.fenchel, gcde.t1.donsker),
                               otherparameter = parameter.t1, 
                               extra.names = c("GCDE(primal)", "GCDE(fenchel)", "GCDE(donsker)"))
    table.t1
}

clrep <- 10
cl <- makeCluster(clrep)
clusterEvalQ(cl, source("D:/Dropbox/UIowa\ Study\ Materials/Ph.D\ research/GCDE/GCDE.R"))
show.diagnostic.pic <- FALSE
clusterExport(cl, c("show.diagnostic.pic", "savepath", "modelname"),
              envir=environment())
t1.tables <- parLapply(cl, 1:clrep, t1fun)
stopCluster(cl)
save.image(paste0(savepath, modelname, ".RData"))


#process the tables from 10 simulation replications here:
t1.summary.table <- gcde.summary.table(t1.tables)
print(xtable(t1.summary.table, caption = "MSPE Comparison", digits = 3),sanitize.text.function=function(x){x})



#Compasiron of different GCDE fit:
#This is used for tuning the parameters.
check.t1 <- generate.table(s = 2000, fit.gcde = gcde.t1.old, 
                           fit.gcde.extra = list(gcde.t1.old.big, gcde.t1.new, gcde.t1.new.big,
                                                 gcde.t1.dual, gcde.t1.dual.big),
                           extra.names = c("old", "old-3layer", "new", "new-3layer", "dual", "dual-3layer"))
check.t1
print(xtable(check.t1, caption = "MSPE Comparison", digits = 3),sanitize.text.function=function(x){x})
ks.plot(fit.list = list(pri = gcde.t1.primal, fen = gcde.t1.fenchel, don = gcde.t1.donsker), ite.factor = 1500,
        model.name = "T1")



#Plot comparison of all methods
set.seed(2)
plot.x <- c(rnorm(data.t1$dim.x))
plot.cde.comparison(x = plot.x, n.pic = 1, modelname = "T1", 
                    fit.gcde = gcde.t1, otherparameter = parameter.t1)







#######################################################################################################
############################     Model T2           ###################################################
#######################################################################################################
source("D:/Dropbox/UIowa\ Study\ Materials/Ph.D\ research/GCDE/GCDE.R")
show.diagnostic.pic <- FALSE
modelname <- "t2"
savepath <- "D:/Dropbox/UIowa Study Materials/Ph.D research/GCDE/code_simulation_0/trained/"

t2fun <- function(r){
    seed <- 131 * r
    cat("Start the", r, "simulation replicates...")
    data.t2 <- data.generate(n = 5000,
                             x = NULL,
                             dim.x = 5,
                             dim.y = 1,
                             dim.epsilon = 1,
                             epsilon.dist = "normal",
                             model = "add-hete",
                             seed = seed,
                             testdata = TRUE
    )
    
    gcde.t2.primal <- fitgcde(data.t2$x, data.t2$y, data.t2$test.x, data.t2$test.y,
                              diagnostic.info    = data.t2$diagnostic.info,
                              show.diagnostic.pic= show.diagnostic.pic,
                              v.label            = "new",
                              kl.dual            = FALSE,
                              iterations         = 1500,
                              dim.eta            = 3,
                              eta.dist           = "normal",
                              gene.nodes.num     = c(50),
                              dis.nodes.num      = c(50,25),
                              disc.train.type    = "fit",            
                              step.disc          = 10,                  
                              step.gene          = 5,                   
                              display            = 300,               
                              batch_size         = 128,
                              gene.learning.rate = 0.001,
                              disc.learning.rate = 0.0001,
                              gene.lr.decay      = 0,              
                              disc.lr.decay      = 0, 
                              gene.regularizer.1 = regularizer_l2(l = 0.000005),   
                              gene.regularizer.o = regularizer_l2(l = 0.000005),
                              disc.regularizer   = regularizer_l2(l = 0.0000),   
                              activation.disc    = "relu",
                              activation.gene    = "relu",
                              lambda.kl          = 20,
                              lambda.mse         = 0,       
                              m                  = 1, vanillaNN = FALSE,)
    
    gcde.t2.fenchel <- fitgcde(data.t2$x, data.t2$y, data.t2$test.x, data.t2$test.y,
                               diagnostic.info    = data.t2$diagnostic.info,
                               show.diagnostic.pic= show.diagnostic.pic,
                               v.label            = "new",
                               kl.dual            = TRUE,
                               iterations         = 7500,
                               dim.eta            = 3,
                               eta.dist           = "normal",
                               gene.nodes.num     = c(50),
                               dis.nodes.num      = c(50,25),
                               disc.train.type    = "batch",            
                               step.disc          = 5,                  
                               step.gene          = 1,                   
                               display            = 1500,               
                               batch_size         = 128,
                               gene.learning.rate = 0.001,
                               disc.learning.rate = 0.0005,
                               gene.lr.decay      = 0,              
                               disc.lr.decay      = 0, 
                               gene.regularizer.1 = regularizer_l2(l = 0.000005),   
                               gene.regularizer.o = regularizer_l2(l = 0.000005),
                               disc.regularizer   = regularizer_l2(l = 0.0000),   
                               activation.disc    = "relu",
                               activation.gene    = "relu",
                               lambda.kl          = 20,
                               lambda.mse         = 0,       
                               m                  = 1, vanillaNN = FALSE)
    gcde.t2.donsker <- fitgcde(data.t2$x, data.t2$y, data.t2$test.x, data.t2$test.y,
                               diagnostic.info    = data.t2$diagnostic.info,
                               show.diagnostic.pic= show.diagnostic.pic,
                               v.label            = "new",
                               kl.dual            = TRUE,
                               dual.form          = "donsker",
                               iterations         = 7500,
                               dim.eta            = 3,
                               eta.dist           = "normal",
                               gene.nodes.num     = c(50),
                               dis.nodes.num      = c(50,25),
                               disc.train.type    = "batch",            
                               step.disc          = 5,                  
                               step.gene          = 1,                   
                               display            = 1500,               
                               batch_size         = 128,
                               gene.learning.rate = 0.001,
                               disc.learning.rate = 0.0005,
                               gene.lr.decay      = 0,              
                               disc.lr.decay      = 0, 
                               gene.regularizer.1 = regularizer_l2(l = 0.000005),   
                               gene.regularizer.o = regularizer_l2(l = 0.000005),
                               disc.regularizer   = regularizer_l2(l = 0.0000),   
                               activation.disc    = "relu",
                               activation.gene    = "relu",
                               lambda.kl          = 20,
                               lambda.mse         = 0,       
                               m                  = 1, vanillaNN = FALSE)
    
    gcde.save.model.weights(list(gcde.t2.primal, gcde.t2.fenchel, gcde.t2.donsker), modelname,
                            savepath, r)
    
    # generate full comparison table
    parameter.t2 <- simulation.tune(fit.gcde = gcde.t2.primal, nnkcde.h.grid = c(0.8, 1.5, 3 ,5 ),
                                    nnkcde.k.grid = c(200, 500, 1000),
                                    nimax = 40)
    
    table.t2 <- generate.table(s = 2000, fit.gcde = gcde.t2.primal, 
                               fit.gcde.extra = list(gcde.t2.fenchel, gcde.t2.donsker),
                               otherparameter = parameter.t2, 
                               extra.names = c("GCDE(primal)", "GCDE(fenchel)", "GCDE(donsker)"))
    table.t2
}

clrep <- 10
cl <- makeCluster(clrep)
clusterEvalQ(cl, source("D:/Dropbox/UIowa\ Study\ Materials/Ph.D\ research/GCDE/GCDE.R"))
show.diagnostic.pic <- FALSE
clusterExport(cl, c("show.diagnostic.pic", "savepath", "modelname"),
              envir=environment())
t2.tables <- parLapply(cl, 1:clrep, t2fun)
stopCluster(cl)
save.image(paste0(savepath, modelname, ".RData"))



#process the tables from 10 simulation replications here:
t2.summary.table <- gcde.summary.table(t2.tables)
print(xtable(t2.summary.table, caption = "MSPE Comparison", digits = 3),sanitize.text.function=function(x){x})




#Compare different GCDE fit 
check.t2 <- generate.table(s = 2000, fit.gcde = gcde.t2.primal, 
                           fit.gcde.extra = list(gcde.t2.fenchel, gcde.t2.donsker),
                           extra.names = c("pri","fenchel","donsker"))
check.t2
print(xtable(check.t2, caption = "MSPE Comparison", digits = 3),sanitize.text.function=function(x){x})
ks.plot(fit.list = list(pri = gcde.t2.primal, fen = gcde.t2.fenchel, don = gcde.t2.donsker), ite.factor = 1500,
        model.name = "T2")


#Comparison of all methods
set.seed(2)
plot.x <- c(rnorm(data.t2$dim.x))
plot.cde.comparison(x = plot.x, n.pic = 1, modelname = "T2", 
                    fit.gcde = gcde.t2, otherparameter = parameter.t2)


#######################################################################################################
############################     Model T3 (previously labeled as T5)           ########################
#######################################################################################################
source("D:/Dropbox/UIowa\ Study\ Materials/Ph.D\ research/GCDE/GCDE.R")
show.diagnostic.pic <- FALSE
modelname <- "t3"
savepath <- "D:/Dropbox/UIowa Study Materials/Ph.D research/GCDE/code_simulation_0/trained/"


t3fun <- function(r){
    seed <- 131 * r
    cat("Start the", r, "simulation replicates...")
    data.t3 <- data.generate(n = 5000,
                             x = NULL,
                             dim.x = 30,
                             dim.y = 1,
                             dim.epsilon = 1,
                             epsilon.dist = "mixture",
                             model = "mul-simple-moved",
                             seed = seed,
                             testdata = TRUE
    )
    
    gcde.t3.primal <- fitgcde(data.t3$x, data.t3$y, data.t3$test.x, data.t3$test.y,
                              diagnostic.info    = data.t3$diagnostic.info,
                              show.diagnostic.pic= show.diagnostic.pic,
                              v.label            = "new",
                              kl.dual            = FALSE,
                              iterations         = 1500,
                              dim.eta            = 3,
                              eta.dist           = "normal",
                              gene.nodes.num     = c(50),
                              dis.nodes.num      = c(50,25),
                              disc.train.type    = "fit",            
                              step.disc          = 10,                
                              step.gene          = 5,                 
                              display            = 300,             
                              batch_size         = 128,
                              gene.learning.rate = 0.001,
                              disc.learning.rate = 0.0001,
                              gene.lr.decay      = 0,              
                              disc.lr.decay      = 0, 
                              gene.regularizer.1 = regularizer_l2(l = 0.0005),   
                              gene.regularizer.o = regularizer_l2(l = 0.0005),
                              disc.regularizer   = regularizer_l2(l = 0.0000),   
                              activation.disc    = "relu",
                              activation.gene    = "relu",
                              lambda.kl          = 20,
                              lambda.mse         = 0,       
                              m                  = 1, vanillaNN = FALSE)
    
    gcde.t3.fenchel <- fitgcde(data.t3$x, data.t3$y, data.t3$test.x, data.t3$test.y,
                            diagnostic.info    = data.t3$diagnostic.info,
                            show.diagnostic.pic= show.diagnostic.pic,
                            v.label            = "new",
                            kl.dual            = TRUE,
                            iterations         = 7500,
                            dim.eta            = 3,
                            eta.dist           = "normal",
                            gene.nodes.num     = c(50),
                            dis.nodes.num      = c(50,25),
                            disc.train.type    = "batch",            
                            step.disc          = 5,                
                            step.gene          = 1,                 
                            display            = 1500,             
                            batch_size         = 128,
                            gene.learning.rate = 0.001,
                            disc.learning.rate = 0.0001,
                            gene.lr.decay      = 0,              
                            disc.lr.decay      = 0, 
                            gene.regularizer.1 = regularizer_l2(l = 0.0005),   
                            gene.regularizer.o = regularizer_l2(l = 0.0005),
                            disc.regularizer   = regularizer_l2(l = 0.0000),   
                            activation.disc    = "relu",
                            activation.gene    = "relu",
                            lambda.kl          = 20,
                            lambda.mse         = 0,       
                            m                  = 1, vanillaNN = FALSE)
    gcde.t3.donsker <- fitgcde(data.t3$x, data.t3$y, data.t3$test.x, data.t3$test.y,
                               diagnostic.info    = data.t3$diagnostic.info,
                               show.diagnostic.pic= show.diagnostic.pic,
                               v.label            = "new",
                               kl.dual            = TRUE,
                               dual.form          = "donsker",
                               iterations         = 7500,
                               dim.eta            = 3,
                               eta.dist           = "normal",
                               gene.nodes.num     = c(50),
                               dis.nodes.num      = c(50,25),
                               disc.train.type    = "batch",            
                               step.disc          = 5,                
                               step.gene          = 1,                 
                               display            = 1500,             
                               batch_size         = 128,
                               gene.learning.rate = 0.001,
                               disc.learning.rate = 0.0001,
                               gene.lr.decay      = 0,              
                               disc.lr.decay      = 0, 
                               gene.regularizer.1 = regularizer_l2(l = 0.0005),   
                               gene.regularizer.o = regularizer_l2(l = 0.0005),
                               disc.regularizer   = regularizer_l2(l = 0.0000),   
                               activation.disc    = "relu",
                               activation.gene    = "relu",
                               lambda.kl          = 20,
                               lambda.mse         = 0,       
                               m                  = 1, vanillaNN = FALSE)
    
    gcde.save.model.weights(list(gcde.t3.primal, gcde.t3.fenchel, gcde.t3.donsker), modelname,
                            savepath, r)
    
    # generate full comparison table
    parameter.t3 <- simulation.tune(fit.gcde = gcde.t3.primal, nnkcde.h.grid = c( 5 ,10, 30, 40 ),
                                    nnkcde.k.grid = c(200, 500, 1000),
                                    nimax = 40)
    
    table.t3 <- generate.table(s = 2000, fit.gcde = gcde.t3.primal,
                               fit.gcde.extra = list(gcde.t3.fenchel, gcde.t3.donsker),
                               otherparameter = parameter.t3, 
                               extra.names = c("GCDE(primal)", "GCDE(fenchel)", "GCDE(donsker)"))
    table.t3
}

clrep <- 10
cl <- makeCluster(clrep)
clusterEvalQ(cl, source("D:/Dropbox/UIowa\ Study\ Materials/Ph.D\ research/GCDE/GCDE.R"))
show.diagnostic.pic <- FALSE
clusterExport(cl, c("show.diagnostic.pic", "savepath", "modelname"),
              envir=environment())
t3.tables <- parLapply(cl, 1:clrep, t3fun)
stopCluster(cl)
save.image(paste0(savepath, modelname, ".RData"))


#process the tables from 10 simulation replications here:
t3.summary.table <- gcde.summary.table(t3.tables)
print(xtable(t3.summary.table, caption = "MSPE Comparison", digits = 3),sanitize.text.function=function(x){x})


#Compare different GCDE fit 
check.t3 <- generate.table(s = 2000, fit.gcde = gcde.t3.primal, 
                           fit.gcde.extra = list(gcde.t3.fenchel, gcde.t3.donsker),
                           extra.names = c("pri","fenchel","donsker"))
check.t3
print(xtable(check.t3, caption = "MSPE Comparison", digits = 3),sanitize.text.function=function(x){x})
ks.plot(fit.list = list(pri = gcde.t3.primal, fen = gcde.t3.fenchel, don = gcde.t3.donsker), ite.factor = 1500,
        model.name = "T3")




set.seed(2)
plot.x <- c(1, rnorm(29))
plot.cde.comparison(x = plot.x, n.pic = 1, modelname = "T3", 
                    fit.gcde = gcde.t3, otherparameter = parameter.t3)





#######################################################################################################
############################     Model T4(previously labeled as T6)            ########################
#######################################################################################################
source("D:/Dropbox/UIowa\ Study\ Materials/Ph.D\ research/GCDE/GCDE.R")
show.diagnostic.pic <- FALSE
modelname <- "t4"
savepath <- "D:/Dropbox/UIowa Study Materials/Ph.D research/GCDE/code_simulation_0/trained/"


t4fun <- function(r){
    seed <- 131 * r
    cat("Start the", r, "simulation replicates...")
    data.t4 <- data.generate(n = 5000,
                             x = NULL,
                             dim.x = 5,
                             dim.y = 1,
                             dim.epsilon = 1,
                             epsilon.dist = "normal",
                             model = "add-toy",
                             seed = seed,
                             testdata = TRUE
    )
    
    
    gcde.t4.primal <- fitgcde(data.t4$x, data.t4$y, data.t4$test.x, data.t4$test.y,
                              diagnostic.info    = data.t4$diagnostic.info,
                              show.diagnostic.pic= show.diagnostic.pic,
                              v.label            = "new",
                              kl.dual            = FALSE,
                              iterations         = 1500,
                              dim.eta            = 4,
                              eta.dist           = "normal",
                              gene.nodes.num     = c(40,15),
                              dis.nodes.num      = c(50,25),
                              disc.train.type    = "fit",            
                              step.disc          = 10,                  
                              step.gene          = 5,                  
                              display            = 300,               
                              batch_size         = 128,
                              gene.learning.rate = 0.001,
                              disc.learning.rate = 0.0001,
                              gene.lr.decay      = 0,              
                              disc.lr.decay      = 0, 
                              gene.regularizer.1 = regularizer_l2(l = 0.00001),   
                              gene.regularizer.o = regularizer_l2(l = 0.0000),
                              disc.regularizer   = regularizer_l2(l = 0.0000),  
                              activation.disc    = "relu",
                              activation.gene    = "relu",
                              lambda.kl          = 20,
                              lambda.mse         = 0,       
                              m                  = 1, vanillaNN = FALSE)
    
    gcde.t4.fenchel <- fitgcde(data.t4$x, data.t4$y, data.t4$test.x, data.t4$test.y,
                            diagnostic.info    = data.t4$diagnostic.info,
                            show.diagnostic.pic= show.diagnostic.pic,
                            v.label            = "new",
                            kl.dual            = TRUE,
                            iterations         = 7500,
                            dim.eta            = 4,
                            eta.dist           = "normal",
                            gene.nodes.num     = c(40,15),
                            dis.nodes.num      = c(50,25),
                            disc.train.type    = "batch",            
                            step.disc          = 5,                  
                            step.gene          = 1,                  
                            display            = 1500,               
                            batch_size         = 128,
                            gene.learning.rate = 0.001,
                            disc.learning.rate = 0.0001,
                            gene.lr.decay      = 0,              
                            disc.lr.decay      = 0, 
                            gene.regularizer.1 = regularizer_l2(l = 0.00001),   
                            gene.regularizer.o = regularizer_l2(l = 0.00000),
                            disc.regularizer   = regularizer_l2(l = 0.00000),  
                            activation.disc    = "relu",
                            activation.gene    = "relu",
                            lambda.kl          = 20,
                            lambda.mse         = 0,       
                            m                  = 1, vanillaNN = FALSE)
    gcde.t4.donsker <- fitgcde(data.t4$x, data.t4$y, data.t4$test.x, data.t4$test.y,
                            diagnostic.info    = data.t4$diagnostic.info,
                            show.diagnostic.pic= show.diagnostic.pic,
                            v.label            = "new",
                            kl.dual            = TRUE,
                            dual.form          = "donsker",
                            iterations         = 7500,
                            dim.eta            = 4,
                            eta.dist           = "normal",
                            gene.nodes.num     = c(40,15),
                            dis.nodes.num      = c(50,25),
                            disc.train.type    = "batch",            
                            step.disc          = 5,                  
                            step.gene          = 1,                  
                            display            = 1500,               
                            batch_size         = 128,
                            gene.learning.rate = 0.001,
                            disc.learning.rate = 0.0001,
                            gene.lr.decay      = 0,              
                            disc.lr.decay      = 0, 
                            gene.regularizer.1 = regularizer_l2(l = 0.00001),   
                            gene.regularizer.o = regularizer_l2(l = 0.00000),
                            disc.regularizer   = regularizer_l2(l = 0.00000),  
                            activation.disc    = "relu",
                            activation.gene    = "relu",
                            lambda.kl          = 20,
                            lambda.mse         = 0,       
                            m                  = 1, vanillaNN = FALSE)
    gcde.save.model.weights(list(gcde.t4.primal, gcde.t4.fenchel, gcde.t4.donsker), modelname,
                            savepath, r)
    
    # generate full comparison table
    parameter.t4 <- simulation.tune(fit.gcde = gcde.t4.primal, nnkcde.h.grid = c( 0.5, 1.5, 2, 3),
                                    nnkcde.k.grid = c(200, 300, 500, 1000),
                                    nimax = 40)
    
    table.t4 <- generate.table(s = 2000, fit.gcde = gcde.t4.primal, 
                               fit.gcde.extra = list(gcde.t4.fenchel, gcde.t4.donsker),
                               otherparameter = parameter.t4, 
                               extra.names = c("GCDE(primal)", "GCDE(fenchel)", "GCDE(donsker)"))
    table.t4
}

clrep <- 10
cl <- makeCluster(clrep)
clusterEvalQ(cl, source("D:/Dropbox/UIowa\ Study\ Materials/Ph.D\ research/GCDE/GCDE.R"))
show.diagnostic.pic <- FALSE
clusterExport(cl, c("show.diagnostic.pic", "savepath", "modelname"),
              envir=environment())
t4.tables <- parLapply(cl, 1:clrep, t4fun)
stopCluster(cl)
save.image(paste0(savepath, modelname, ".RData"))


#process the tables from 10 simulation replications here:
t4.summary.table <- gcde.summary.table(t4.tables)
print(xtable(t4.summary.table, caption = "MSPE Comparison", digits = 3),sanitize.text.function=function(x){x})






#check perfoamance for GCDE only
check.t4 <- generate.table(s = 2000, fit.gcde = gcde.t4.primal, 
                           fit.gcde.extra = list(gcde.t4.fenchel, gcde.t4.donsker),
                           extra.names = c("primal","Fenchel","donsker"))
check.t4
knitr::kable(check.t4, digits = 3)
print(xtable(check.t4, caption = "MSPE Comparison", digits = 3),sanitize.text.function=function(x){x})
ks.plot(fit.list = list(primal = gcde.t4.primal, fenc = gcde.t4.fenchel, don = gcde.t4.donsker), 
        ite.factor = 1500,
        model.name = "t4")



#Whole comparison
set.seed(12)
plot.x <- c(rnorm(5))
plot.cde.comparison(x = plot.x, n.pic = 1, modelname = "t4", 
                    fit.gcde = gcde.t4, otherparameter = parameter.t4)








#######################################################################################################
############################     Model T5(previously labeled as T3)          ##########################
#######################################################################################################
source("D:/Dropbox/UIowa\ Study\ Materials/Ph.D\ research/GCDE/GCDE.R")
show.diagnostic.pic <- FALSE
modelname <- "t5"
savepath <- "D:/Dropbox/UIowa Study Materials/Ph.D research/GCDE/code_simulation_0/trained/"


t5fun <- function(r){
    seed <- 131 * r
    cat("Start the", r, "simulation replicates...")
    data.t5 <- data.generate(n = 5000,
                             x = NULL,
                             dim.x = 5,
                             dim.y = 1,
                             dim.epsilon = 1,
                             epsilon.dist = "normal",
                             model = "mul-simple-moved",
                             seed = seed,
                             testdata = TRUE
    )
    
    gcde.t5.primal <- fitgcde(data.t5$x, data.t5$y, data.t5$test.x, data.t5$test.y,
                       diagnostic.info    = data.t5$diagnostic.info,
                       show.diagnostic.pic= show.diagnostic.pic,
                       v.label            = "new",
                       kl.dual            = FALSE,
                       iterations         = 1500,
                       dim.eta            = 3,
                       eta.dist           = "normal",
                       gene.nodes.num     = c(50),
                       dis.nodes.num      = c(50,25),
                       disc.train.type    = "fit",            
                       step.disc          = 10,                  
                       step.gene          = 5,                   
                       display            = 300,               
                       batch_size         = 128,
                       gene.learning.rate = 0.001,
                       disc.learning.rate = 0.0001,
                       gene.lr.decay      = 0,  #No decay is better than 1/500.            
                       disc.lr.decay      = 0, 
                       gene.regularizer.1 = regularizer_l2(l = 0.000005),   
                       gene.regularizer.o = regularizer_l2(l = 0.000005),
                       disc.regularizer   = regularizer_l2(l = 0.0000),   
                       activation.disc    = "relu",
                       activation.gene    = "relu",
                       lambda.kl          = 20,
                       lambda.mse         = 0,       
                       m                  = 1, vanillaNN = FALSE)
    
    gcde.t5.fenchel <- fitgcde(data.t5$x, data.t5$y, data.t5$test.x, data.t5$test.y,
                            diagnostic.info    = data.t5$diagnostic.info,
                            show.diagnostic.pic= show.diagnostic.pic,
                            v.label            = "new",
                            kl.dual            = TRUE,
                            iterations         = 7500,
                            dim.eta            = 3,
                            eta.dist           = "normal",
                            gene.nodes.num     = c(50),
                            dis.nodes.num      = c(50,25),
                            disc.train.type    = "batch",            
                            step.disc          = 5,                  
                            step.gene          = 1,                   
                            display            = 1500,               
                            batch_size         = 128,
                            gene.learning.rate = 0.001,
                            disc.learning.rate = 0.0005,
                            gene.lr.decay      = 0,              
                            disc.lr.decay      = 0, 
                            gene.regularizer.1 = regularizer_l2(l = 0.000005),   
                            gene.regularizer.o = regularizer_l2(l = 0.000005),
                            disc.regularizer   = regularizer_l2(l = 0.0000),   
                            activation.disc    = "relu",
                            activation.gene    = "relu",
                            lambda.kl          = 20,
                            lambda.mse         = 0,       
                            m                  = 1, vanillaNN = FALSE)
    gcde.t5.donsker <- fitgcde(data.t5$x, data.t5$y, data.t5$test.x, data.t5$test.y,
                            diagnostic.info    = data.t5$diagnostic.info,
                            show.diagnostic.pic= show.diagnostic.pic,
                            v.label            = "new",
                            kl.dual            = TRUE,
                            dual.form          = "donsker",
                            iterations         = 7500,
                            dim.eta            = 3,
                            eta.dist           = "normal",
                            gene.nodes.num     = c(50),
                            dis.nodes.num      = c(50,25),
                            disc.train.type    = "batch",            
                            step.disc          = 5,                  
                            step.gene          = 1,                   
                            display            = 1500,               
                            batch_size         = 128,
                            gene.learning.rate = 0.001,
                            disc.learning.rate = 0.0005,
                            gene.lr.decay      = 0,              
                            disc.lr.decay      = 0, 
                            gene.regularizer.1 = regularizer_l2(l = 0.000005),   
                            gene.regularizer.o = regularizer_l2(l = 0.000005),
                            disc.regularizer   = regularizer_l2(l = 0.0000),   
                            activation.disc    = "relu",
                            activation.gene    = "relu",
                            lambda.kl          = 20,
                            lambda.mse         = 0,       
                            m                  = 1, vanillaNN = FALSE)
    
    gcde.save.model.weights(list(gcde.t5.primal, gcde.t5.fenchel, gcde.t5.donsker), modelname,
                            savepath, r)
    
    # generate full comparison table
    parameter.t5 <- simulation.tune(fit.gcde = gcde.t5.primal, nnkcde.h.grid = c(1.5, 3 ,5 ,10,15),
                                    nnkcde.k.grid = c(200, 500, 1000),
                                    nimax = 40)
    
    table.t5 <- generate.table(s = 2000, fit.gcde = gcde.t5.primal, 
                               fit.gcde.extra = list(gcde.t5.fenchel, gcde.t5.donsker),
                               otherparameter = parameter.t5, 
                               extra.names = c("GCDE(primal)", "GCDE(fenchel)", "GCDE(donsker)"))
    table.t5
}


clrep <- 10
cl <- makeCluster(clrep)
clusterEvalQ(cl, source("D:/Dropbox/UIowa\ Study\ Materials/Ph.D\ research/GCDE/GCDE.R"))
show.diagnostic.pic <- FALSE
clusterExport(cl, c("show.diagnostic.pic", "savepath", "modelname"),
              envir=environment())
t5.tables <- parLapply(cl, 1:clrep, t5fun)
stopCluster(cl)
save.image(paste0(savepath, modelname, ".RData"))


#process the tables from 10 simulation replications here:
t5.summary.table <- gcde.summary.table(t5.tables)
print(xtable(t5.summary.table, caption = "MSPE Comparison", digits = 3),sanitize.text.function=function(x){x})



#Plot comparison
set.seed(2)
plot.x <- c(rnorm(data.t5$dim.x))
plot.cde.comparison(x = plot.x, n.pic = 1, modelname = "t5", 
                    fit.gcde = gcde.t5, otherparameter = parameter.t5)



#######################################################################################################
############################     Model T6(previously labeled as T4)         ###########################
#######################################################################################################
source("D:/Dropbox/UIowa\ Study\ Materials/Ph.D\ research/GCDE/GCDE.R")
show.diagnostic.pic <- FALSE
modelname <- "t6"
savepath <- "D:/Dropbox/UIowa Study Materials/Ph.D research/GCDE/code_simulation_0/trained/"


t6fun <- function(r){
    seed <- 131 * r
    cat("Start the", r, "simulation replicates...")
    data.t6 <- data.generate(n = 5000,
                             x = NULL,
                             dim.x = 5,
                             dim.y = 1,
                             dim.epsilon = 1,
                             epsilon.dist = "mixture",
                             model = "mul-simple-moved",
                             seed = seed,
                             testdata = TRUE
    )
    
    gcde.t6.primal <- fitgcde(data.t6$x, data.t6$y, data.t6$test.x, data.t6$test.y,
                           diagnostic.info    = data.t6$diagnostic.info,
                           show.diagnostic.pic= show.diagnostic.pic,
                           v.label            = "new",
                           kl.dual            = FALSE,
                           iterations         = 1500,
                           dim.eta            = 3,
                           eta.dist           = "normal",
                           gene.nodes.num     = c(50),
                           dis.nodes.num      = c(50,25),
                           disc.train.type    = "fit",            
                           step.disc          = 10,                  
                           step.gene          = 5,                   
                           display            = 300,               
                           batch_size         = 128,
                           gene.learning.rate = 0.001,
                           disc.learning.rate = 0.0001,
                           gene.lr.decay      = 0,              
                           disc.lr.decay      = 0, 
                           gene.regularizer.1 = regularizer_l2(l = 0.00001),   
                           gene.regularizer.o = regularizer_l2(l = 0.00001),
                           disc.regularizer   = regularizer_l2(l = 0.0000),   
                           activation.disc    = "relu",
                           activation.gene    = "relu",
                           lambda.kl          = 20,
                           lambda.mse         = 0,       
                           m                  = 1, vanillaNN = FALSE)
    
    gcde.t6.fenchel <- fitgcde(data.t6$x, data.t6$y, data.t6$test.x, data.t6$test.y,
                            diagnostic.info    = data.t6$diagnostic.info,
                            show.diagnostic.pic= show.diagnostic.pic,
                            v.label            = "new",
                            kl.dual            = TRUE,
                            iterations         = 7500,
                            dim.eta            = 3,
                            eta.dist           = "normal",
                            gene.nodes.num     = c(50),
                            dis.nodes.num      = c(50,25),
                            disc.train.type    = "batch",            
                            step.disc          = 4,                  
                            step.gene          = 1,                   
                            display            = 1500,               
                            batch_size         = 128,
                            gene.learning.rate = 0.001,
                            disc.learning.rate = 0.0005,
                            gene.lr.decay      = 0,              
                            disc.lr.decay      = 0, 
                            gene.regularizer.1 = regularizer_l2(l = 0.00001),   
                            gene.regularizer.o = regularizer_l2(l = 0.00001),
                            disc.regularizer   = regularizer_l2(l = 0.0000),   
                            activation.disc    = "relu",
                            activation.gene    = "relu",
                            lambda.kl          = 20,
                            lambda.mse         = 0,       
                            m                  = 1, vanillaNN = FALSE)
    gcde.t6.donsker <- fitgcde(data.t6$x, data.t6$y, data.t6$test.x, data.t6$test.y,
                            diagnostic.info    = data.t6$diagnostic.info,
                            show.diagnostic.pic= show.diagnostic.pic,
                            v.label            = "new",
                            kl.dual            = TRUE,
                            dual.form          = "donsker",
                            iterations         = 7500,
                            dim.eta            = 3,
                            eta.dist           = "normal",
                            gene.nodes.num     = c(50),
                            dis.nodes.num      = c(50,25),
                            disc.train.type    = "batch",            
                            step.disc          = 4,                  
                            step.gene          = 1,                   
                            display            = 1500,               
                            batch_size         = 128,
                            gene.learning.rate = 0.001,
                            disc.learning.rate = 0.0005,
                            gene.lr.decay      = 0,              
                            disc.lr.decay      = 0, 
                            gene.regularizer.1 = regularizer_l2(l = 0.00001),   
                            gene.regularizer.o = regularizer_l2(l = 0.00001),
                            disc.regularizer   = regularizer_l2(l = 0.0000),   
                            activation.disc    = "relu",
                            activation.gene    = "relu",
                            lambda.kl          = 20,
                            lambda.mse         = 0,       
                            m                  = 1, vanillaNN = FALSE)
    
    gcde.save.model.weights(list(gcde.t6.primal, gcde.t6.fenchel, gcde.t6.donsker), modelname,
                            savepath, r)
    # generate full comparison table
    parameter.t6 <- simulation.tune(fit.gcde = gcde.t6.primal, nnkcde.h.grid = c(5 ,10, 30),
                                    nnkcde.k.grid = c(200, 500, 1000),
                                    nimax = 40)
    
    table.t6 <- generate.table(s = 2000, fit.gcde = gcde.t6.primal, 
                               fit.gcde.extra = list(gcde.t6.fenchel, gcde.t6.donsker),
                               otherparameter = parameter.t6, 
                               extra.names = c("GCDE(primal)", "GCDE(fenchel)", "GCDE(donsker)"))
    table.t6
}

clrep <- 10
cl <- makeCluster(clrep)
clusterEvalQ(cl, source("D:/Dropbox/UIowa\ Study\ Materials/Ph.D\ research/GCDE/GCDE.R"))
show.diagnostic.pic <- FALSE
clusterExport(cl, c("show.diagnostic.pic", "savepath", "modelname"),
              envir=environment())
t6.tables <- parLapply(cl, 1:clrep, t6fun)
stopCluster(cl)
save.image(paste0(savepath, modelname, ".RData"))



#process the tables from 10 simulation replications here:
t6.summary.table <- gcde.summary.table(t6.tables)
print(xtable(t6.summary.table, caption = "MSPE Comparison", digits = 3),sanitize.text.function=function(x){x})


#check perfoamance for GCDE only
check.t6 <- generate.table(s = 2000, fit.gcde = gcde.t6.new, 
                           fit.gcde.extra = list(gcde.t6.new.leaky, gcde.t6.dual, gcde.t6.dual.leaky),
                           extra.names = c("primal", "primal leaky", "dual", "dual leaky"))
check.t6
print(xtable(check.t6, caption = "MSPE Comparison", digits = 3),sanitize.text.function=function(x){x})
ks.plot(fit.list = list(primal = gcde.t6.new, primalleaky = gcde.t6.new.leaky, 
                        dual = gcde.t6.dual, dualleaky = gcde.t6.dual.leaky), ite.factor = 1500,
        model.name = "t6")



#Generate comparison table
set.seed(2)
plot.x <- c(rnorm(data.t6$dim.x))
plot.cde.comparison(x = plot.x, n.pic = 1, modelname = "t6", 
                    fit.gcde = gcde.t6, otherparameter = parameter.t6)








