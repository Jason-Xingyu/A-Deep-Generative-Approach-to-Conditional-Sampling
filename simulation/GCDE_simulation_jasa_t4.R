my_code_path <- "../../GCDE/code"
sapply(list.files(my_code_path, full.names = TRUE), source)


# model t4 
##############################################################################
####                          Conpare different n
#############################################################################
show.diagnostic.pic <- FALSE
modelname <- "t4"
savepath <- "./trained_jasa/"

t4fun <- function(r){
    seed <- 131 * r
    cat("Start the", r, "simulation replicates...")
    exp.name <- c("n1000", "n2500", "n5000", "n7500", "n10000")
    data.t4.1000 <- data.generate(n = 1000,
                                  x = NULL,
                                  dim.x = 5,
                                  dim.y = 1,
                                  dim.epsilon = 1,
                                  epsilon.dist = "normal",
                                  model = "add-toy",
                                  seed = seed,
                                  testdata = TRUE
    )
    
    data.t4.2500 <- data.generate(n = 2500,
                                  x = NULL,
                                  dim.x = 5,
                                  dim.y = 1,
                                  dim.epsilon = 1,
                                  epsilon.dist = "normal",
                                  model = "add-toy",
                                  seed = seed,
                                  testdata = TRUE
    )
    
    data.t4.5000 <- data.generate(n = 5000,
                                  x = NULL,
                                  dim.x = 5,
                                  dim.y = 1,
                                  dim.epsilon = 1,
                                  epsilon.dist = "normal",
                                  model = "add-toy",
                                  seed = seed,
                                  testdata = TRUE
    )
    
    data.t4.7500 <- data.generate(n = 7500,
                                  x = NULL,
                                  dim.x = 5,
                                  dim.y = 1,
                                  dim.epsilon = 1,
                                  epsilon.dist = "normal",
                                  model = "add-toy",
                                  seed = seed,
                                  testdata = TRUE
    )
    
    data.t4.10000 <- data.generate(n = 10000,
                                   x = NULL,
                                   dim.x = 5,
                                   dim.y = 1,
                                   dim.epsilon = 1,
                                   epsilon.dist = "normal",
                                   model = "add-toy",
                                   seed = seed,
                                   testdata = TRUE
    )
    
    gcde.t4.n1000 <- fitgcde(data.t4.1000$x, data.t4.1000$y, data.t4.1000$test.x, data.t4.1000$test.y,
                             #load.disc = paste0(savepath,modelname,"/",modelname,"_",exp.name[1],
                             #                   "_discriminator_", r, ".h5"),
                             #load.gene = paste0(savepath,modelname,"/",modelname,"_",exp.name[1],
                             #                   "_generator_", r, ".h5"),
                             diagnostic.info    = data.t4.1000$diagnostic.info,
                             show.diagnostic.pic= show.diagnostic.pic,
                             v.label            = "new",
                             kl.dual            = TRUE,
                             dual.form          = "fenchel",
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
                             gene.regularizer.o = regularizer_l2(l = 0.0000),
                             disc.regularizer   = regularizer_l2(l = 0.0000),   
                             activation.disc    = "relu",
                             activation.gene    = "relu",
                             lambda.kl          = 20,
                             lambda.mse         = 0,       
                             m                  = 1, vanillaNN = FALSE)
    
    gcde.t4.n2500 <- fitgcde(data.t4.2500$x, data.t4.2500$y, data.t4.2500$test.x, data.t4.2500$test.y,
                             #load.disc = paste0(savepath,modelname,"/",modelname,"_",exp.name[2],
                             #                   "_discriminator_", r, ".h5"),
                             #load.gene = paste0(savepath,modelname,"/",modelname,"_",exp.name[2],
                             #                   "_generator_", r, ".h5"),
                             diagnostic.info    = data.t4.2500$diagnostic.info,
                             show.diagnostic.pic= show.diagnostic.pic,
                             v.label            = "new",
                             kl.dual            = TRUE,
                             dual.form          = "fenchel",
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
                             gene.regularizer.o = regularizer_l2(l = 0.0000),
                             disc.regularizer   = regularizer_l2(l = 0.0000),   
                             activation.disc    = "relu",
                             activation.gene    = "relu",
                             lambda.kl          = 20,
                             lambda.mse         = 0,       
                             m                  = 1, vanillaNN = FALSE)
    
    gcde.t4.n5000 <- fitgcde(data.t4.5000$x, data.t4.5000$y, data.t4.5000$test.x, data.t4.5000$test.y,
                             #load.disc = paste0(savepath,modelname,"/",modelname,"_",exp.name[3],
                             #                   "_discriminator_", r, ".h5"),
                             #load.gene = paste0(savepath,modelname,"/",modelname,"_",exp.name[3],
                             #                   "_generator_", r, ".h5"),
                             diagnostic.info    = data.t4.5000$diagnostic.info,
                             show.diagnostic.pic= show.diagnostic.pic,
                             v.label            = "new",
                             kl.dual            = TRUE,
                             dual.form          = "fenchel",
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
                             gene.regularizer.o = regularizer_l2(l = 0.0000),
                             disc.regularizer   = regularizer_l2(l = 0.0000),   
                             activation.disc    = "relu",
                             activation.gene    = "relu",
                             lambda.kl          = 20,
                             lambda.mse         = 0,       
                             m                  = 1, vanillaNN = FALSE)
    
    gcde.t4.n7500 <- fitgcde(data.t4.7500$x, data.t4.7500$y, data.t4.7500$test.x, data.t4.7500$test.y,
                             #load.disc = paste0(savepath,modelname,"/",modelname,"_",exp.name[4],
                             #                   "_discriminator_", r, ".h5"),
                             #load.gene = paste0(savepath,modelname,"/",modelname,"_",exp.name[4],
                             #                   "_generator_", r, ".h5"),
                             diagnostic.info    = data.t4.7500$diagnostic.info,
                             show.diagnostic.pic= show.diagnostic.pic,
                             v.label            = "new",
                             kl.dual            = TRUE,
                             dual.form          = "fenchel",
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
                             gene.regularizer.o = regularizer_l2(l = 0.0000),
                             disc.regularizer   = regularizer_l2(l = 0.0000),   
                             activation.disc    = "relu",
                             activation.gene    = "relu",
                             lambda.kl          = 20,
                             lambda.mse         = 0,       
                             m                  = 1, vanillaNN = FALSE)
    
    gcde.t4.n10000 <- fitgcde(data.t4.10000$x, data.t4.10000$y, data.t4.10000$test.x, data.t4.10000$test.y,
                              #load.disc = paste0(savepath,modelname,"/",modelname,"_",exp.name[5],
                              #                   "_discriminator_", r, ".h5"),
                              #load.gene = paste0(savepath,modelname,"/",modelname,"_",exp.name[5],
                              #                   "_generator_", r, ".h5"),
                              diagnostic.info    = data.t4.10000$diagnostic.info,
                              show.diagnostic.pic= show.diagnostic.pic,
                              v.label            = "new",
                              kl.dual            = TRUE,
                              dual.form          = "fenchel",
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
                              gene.regularizer.o = regularizer_l2(l = 0.0000),
                              disc.regularizer   = regularizer_l2(l = 0.0000),   
                              activation.disc    = "relu",
                              activation.gene    = "relu",
                              lambda.kl          = 20,
                              lambda.mse         = 0,       
                              m                  = 1, vanillaNN = FALSE)
    
    
    gcde.save.model.weights.jasa(list(gcde.t4.n1000, gcde.t4.n2500, gcde.t4.n5000, gcde.t4.n7500,
                                      gcde.t4.n10000), modelname,exp.name,
                                 savepath, r)
    
    
    table.t4 <- generate.table(s = 2000, fit.gcde = gcde.t4.n1000, 
                               fit.gcde.extra = list(gcde.t4.n2500,gcde.t4.n5000, gcde.t4.n7500,
                                                     gcde.t4.n10000),
                               otherparameter = NULL, 
                               extra.names = exp.name)
    table.t4
}

clrep <- 10
cl <- makeCluster(clrep)
show.diagnostic.pic <- FALSE
clusterExport(cl, c("show.diagnostic.pic", "savepath", "modelname","my_code_path"),
              envir=environment())
clusterEvalQ(cl, sapply(list.files(my_code_path, full.names = TRUE), source))

t4.tables <- parLapply(cl, 1:clrep, t4fun)
stopCluster(cl)
save.image(paste0(savepath, modelname, "_n.RData"))

modelname <- "t4"
load(paste0(savepath, modelname, "_n.RData"))


#process the tables from 10 simulation replications here:
t4.summary.table <- gcde.summary.table(t4.tables)
print(xtable(t4.summary.table, caption = "MSPE Comparison", digits = 3),sanitize.text.function=function(x){x})


#Plot
pdf(file = paste0("./trained_jasa/", modelname, "_n_mse.pdf"), width = 6,height = 6)
t4.mse <- t4.summary.table[1,c(1,3,5,7)]
t4.mse <- log(t4.mse)
n.x <- log(c(1000,2500,5000,7500))
plot(t4.mse~n.x, type = "o", xlab = "log(n)",ylab ="log(MSE)", ylim = c(min(t4.mse),max(t4.mse)),
     main = expression("T"[4]))
#for(i in 2:7) lines(t4.mse[i,]~n.x, type = "o", col = i)
#legend("topright", lty = rep(1,7), col = 1:7, legend = c(rownames(t4.mse)[1:2], expression(tau == "0.05"),
#                                                         expression(tau == "0.25"),expression(tau == "0.5"),
#                                                         expression(tau == "0.75"),expression(tau == "0.95")))
dev.off()

#LS reg
ls.fit.t4 <- lm(t4.mse~n.x)


################################################################################################
################################################################################################
####                            Compare different width depth
################################################################################################
################################################################################################
my_code_path <- "../../GCDE/code"
sapply(list.files(my_code_path, full.names = TRUE), source)
show.diagnostic.pic <- FALSE
modelname <- "t4"
savepath <- "./trained_jasa/"

t4fun <- function(r){
    seed <- 131 * r
    if(r==9) seed <- 131*19
    cat("Start the", r, "simulation replicates...")
    exp.name <- c("wdhalf", "wdnormal","wddouble")
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
    
    
    gcde.t4.wdhalf <- fitgcde(data.t4$x, data.t4$y, data.t4$test.x, data.t4$test.y,
                              #load.disc = paste0(savepath,modelname,"/",modelname,"_",exp.name[1],
                              #                   "_discriminator_", r, ".h5"),
                              #load.gene = paste0(savepath,modelname,"/",modelname,"_",exp.name[1],
                              #                   "_generator_", r, ".h5"),
                              diagnostic.info    = data.t4$diagnostic.info,
                              show.diagnostic.pic= show.diagnostic.pic,
                              v.label            = "new",
                              kl.dual            = TRUE,
                              dual.form          = "fenchel",
                              iterations         = 7500,
                              dim.eta            = 4,
                              eta.dist           = "normal",
                              gene.nodes.num     = c(20,8),
                              dis.nodes.num      = c(25,13),
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
                              gene.regularizer.o = regularizer_l2(l = 0.0000),
                              disc.regularizer   = regularizer_l2(l = 0.0000),   
                              activation.disc    = "relu",
                              activation.gene    = "relu",
                              lambda.kl          = 20,
                              lambda.mse         = 0,       
                              m                  = 1, vanillaNN = FALSE)
    
    gcde.t4.wdnormal <- fitgcde(data.t4$x, data.t4$y, data.t4$test.x, data.t4$test.y,
                                #load.disc = paste0(savepath,modelname,"/",modelname,"_",exp.name[2],
                                #                   "_discriminator_", r, ".h5"),
                                #load.gene = paste0(savepath,modelname,"/",modelname,"_",exp.name[2],
                                #                   "_generator_", r, ".h5"),
                                diagnostic.info    = data.t4$diagnostic.info,
                                show.diagnostic.pic= show.diagnostic.pic,
                                v.label            = "new",
                                kl.dual            = TRUE,
                                dual.form          = "fenchel",
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
                                gene.regularizer.o = regularizer_l2(l = 0.0000),
                                disc.regularizer   = regularizer_l2(l = 0.0000),   
                                activation.disc    = "relu",
                                activation.gene    = "relu",
                                lambda.kl          = 20,
                                lambda.mse         = 0,       
                                m                  = 1, vanillaNN = FALSE)
    gcde.t4.wddouble <- fitgcde(data.t4$x, data.t4$y, data.t4$test.x, data.t4$test.y,
                                #load.disc = paste0(savepath,modelname,"/",modelname,"_",exp.name[3],
                                #                   "_discriminator_", r, ".h5"),
                                #load.gene = paste0(savepath,modelname,"/",modelname,"_",exp.name[3],
                                #                   "_generator_", r, ".h5"),
                                diagnostic.info    = data.t4$diagnostic.info,
                                show.diagnostic.pic= show.diagnostic.pic,
                                v.label            = "new",
                                kl.dual            = TRUE,
                                dual.form          = "fenchel",
                                iterations         = 7500,
                                dim.eta            = 4,
                                eta.dist           = "normal",
                                gene.nodes.num     = c(80,30,30),
                                dis.nodes.num      = c(100,50,25),
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
                                gene.regularizer.o = regularizer_l2(l = 0.0000),
                                disc.regularizer   = regularizer_l2(l = 0.0000),   
                                activation.disc    = "relu",
                                activation.gene    = "relu",
                                lambda.kl          = 20,
                                lambda.mse         = 0,       
                                m                  = 1, vanillaNN = FALSE)
    
    
    gcde.save.model.weights.jasa(list(gcde.t4.wdhalf, gcde.t4.wdnormal, gcde.t4.wddouble), modelname,exp.name,
                                 savepath, r)
    
    
    table.t4 <- generate.table(s = 2000, fit.gcde = gcde.t4.wdhalf, 
                               fit.gcde.extra = list(gcde.t4.wdnormal, gcde.t4.wddouble),
                               otherparameter = NULL, 
                               extra.names = exp.name)
    table.t4
}

clrep <- 10
cl <- makeCluster(clrep)
show.diagnostic.pic <- FALSE
clusterExport(cl, c("show.diagnostic.pic", "savepath", "modelname","my_code_path"),
              envir=environment())
clusterEvalQ(cl, sapply(list.files(my_code_path, full.names = TRUE), source))

t4.tables <- parLapply(cl, 1:clrep, t4fun)
stopCluster(cl)
save.image(paste0(savepath, modelname, "wd.RData"))



#process the tables from 10 simulation replications here:
t4.summary.table <- gcde.summary.table(t4.tables)
print(xtable(t4.summary.table, caption = "MSPE Comparison", digits = 3),sanitize.text.function=function(x){x})




################################################################################################
################################################################################################
####                            Compare different m
################################################################################################
################################################################################################
my_code_path <- "../../GCDE/code"
sapply(list.files(my_code_path, full.names = TRUE), source)
show.diagnostic.pic <- FALSE
modelname <- "t4"
savepath <- "./trained_jasa/"

t4fun <- function(r){
    seed <- 131 * r
    cat("Start the", r, "simulation replicates...")
    exp.name <- c("m1", "m3", "m5","m10")
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
    
    
    gcde.t4.m1 <- fitgcde(data.t4$x, data.t4$y, data.t4$test.x, data.t4$test.y,
                          #load.disc = paste0(savepath,modelname,"/",modelname,"_",exp.name[1],
                          #                   "_discriminator_", r, ".h5"),
                          #load.gene = paste0(savepath,modelname,"/",modelname,"_",exp.name[1],
                          #                   "_generator_", r, ".h5"),
                          diagnostic.info    = data.t4$diagnostic.info,
                          show.diagnostic.pic= show.diagnostic.pic,
                          v.label            = "new",
                          kl.dual            = TRUE,
                          dual.form          = "fenchel",
                          iterations         = 7500,
                          dim.eta            = 1,
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
                          gene.regularizer.o = regularizer_l2(l = 0.0000),
                          disc.regularizer   = regularizer_l2(l = 0.0000),   
                          activation.disc    = "relu",
                          activation.gene    = "relu",
                          lambda.kl          = 20,
                          lambda.mse         = 0,       
                          m                  = 1, vanillaNN = FALSE)
    
    gcde.t4.m2 <- fitgcde(data.t4$x, data.t4$y, data.t4$test.x, data.t4$test.y,
                          #load.disc = paste0(savepath,modelname,"/",modelname,"_",exp.name[2],
                          #                   "_discriminator_", r, ".h5"),
                          #load.gene = paste0(savepath,modelname,"/",modelname,"_",exp.name[2],
                          #                   "_generator_", r, ".h5"),
                          diagnostic.info    = data.t4$diagnostic.info,
                          show.diagnostic.pic= show.diagnostic.pic,
                          v.label            = "new",
                          kl.dual            = TRUE,
                          dual.form          = "fenchel",
                          iterations         = 7500,
                          dim.eta            = 3,
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
                          gene.regularizer.o = regularizer_l2(l = 0.0000),
                          disc.regularizer   = regularizer_l2(l = 0.0000),   
                          activation.disc    = "relu",
                          activation.gene    = "relu",
                          lambda.kl          = 20,
                          lambda.mse         = 0,       
                          m                  = 1, vanillaNN = FALSE)
    gcde.t4.m3 <- fitgcde(data.t4$x, data.t4$y, data.t4$test.x, data.t4$test.y,
                          #load.disc = paste0(savepath,modelname,"/",modelname,"_",exp.name[3],
                          #                   "_discriminator_", r, ".h5"),
                          #load.gene = paste0(savepath,modelname,"/",modelname,"_",exp.name[3],
                          #                   "_generator_", r, ".h5"),
                          diagnostic.info    = data.t4$diagnostic.info,
                          show.diagnostic.pic= show.diagnostic.pic,
                          v.label            = "new",
                          kl.dual            = TRUE,
                          dual.form          = "fenchel",
                          iterations         = 7500,
                          dim.eta            = 5,
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
                          gene.regularizer.o = regularizer_l2(l = 0.0000),
                          disc.regularizer   = regularizer_l2(l = 0.0000),   
                          activation.disc    = "relu",
                          activation.gene    = "relu",
                          lambda.kl          = 20,
                          lambda.mse         = 0,       
                          m                  = 1, vanillaNN = FALSE)
    
    gcde.t4.m4 <- fitgcde(data.t4$x, data.t4$y, data.t4$test.x, data.t4$test.y,
                          #load.disc = paste0(savepath,modelname,"/",modelname,"_",exp.name[4],
                          #                   "_discriminator_", r, ".h5"),
                          #load.gene = paste0(savepath,modelname,"/",modelname,"_",exp.name[4],
                          #                   "_generator_", r, ".h5"),
                          diagnostic.info    = data.t4$diagnostic.info,
                          show.diagnostic.pic= show.diagnostic.pic,
                          v.label            = "new",
                          kl.dual            = TRUE,
                          dual.form          = "fenchel",
                          iterations         = 7500,
                          dim.eta            = 10,
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
                          gene.regularizer.o = regularizer_l2(l = 0.0000),
                          disc.regularizer   = regularizer_l2(l = 0.0000),   
                          activation.disc    = "relu",
                          activation.gene    = "relu",
                          lambda.kl          = 20,
                          lambda.mse         = 0,       
                          m                  = 1, vanillaNN = FALSE)
    
    
    gcde.save.model.weights.jasa(list(gcde.t4.m1, gcde.t4.m2, gcde.t4.m3, gcde.t4.m4), modelname,exp.name,
                                 savepath, r)
    
    
    table.t4 <- generate.table(s = 2000, fit.gcde = gcde.t4.m1, 
                               fit.gcde.extra = list(gcde.t4.m2, gcde.t4.m3, gcde.t4.m4),
                               otherparameter = NULL, 
                               extra.names = exp.name)
    table.t4
}

clrep <- 10
cl <- makeCluster(clrep)
show.diagnostic.pic <- FALSE
clusterExport(cl, c("show.diagnostic.pic", "savepath", "modelname","my_code_path"),
              envir=environment())
clusterEvalQ(cl, sapply(list.files(my_code_path, full.names = TRUE), source))

t4.tables <- parLapply(cl, 1:clrep, t4fun)
stopCluster(cl)
save.image(paste0(savepath, modelname, "m.RData"))



#process the tables from 10 simulation replications here:
t4.summary.table <- gcde.summary.table(t4.tables)
print(xtable(t4.summary.table, caption = "MSPE Comparison", digits = 3),sanitize.text.function=function(x){x})

