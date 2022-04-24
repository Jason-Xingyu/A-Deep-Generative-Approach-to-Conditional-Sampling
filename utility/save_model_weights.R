gcde.save.model.weights <- function(gcde.list, modelname, folderpath, r){
    save_model_weights_hdf5(gcde.list[[1]]$generator, paste0(folderpath, modelname, "/", modelname,
                                                            "_primal_generator_", r, ".h5"))
    save_model_weights_hdf5(gcde.list[[1]]$discriminator, paste0(folderpath, modelname, "/", modelname,
                                                                 "_primal_discriminator_", r, ".h5"))
    save_model_weights_hdf5(gcde.list[[2]]$generator, paste0(folderpath, modelname, "/", modelname,
                                                             "_fenchel_generator_", r, ".h5"))
    save_model_weights_hdf5(gcde.list[[2]]$discriminator, paste0(folderpath, modelname, "/", modelname,
                                   "_fenchel_discriminator_", r, ".h5"))
    save_model_weights_hdf5(gcde.list[[3]]$generator, paste0(folderpath, modelname, "/", modelname,
                                                             "_donsker_generator_", r, ".h5"))
    save_model_weights_hdf5(gcde.list[[3]]$discriminator, 
                            paste0(folderpath, modelname, "/", modelname,
                                   "_donsker_discriminator_", r, ".h5"))
}

#For use in JASA response
gcde.save.model.weights.jasa <- function(gcde.list, modelname, exp.name, folderpath, r){
    nn <- length(gcde.list)
    for (i in 1:nn) {
        save_model_weights_hdf5(gcde.list[[i]]$generator, paste0(folderpath, modelname, "/", modelname,
                                                                 "_",exp.name[i],
                                                                 "_generator_", r, ".h5"))
        save_model_weights_hdf5(gcde.list[[i]]$discriminator, paste0(folderpath, modelname, "/", modelname,
                                                                     "_",exp.name[i],
                                                                     "_discriminator_", r, ".h5"))
    }
}

gcde.save.model.weights.wgangp <- function(gcde.fit, modelname, folderpath, r){
    save_model_weights_hdf5(gcde.fit$generator, paste0(folderpath, modelname, "/", modelname,
                                                             "_","wgangp",
                                                             "_generator_", r, ".h5"))
    save_model_weights_hdf5(gcde.fit$discriminator, paste0(folderpath, modelname, "/", modelname,
                                                                 "_","wgangp",
                                                                 "_discriminator_", r, ".h5"))
}



gcde.save.one.model.weights <- function(gcde.fit, modelname, folderpath){
    save_model_weights_hdf5(gcde.fit$generator, paste0(folderpath,"/", modelname, "_generator.h5"))
    save_model_weights_hdf5(gcde.fit$discriminator, paste0(folderpath,"/", modelname, "_discriminator.h5"))
}
