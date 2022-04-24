gcde.summary.table <- function(table.list){
    num <- length(table.list)
    m <- ncol(table.list[[1]])
    n <- nrow(table.list[[1]])
    summarytable <- matrix(0, nrow = n, ncol = 2*m)
    for (i in 1:m) {
        coltable <- matrix(0, nrow = n, ncol = num)
        for(j in 1:num){
            coltable[ ,j] <- table.list[[j]][ ,i]
        }
        summarytable[ ,(2*i-1)] <- apply(coltable, 1, mean)
        summarytable[ ,(2 * i)] <- apply(coltable, 1, sd)/sqrt(num)
    }
    rownames(summarytable) <- rownames(table.list[[1]])
    colna <- rep("sd", 2*m)
    colna[seq(1, 2*m-1, by = 2)] <- colnames(table.list[[1]])
    colnames(summarytable) <- colna
    
    summarytable
}


