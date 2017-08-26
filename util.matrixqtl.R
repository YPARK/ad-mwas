library(MatrixEQTL)

make.s.data <- function(mat, row.names, col.names) {
    stopifnot(is.matrix(mat))
    ret <- SlicedData$new()
    ret$CreateFromMatrix(mat)
    rownames(ret) <- row.names
    colnames(ret) <- col.names
    return(ret)
}

run.matrix.qtl <- function(X, Y, snp.names, out.names, me.out.file) {

    stopifnot(nrow(X) == nrow(Y))
    n <- nrow(X)
    individuals <- as.character(1:n)

    y.data <- make.s.data(t(Y), out.names, individuals)
    x.data <- make.s.data(t(X), snp.names, individuals)
    c.data <- SlicedData$new()

    me <- Matrix_eQTL_engine(snps = x.data,
                             gene = y.data,
                             cvrt = c.data,
                             pvOutputThreshold = 1,
                             noFDRsaveMemory = TRUE,
                             output_file_name = me.out.file)

    return(me)
}
