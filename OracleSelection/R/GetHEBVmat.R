## Get HEBVs
#-------------------------------------------------------------------------
#Copyright (c) 2020 Antoine Allier
#-------------------------------------------------------------------------

# This code runs iterations of long-term genomic selection simulations to 

GetHEBVmat = function(ObjectGeno,ObjectBeta,ObjectMap,
                      WindowSize,StepSize,lambda = StepSize/WindowSize){
  # Matrix with loci effects for each individual in line and loci in column
  GEBVmat <- ObjectGeno*matrix(ObjectBeta,ncol = ncol(ObjectGeno),
                               nrow = nrow(ObjectGeno),byrow = TRUE)
  # Affect loci to haplotypes
  rownames(ObjectMap) <- ObjectMap$MARKER
  Haplo <- do.call(rbind,lapply(unique(ObjectMap$CHROMOSOME),function(chr_tmp){
    Lchr <- ObjectMap[intersect(ObjectMap$MARKER[ObjectMap$CHROMOSOME==chr_tmp],
                                colnames(GEBVmat)),"MARKER"]
    LHchr <- rollapply(Lchr,width = WindowSize,by = StepSize,align = "left",
                       FUN = function(d) {return(d)})
    rownames(LHchr) <- paste0("Hap_",seq(1,nrow(LHchr)),"_CHR_",chr_tmp)
    return(LHchr)
  }))
  # Design matrix Z affecting loci to haplotypes
  Z <- matrix(0,nrow = length(intersect(ObjectMap$MARKER,colnames(GEBVmat))),
              ncol = nrow(Haplo))
  rownames(Z) <- intersect(ObjectMap$MARKER,colnames(GEBVmat))
  colnames(Z) <- rownames(Haplo)
  invisible(lapply(1:ncol(Z),function(i){
    Z[Haplo[i,],i] <<- 1
  }))
  # Compute the HEBV matrix and mean position of each haplotypes
  HEBV <- (GEBVmat[,rownames(Z)]%*%Z)*lambda
  POS <- data.frame(HAP = colnames(Z),
                    CHR = sub(".*CHR_","",colnames(Z)),
                    POSITION = t(ObjectMap[rownames(Z),"POSITION"]%*%Z)/apply(Z,2,sum),
                    stringsAsFactors = FALSE)
  return(list(HEBV = HEBV,
              POSITION = POS))
}