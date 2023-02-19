## Calculate H-score
#-------------------------------------------------------------------------
#Copyright (c) 2019 Antoine Allier
#-------------------------------------------------------------------------

ComputeH = function(ObjectHEBVmat,ListPopE,ListPopD){
  # Evaluate only the elite haplotypes
  if(is.null(ListPopD)) {
    tmp <- apply(ObjectHEBVmat[ListPopE,],2,max) #maximum of the columns
    output <- data.frame(H = sum(tmp),
                         LINE = "PopE",
                         stringsAsFactors = FALSE)
  } else {
    # Evaluate the elite and donor haplotypes
    output <- do.call(rbind,lapply(ListPopD,function(i){
      tmp <- apply(ObjectHEBVmat[c(ListPopE,i),],2,max)
      return(data.frame(H = sum(tmp),
                        LINE = i,
                        stringsAsFactors = FALSE))
    }))
  }
  return(output)
}