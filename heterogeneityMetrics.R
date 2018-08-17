#####helper functions############################
# preprocess mutation data for function "calculatePairwiseParam"
preprocessMuts <- function(data) {
  # first get rid of clonal mutations
  cols.ccf = grepl("ccf$",  colnames(data))
  test <- apply(data[,cols.ccf], 1, function(x) all(x>0.5))
  if (length(which(test))>0) {data <- data[-which(test),]}
  
  # also get rid of all mutations that are 0 across the samples we are interested in
  test <- apply(data[,cols.ccf], 1, function(x) all(x==0))
  if (length(which(test))>0) {data <- data[-which(test),]}
  
  # get rid of mutations that don't have a depth of 20 across all samples
  cols.refc = grepl("refc$",  colnames(data))
  cols.altc = grepl("altc$",  colnames(data))
  total.depth = data[,cols.refc]+ data[,cols.altc]
  test <- apply(total.depth, 1, function(x) all(x>=20))
  data <- data[which(test),]
  return(data)
}

### keep clonal mutations for these calculations
preprocessHFR <- function(data) {
  # get rid of mutations that don't have a depth of 20 across all samples
  cols.refc = grepl("refc$",  colnames(data))
  cols.altc = grepl("altc$",  colnames(data))
  total.depth = data[,cols.refc]+ data[,cols.altc]
  test <- apply(total.depth, 1, function(x) all(x>=20))
  data <- data[which(test),]
  return(data)
}

# one of the heterogeneity metrics
fst.hudson <- function(af, minAF=0.08) {  
  mafis = which(grepl("maf", colnames(af)))
  keep = as.vector(apply(af, 1, function(x, mafis) {
    maxmaf = max(as.numeric(x[mafis]))
    if (maxmaf > minAF) { TRUE
    } else {FALSE}
  }, mafis=mafis))     #filter data
  af = af[keep,]
  
  Ns = c()
  Ds = c()
  for(k in 1:nrow(af)) {
    n1 = af$depth1[k]
    n2 = af$depth2[k]
    p1 = af$maf1[k]
    p2 = af$maf2[k]
    N = (p1-p2)^2-(p1*(1-p1))/(n1-1)-(p2*(1-p2))/(n2-1)     # covariance
    D = p1*(1-p2)+p2*(1-p1)                                 # standard deviations
    Ns = c(Ns, N)
    Ds = c(Ds, D)
  }
  Fst.h = mean(Ns)/mean(Ds)
  return(Fst.h)
}

### calculate and output all summary statistics (heterogeneity metrics) 
# using maf and dpeth information
# sampAB is a matrix containing information for both samples
# snA and snB are the sample names for samples A and B
subclonalMut_otherCancer <- function(sampAB, snA, snB, minAF=0.08, statsAF=0.08, highAF=0.2, ssAF=0)  { 
  mafaAi = sampAB[,grep(paste(snA, "ccf$", sep=""), colnames(sampAB))]/2
  mafaBi = sampAB[,grep(paste(snB, "ccf$", sep=""), colnames(sampAB))]/2
  nbAi = 1
  nbBi = 1
  depthAi = sampAB[,grep(paste(snA, "refc", sep=""), colnames(sampAB))] + 
    sampAB[,grep(paste(snA, "altc", sep=""), colnames(sampAB))]
  depthBi = sampAB[,grep(paste(snB, "refc", sep=""), colnames(sampAB))] + 
    sampAB[,grep(paste(snB, "altc", sep=""), colnames(sampAB))]
  
  # subclonal mutations, and exclude cases where one sample had LOH
  subAi = which(mafaAi > minAF & ((mafaBi == 0 & (nbBi != 0 | nbAi == 0)) | mafaBi != 0) )
  mutsA = mafaAi[subAi]
  # sample specific mutations
  ssAi  = intersect(subAi, which( mafaAi > minAF & mafaBi <= ssAF ))
  
  subBi = which(mafaBi > minAF & ((mafaAi == 0 & (nbAi != 0 | nbBi == 0)) | mafaAi != 0) )
  mutsB = mafaBi[subBi]
  ssBi  = intersect(subBi, which( mafaBi > minAF & mafaAi <= ssAF ))
  
  KSD = as.numeric(ks.test( mutsA[which(mutsA > statsAF)], mutsB[which(mutsB > statsAF)] )$statistic)
  allSubRows = union(subAi,subBi)
  
  # for FST
  mutsSub = data.frame( maf1 = mafaAi[allSubRows], depth1=depthAi[allSubRows], maf2 = mafaBi[allSubRows], depth2=depthBi[allSubRows] )
  FST = fst.hudson(mutsSub, minAF=statsAF)
  
  # for other stats
  mutsA2 = mafaAi[intersect(subAi, which( mafaAi > statsAF ))]
  mutsAh2 = mafaAi[intersect(subAi, which( mafaAi > highAF ))]
  
  mutsASp2 = mafaAi[intersect(subAi, which( mafaAi > statsAF & mafaBi == 0))]
  mutsASph2 = mafaAi[intersect(subAi, which( mafaAi > highAF & mafaBi == 0))]
  
  mutsB2 = mafaBi[intersect(subBi, which(mafaBi > statsAF ))]
  mutsBh2 = mafaBi[intersect(subBi, which(mafaBi > highAF ))]
  
  mutsBSp2 = mafaBi[intersect(subBi, which(mafaBi > statsAF & mafaAi == 0))]
  mutsBSph2 = mafaBi[intersect(subBi, which(mafaBi > highAF & mafaAi == 0))]
  
  ratioHighSubA=length(mutsAh2)/length(mutsA2)
  ratioHighSubB=length(mutsBh2)/length(mutsB2)
  ratioHighSsA=length(mutsASph2)/length(mutsASp2)
  ratioHighSsB=length(mutsBSph2)/length(mutsBSp2)
  ratioSsA=length(mutsASp2)/length(allSubRows)
  ratioSsB=length(mutsBSp2)/length(allSubRows)
  ratioHRsA=length(mutsAh2)/length(allSubRows)
  ratioHRsB=length(mutsBh2)/length(allSubRows)
  
  # list for output    
  muts = list(rHighSub=mean(na.omit(c(ratioHighSubA,ratioHighSubB))), 
              rHighSs=mean(na.omit(c(ratioHighSsA,ratioHighSsB))),
              rSs=mean(c(ratioSsA,ratioSsB)), rHRs=mean(c(ratioHRsA,ratioHRsB)),
              FST=FST, KSD=KSD)
  return(muts)
}


##################### main functions #######################
# pairwise heterogeneity metrics used for ABC inference
# takes in the raw data (ccfs, mafs, and depths of coverage) as input
calculatePairwiseParams <- function(data){
  data = preprocessMuts(data)
  col.names = colnames(data)
  cols.ccf = grepl("ccf$",  col.names)
  data.ccf = data[,c(which(cols.ccf))];
  
  #### calculate stats in a pairwise manner (and then average) using function
  num_samples = ncol(data.ccf)
  comb_matrix = combn(num_samples,2)
  # create results vectors for all of the other stats
  result = matrix(nrow = ncol(comb_matrix), ncol=6) # because 6 parametes
  sample_names = substring(colnames(data.ccf), 0, nchar(colnames(data.ccf))-3)
  ## all possible combinations
  for (i in 1:ncol(comb_matrix)) {
    result[i,] = unlist(subclonalMut_otherCancer(data, sample_names[comb_matrix[1,i]], sample_names[comb_matrix[2,i]]))
  }
  # take the mean across all pairwise combos to come up with final values
  result_final=colMeans(result)
  return(result_final)
}

# v is a vector of the indices of the two CCFs to compare
# (or mafas, in which case need to change from 0.5 to 0.25 and 0.05 to 0.025)
calculateHFR <- function(data, v) {
  n <- 0
  totalhet <- 0
  for (i in 1:(length(v)-1)) {
    for (j in c((i+1):length(v))) {
      a <- subset(data, data[,v[i]] > 0.5 & data[,v[j]] < 0.05)
      b <- subset(data, data[,v[j]] > 0.5 & data[,v[i]] < 0.05)
      c <- subset(data, data[,v[i]] > 0.5 & data[,v[j]] > 0.5)
      het1 <- dim(a)[1] / (dim(a)[1]+dim(c)[1])
      het2 <- dim(b)[1]/ (dim(b)[1]+dim(c)[1])
      totalhet <- totalhet+het1+het2
      n <- n+2
    }
  }
  return(totalhet/n)
}

# v is a vector of the indices of the CCFs; the first one is treated as "pre" and the others as "post"
# (or mafas, in which case need to change from 0.5 to 0.25 and 0.1 to 0.05)
## this is for temporal HFR
calculatetHFR <- function(data, v) {
  clonal <- 0.5
  rare <- 0.1
  a <- data
  for (i in 2:(length(v))) {
    a <- subset(a, a[,v[i]] > clonal)
  }
  b <- subset(a, a[,v[1]] < rare)
  c <- subset(a, a[,v[1]] > clonal)
  return(dim(b)[1]/(dim(b)[1] + dim(c)[1]))
}
