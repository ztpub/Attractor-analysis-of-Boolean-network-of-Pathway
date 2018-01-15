library(GSVA)


setwd("G:\\TMP\\test\\PNA")
dat = read.csv("genenameDUKE.csv",header=FALSE, row.names=1)
# 

#################################################################################


ztexp = as.matrix(dat)


library(KEGGREST)
kL = keggLink("pathway", "hsa")
ptn = c()
for(i in 1:length(kL))
  ptn = c(ptn,kL[[i]])
ptn = unique(ptn)

ZT = c()
ZT$pathwaynames = c()
ZT$gs = c()
for(i in 1:length(ptn))
{
  print(i)
  print(length(ptn))
  kG = keggGet(ptn[i])
  gs = kG[[1]]$GENE
  tmp = c()
  if(length(gs)==0) 
    next
  for(j in 1:(length(gs)/2))
  {
    ss = gs[2*j]
    print(ss)
    xx = strsplit(ss,";")
    tmp = c(tmp,xx[[1]][1])
  }
  ZT$pathwaynames[i] = kG[[1]]$NAME
  ZT$gs[[kG[[1]]$NAME]] = tmp
}

PV = gsva(ztexp, ZT$gs,method=c("gsva"))
write.table(PV$es.obs,".\\gsva_bibduke.txt")
