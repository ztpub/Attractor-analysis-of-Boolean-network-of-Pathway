library(BoolNet)


setwd("H:\\TMP\\test\\R1\\")
setwd("C:\\Users\\Administrator\\Desktop\\test\\R1\\")


getNetMatrix <- function(network)
{
    library('stringr')
  
    genes = network$genes
    ints <- network$interactions
    
    gN = length(genes)
    
    Mat = matrix((1:gN*gN)*0, nrow=gN, ncol=gN) 
    
    for(i in 1:gN)
    {
        ele <- ints[[i]]
        ss <- ele$input
        #print(ss)
        
        for(j in 1:length(ss))
        {
            Mat[ss[j],i] = 1
        }
    
    }
    
    names(Mat) <- genes
    MatVec <- stack(Mat)
    
    
    return(MatVec[,1])
}


testBasinSizes <- function(network, accumulate=TRUE, params)
{
    attr <- getAttractors(network)
    basinSizes <- sapply(attr$attractors, function(a)
    {
       a$basinSize
    })
   
    if (accumulate)
        return(mean(basinSizes))
    else
        return(basinSizes)
}

# 
myAttractorMat <- function(fi1,fi2){
  #################################################################################

  data= read.csv("preboolnet.csv",header=FALSE)
  bin <- binarizeTimeSeries(data, method = "scanStatistic", windowSize = 0.005*fi1, sign.level = 0.01*fi2)
  sn=20;#样本个数
  tn=15;
  sx=c(3,4,5,6,7,9,11,15,16,20)#SX10
  asx=c(2,8,13,14,18,19)#6asx
  
  t=c(1:tn)
  tp=matrix(nrow=tn,ncol=sn)
  for (j in 1:tn)
  {
    s=c()
    for (i in 1:sn)
    {
      temp=(i-1)*tn+t(j)
      s=c(s,temp)
    }
    tp[j,]=s#tp=15*20
  }
  
  test1 <-c()
  nu<-c()
  
  netMat <- c()    

  #par(mfrow=c(1,4))
  for(i in 1:sn)
  {
    #print("---------------------------------------")
    sam <- bin$binarizedMeasurements[,tp[,i]]
    net <- reconstructNetwork(sam, method="bestfit", maxK=3, returnPBN=TRUE)
    net1 <- chooseNetwork(net, rep(1, length(net$genes)))
    
    
    
    rmat <- getNetMatrix(net1)
    netMat <- rbind(netMat,rmat)    
    
    attractors <- getAttractors(net1)  
    
    
    for(j in 1:length(attractors$attractors))
    {
      #print(getAttractorSequence(attractors, j))
      #print(attractors$attractors[[j]]$basinSize)
    }
    
    NumAtt = length(attractors$attractors)
    
    MaxAtt = 0
    StateAtt = c()
    for(j in 1:length(attractors$attractors))
    {
      if(attractors$attractors[[j]]$basinSize>MaxAtt)
      {
        MaxAtt = attractors$attractors[[j]]$basinSize
        StateAtt = getAttractorSequence(attractors, j)
      }
    }
    test1<-rbind(test1,StateAtt)
    nu<-c(nu,nrow(StateAtt))
    #print(c(i,NumAtt,MaxAtt,nrow(StateAtt)))#number,nu of attractor, largest states,attracter state
  }
  test1<-t(test1)
  
  
  CA <- dist(netMat,'euclidean')
  #print(CA)
  hc1<-hclust(CA)#,"mcquitty")
  out.id=cutree(hc1,k=2)                 #得到分为3类的数值
  accNet = max((length(which(out.id[sx]==1))+length(which(out.id[asx]==2)))/16,(length(which(out.id[sx]==2))+length(which(out.id[asx]==1)))/16)
  
  
  CA<-matrix(nrow=sn,ncol=sn)
  
  for (k in 1:sn)
  {
    for (l in k:sn)
    {
      
      # test1<-scale(test1)
      if(k==1)
      {
        m=0
      }
      if(k>1)
      {
        m=sum(nu[1:k-1])
      }
      if(l==1)
      {
        n=0
      }
      if(l>1)
      {
        n=sum(nu[1:l-1])
      }
      
      V1 = test1[,(m+1):(m+nu[k])]
      if(nu[k]-1==1)
      {
        V1 = t(t(V1))
      }
      V2 = test1[,(n+1):(n+nu[l])]
      if(nu[l]-1==1)
      {
        V2 = t(t(V2))
      }
      # 
      # if(l==10)
      # {
      #  print(dim(V1))
      #  print(V2)
      #  print(t(t(V2)))
      #  print(ca<-cancor(V1,V2))
      # }
      
      if(is.null(dim(V1)))
      {
          CA[l,k] <- 1
          next
      }
      
      if(is.null(dim(V2)))
      {
          CA[l,k] <- 1
          next
      }
      
      ca<-cancor(V1,V2)
      CA[l,k]<-max(ca$cor)
      
      # tryCatch({
      #   ca<-cancor(V1,V2)
      #   CA[l,k]<-max(ca$cor)
      #   if(l==10)
      #   {
      #     print(CA[l,k])
      #     print(V2)
      #   }
      # },error=function(e){print(dim(V2))})#CA[l,k] <- 1
      
      
      #CA[l,k]<-CA[k,l]
    }  
  }
  
  #print(CA)
  CA<-as.dist(CA)
  #print(CA)
  hc1<-hclust(CA)#,"mcquitty")
  out.id=cutree(hc1,k=2)                 #得到分为3类的数值
  acc = max((length(which(out.id[sx]==1))+length(which(out.id[asx]==2)))/16,(length(which(out.id[sx]==2))+length(which(out.id[asx]==1)))/16)
  

  rM <- c(fi1, fi2, acc, accNet, 0, 0, 0, 0, 0, 0, 0)
  print(rM)

  if(acc>0.9)
  {
    evl1 <- c()
    evl2 <- c()
    evl3 <- c()
    evl4 <- c()
    evl5 <- c()
    evl6 <- c()
    evl7 <- c()


    for(i in 1:sn)
    {
      #print("---------------------------------------")
      sam <- bin$binarizedMeasurements[,tp[,i]]
      net <- reconstructNetwork(sam, method="bestfit", maxK=3, returnPBN=TRUE)
      net1 <- chooseNetwork(net, rep(1, length(net$genes)))
      
      #print(i)
      #print(net1)
      #par(mfrow=c(1,1))
      #plotNetworkWiring(net1)
      
      # compare the in-degrees of the states in the
      # cell cycle network to random networks
      pv <- testNetworkProperties(net1, testFunction="testIndegree", alternative="greater")
      evl1 <- c(evl1,pv$pval)
      # compare the robustness of attractors in the cell cycle network
      # to random networks by perturbing the networks
      pv <- testNetworkProperties(net1, testFunction="testAttractorRobustness",
                                  testFunctionParams=list(perturb="functions", numSamples=10),
                                  alternative="greater")
      evl2 <- c(evl2,pv$pval)  
      
      
      # compare the robustness of attractors in the cell cycle network
      # to random networks by perturbing the state trajectories
      pv <- testNetworkProperties(net1, testFunction="testAttractorRobustness",
                                  testFunctionParams=list(perturb="trajectories", numSamples=10),
                                  alternative="greater")
      evl3 <- c(evl3,pv$pval)  
      
      # compare the robustness of single state transitions in the cell cycle network
      pv <- testNetworkProperties(net1, testFunction="testTransitionRobustness",
                                  testFunctionParams=list(numSamples=10),
                                  alternative="less")
      evl4 <- c(evl4,pv$pval)  
      
      pv <- testNetworkProperties(net1,
                                  numRandomNets=100,
                                  testFunction="testBasinSizes",
                                  xlab="Average size of basins of attraction")
      evl5 <- c(evl5,pv$pval)     
      
      
      pr <- perturbTrajectories(net1,
                                measure="hamming",
                                numSamples=100,
                                flipBits=1)
      evl6 <- c(evl6,pr$value)  
      
      pr <- perturbTrajectories(net1,
                                measure="attractor",
                                numSamples=100,
                                flipBits=1)
      evl7 <- c(evl7,pr$value)    
      
    }

    tc1 <- t.test(evl1[sx],evl1[asx])
    tc2 <- t.test(evl2[sx],evl2[asx])
    tc3 <- t.test(evl3[sx],evl3[asx])
    tc4 <- t.test(evl4[sx],evl4[asx])
    tc5 <- t.test(evl5[sx],evl5[asx])
    tc6 <- t.test(evl6[sx],evl6[asx])
    tc7 <- t.test(evl7[sx],evl7[asx])
    
    rM <- c(fi1,fi2,acc, accNet, tc1$p.value, tc2$p.value, tc3$p.value, tc4$p.value, tc5$p.value, tc6$p.value, tc7$p.value)
    print(rM)
  }
  
  return(rM)
}

#rM = myAttractorMat(33, 10)



Rmat = c()
for(fi1 in 20:40)   #20:40
{
  for(fi2 in 5:15)  # 5:15
  {
    tryCatch({
      #print(c(fi1,fi2))
      Rlist = myAttractorMat(fi1, fi2)
      Rmat = rbind(Rmat, Rlist)
    },error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
  }
}

library(vioplot)
x1 <- Rmat[,3]
x2 <- Rmat[,4]
vioplot(x1, x2, names=c("State", "Structure"), col="gold")
title("Violin Plots of efficiency")

write.table(Rmat,"Rmat.txt")