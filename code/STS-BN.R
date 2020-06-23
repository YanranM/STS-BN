library(Rfast)

####################################
#  Functions we are going to use   #
####################################
collevel<-function(aldata){
  clev=c()
  nc<-ncol(aldata)
  for(i in 1:nc){
    clev=c(clev,length(unique(aldata[,i])))
  }
  return(clev)
}


dof<-function(datafr){
  datafr=data.frame(datafr)
  sn=ncol(datafr)
  cl=collevel(datafr[,1:2])
  df=(cl[1]-1)*(cl[2]-1)
  if(sn>2){
    ###
    datach=data.frame(v1=paste0(datafr[,1],datafr[,2]),v2=do.call(paste0,datafr[c(-1,-2)]))
    chartable=table(datach)
    chartt=colSums(chartable!=0)
    chartt=chartt-cl[1]-cl[2]+1
    chartt[chartt<0]=0
    df=sum(chartt)
    ###
  }
  else{
    tabledata=table(datafr)
    ec=sum(table(datafr)==0)
    df=df-ec
  }
  return(df)
}

######## stop at df=0 ##########
Phase_I<-function(t,aldata,bv1){     #t:target node, aldata:input data, bv1: threshold value 
  #print(Sys.time())
  np=ncol(aldata)
  v<-1:np
  flag=0
  MB=c()
  xi=0
  data_mat=apply(as.matrix(aldata),2,as.numeric)
  while(flag==0){
    g=1
    for(i in 1:np){
      if(i != t & !is.element(i,MB)){
        dc_n=length(MB)+1
        if(dc_n==1){
          gt_df=dof(data_mat[,c(t,i)])
          if(gt_df != 0){
            gt=g2Test_univariate(data_mat[,c(t,i)],c(3,3))
            gt_p=1-pchisq(gt$statistic,gt_df)
          }
          else{
            gt_p=2   ###impossible value
          }
          ###GT####
        }
        else{
          temp_data=data_mat[,c(t,i,MB)]
          gt_df=dof(temp_data)
          if(gt_df!=0){
            gt=g2Test(temp_data,1,2,3:(dc_n+1),rep(3,dc_n+1)) 
            gt_p=1-pchisq(gt$statistic,gt_df)
          }
          else{
            gt_p=2  ###impossible value
          }
          ###GT###
        }
        if(gt_p<g){
          g=gt_p
          xi=i
        }
      }
    }
    if(g<bv1){
      MB=union(xi,MB)
      #add to MT
    }
    else{
      flag=1
    }
  }
  return(MB)
}




#######
nulist<-function(n){
  l=list(c())
  while(n-1>0){
    n=n-1
    l=c(l,list(c()))
  }
  return(l)
}

################################
#         Phase_II             #
################################
biccoe=0.17
s_singlenodeBIC<-function(x,pa){
  if(is.null(pa)){
    q_len=1
    k=length(levels(x))
    k_level=levels(x)
    Nijk=matrix(nrow = k,ncol = q_len)  ###matrix Nijk  
    pd=matrix(nrow = k,ncol = q_len)  #P(D|thetas,S)
    for(i in 1:k){
      Nijk[i]=sum(x==levels(x)[i])
    }
    for(i in 1:k){
      pd[i]=Nijk[i]*log(Nijk[i]/sum(Nijk))
    }
    cs=q_len*(k-1)
    BIC_score=sum(pd)-biccoe*cs*log(length(x))
  }
  ###################################
  else{
    aldata<-data.frame(x,pa)
    k=length(levels(x))
    k_level=levels(x)
    p=ncol(pa)
    if(is.null(p)){
      p=1
      pa_level=list(levels(pa))
    }
    else{
      pa_level=list()
      for(i in 1:p){
        pa_level=c(pa_level,list(levels(pa[,i])))
      }
    }
    pa_q=expand.grid(pa_level)
    q_len=nrow(pa_q)
    Nijk=matrix(nrow = k,ncol = q_len)  ###matrix Nijk
    data_char=do.call(paste0,aldata)
    pa_q_char=do.call(paste0,pa_q)
    
    
    for(i in 1:k){
      #sample size of particular configuration
      for(j in 1:q_len){
        ###character
        temp_char=paste0(k_level[i],pa_q_char[j])
        Nijk[i,j]=sum(data_char==temp_char)
      }
    }
    pd=matrix(nrow = k,ncol = q_len)  #P(D|thetas,S)
    for(i in 1:k){
      for(j in 1:q_len){
        if(Nijk[i,j]==0){
          pd[i,j]=0
        }
        else{
          pd[i,j]=Nijk[i,j]*log(Nijk[i,j]/sum(Nijk[,j]))
        }
      }
    }
    cs=q_len*(k-1)
    BIC_score=sum(pd)-biccoe*cs*log(nrow(aldata))
    
  }
  return(BIC_score)
}

#######
s_overallBIC<-function(ndata,lis){
  s=0
  for(i in 1:length(lis)){
    if(is.null(lis[[i]])){
      s=s+s_singlenodeBIC(ndata[,i],c())
    }
    else{
      s=s+s_singlenodeBIC(ndata[,i],ndata[,lis[[i]]])
    }
    
  }
  return(s)
}


#############
Phase_II<-function(x,pacan,v){
  temp_data=data.frame(pacan,x)
  n=ncol(pacan)
  p=v
  if(is.null(n)){
    n=1
  }
  plis=nulist(n+1)
  plis[n+1]=list(p)
  if(is.null(v)){
    i=0
  }
  else{
    i=max(v)
  }
  temps1=s_overallBIC(temp_data,plis)
  for(q in (i+1):n){
    v2=union(v,q)
    plis[n+1]=list(v2)
    temps2=s_overallBIC(temp_data,plis)
    s2=temps2
    if(temps2>temps1){
      if(q<n){                       
        sp2=Phase_II(x,pacan,v2)
        s2=sp2[[1]]
        v2=sp2[[2]]
      }
    }
    if(s2>temps1){
      temps1=s2
      p=v2
    }
  }
  als1p=list(temps1,p)
  return(als1p)
  
}

##############################
#         Run STS-BN         #
##############################
tdata=read.csv("example.csv",colClasses ="factor")
#In this data, "A" and "B" are disease associated sites.
bv1=0.05
CNS=c()
CNS=Phase_I(1,tdata,bv1)
CNS=rev(CNS)
PCt=Phase_II(tdata[,1],tdata[,CNS],c())
PCname=colnames(tdata)[CNS[PCt[[2]]]]
print(PCname)
