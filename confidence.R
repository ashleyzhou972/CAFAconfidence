getConf<-function(filepath){
  temp = unlist(strsplit(filepath,"/"))
  filename = temp[length(temp)]
  filename = unlist(strsplit(filename,'\\.'))
  filename[1]
  conf = read.table(filepath,sep="\t")
  conf[,1] = as.numeric(conf[,1])
  plot(conf[,1],conf[,2],pch = 20,xlab = "confidence",ylab = "TRUE Positive or not", main = filename)
  tp = conf[conf[,2]=="True",]
  fp = conf[conf[,2]=="False",]
  na = conf[is.na(conf[,2]),]
  par(mfrow = c(2,1))
  hist(tp[,1],breaks = 500, main = filename)
  hist(fp[,1],breaks = 500, main = filename)
  par(mfrow = c(1,1))
  summary(tp[,1])
  summary(fp[,1])
  ks.test(tp[,1],fp[,1],alternative = "less")
  boxplot(tp[,1],fp[,1],names = c("True postives","False postives"),ylab="Confidence",main=filename)
  newcol = unlist(lapply(conf[,2], FUN=conv))
  newconf = cbind(conf,newcol)
  mse = sum((newconf[,1]-newconf[,3])^2)/nrow(conf)
  return (mse)
}

conv<-function(conf){
  if (conf=="True") return(1)
  else return(0)
}


getConf('/home/nzhou/git/CAFAAssess/confidence/117/PaccanaroLab_1_9606_BPO_confdata.txt')
filepath='/home/nzhou/git/CAFAAssess/confidence/117/PaccanaroLab_1_9606_BPO_confdata.txt'


#New 06/07/2016 
#conf<-read.table('/home/nzhou/git/CAFAAssess/confidence/117/PaccanaroLab_1_9606_BPO_confdata.txt',sep='\t')
conv<-function(truth){
  if (truth=="True") return(1)
  else return(0)
}
mse<-function(pred_vec, truth_vec){
  if (length(pred_vec)==length(truth_vec)){
    return(sum((pred_vec-truth_vec)^2)/length(pred_vec))
  }
  else{
    print("vector not of the same length")
  }
}


getMSE<-function(filepath){
  conf = read.table(filepath,sep="\t")
  newcol = unlist(lapply(conf[,3], FUN=conv))
  index = rep(NA,nrow(conf))
  j=1
  index[j]=1
  for (i in 1:nrow(conf)){
    if (conf[index[j],1]!=conf[i,1]){
      j = j+1
      index[j] = i
    }
  }
  index=index[1:j]
  
  l = rep(NA,length(index)-1)
  for (i in 1:(length(index)-1)){
    l[i]=mse(conf[(index[i]:index[i+1]),2],newcol[index[i]:index[i+1]])
  }
  tp = conf[conf[,3]=="True",]
  fp = conf[conf[,3]=="False",]
  na = conf[is.na(conf[,2]),]
  png(paste(filepath,'_hist.png',sep=""))
  par(mfrow = c(2,1))
  hist(tp[,1],breaks = 500, main = filename)
  hist(fp[,1],breaks = 500, main = filename)
  dev.off()
  return(mean(l))
}

tianlab = read.table('/home/nzhou/git/CAFAAssess/confidence/85/tianlab_1_9606_BPO_confdata.txt',sep='\t')

getMSE(conf)
getMSE(tianlab)
naive = read.table('~/git/CAFAAssess/confidence/Naive/BN4S_bpo_confdata.txt',sep='\t')
blast = read.table('~/git/CAFAAssess/confidence/BLAST/BB4S_bpo_confdata.txt',sep='\t')
getMSE(naive)
getMSE(blast)


bpo = c('/117/PaccanaroLab_1_9606_BPO_confdata.txt','/85/tianlab_1_9606_BPO_confdata.txt','/129/tu_1_9606_BPO_confdata.txt','/127/txt_BPO_confdata.txt','/94/txt_BPO_confdata.txt','/83/bonneauLab_1_9606_BPO_confdata.txt','/75/Homo_sapiens_BPO_confdata.txt','/Naive/BN4S_bpo_confdata.txt','/BLAST/BB4S_bpo_confdata.txt')
mfo = c('/129/tu_1_9606_MFO_confdata.txt','/72/EVEX_1_9606_MFO_confdata.txt','/85/tianlab_1_9606_MFO_confdata.txt','/75/Homo_sapiens_MFO_confdata.txt','/127/txt_MFO_confdata.txt','/101/GO2Proto_1_9606_0_MFO_confdata.txt','/102/brennersifter_1_9606_0_MFO_confdata.txt','/94/txt_MFO_confdata.txt','/Naive/BN4S_mfo_confdata.txt','/BLAST/BB4S_mfo_confdata.txt')
cco = c('/72/EVEX_1_9606_CCO_confdata.txt','/85/tianlab_1_9606_CCO_confdata.txt','/129/tu_1_9606_CCO_confdata.txt','/115/Homo_sapiens_CCO_confdata.txt','/83/bonneauLab_1_9606_CCO_confdata.txt','/86/IASLAcademiaSinica_1_9606_CCO_confdata.txt','/127/txt_CCO_confdata.txt','/Naive/BN4S_cco_confdata.txt','/BLAST/BB4S_cco_confdata.txt')

#output in batches
bpo = c('/117/PaccanaroLab_1_9606_BPO_confdata.txt','/85/tianlab_1_9606_BPO_confdata.txt','/129/tu_1_9606_BPO_confdata.txt','/127/txt_BPO_confdata.txt','/94/txt_BPO_confdata.txt','/83/bonneauLab_1_9606_BPO_confdata.txt','/75/Homo_sapiens_BPO_confdata.txt','/Naive/BN4S_bpo_confdata.txt','/BLAST/BB4S_bpo_confdata.txt')
bpo_vec = rep(NA,length(bpo))
i=1
for (name in bpo){
  filepath = paste('/home/nzhou/git/CAFAAssess/confidence',name,sep="")
  a = read.table(filepath)
  bpo_vec[i]= getMSE(a)
  i=i+1
}

mfo = c('/129/tu_1_9606_MFO_confdata.txt','/72/EVEX_1_9606_MFO_confdata.txt','/85/tianlab_1_9606_MFO_confdata.txt','/75/Homo_sapiens_MFO_confdata.txt','/127/txt_MFO_confdata.txt','/101/GO2Proto_1_9606_0_MFO_confdata.txt','/102/brennersifter_1_9606_0_MFO_confdata.txt','/94/txt_MFO_confdata.txt','/Naive/BN4S_mfo_confdata.txt','/BLAST/BB4S_mfo_confdata.txt')
mfo_vec = rep(NA,length(mfo))
i=1
for (name in mfo){
  filepath = paste('/home/nzhou/git/CAFAAssess/confidence',name,sep="")
  a = read.table(filepath)
  mfo_vec[i]= getMSE(a)
  i=i+1
}


cco = c('/72/EVEX_1_9606_CCO_confdata.txt','/85/tianlab_1_9606_CCO_confdata.txt','/129/tu_1_9606_CCO_confdata.txt','/115/Homo_sapiens_CCO_confdata.txt','/83/bonneauLab_1_9606_CCO_confdata.txt','/86/IASLAcademiaSinica_1_9606_CCO_confdata.txt','/127/txt_CCO_confdata.txt','/Naive/BN4S_cco_confdata.txt','/BLAST/BB4S_cco_confdata.txt')
cco_vec = rep(NA,length(cco))
i=1
for (name in cco){
  filepath = paste('/home/nzhou/git/CAFAAssess/confidence',name,sep="")
  a = read.table(filepath)
  cco_vec[i]= getMSE(a)
  i=i+1
}