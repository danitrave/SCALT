#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

FILE <- args[1]      #input file

#load("cellLinesCV_folds.RData")
load(FILE)
df <- data.frame(trainStart=rep(NA,length(Train_Idx)),trainEnd=rep(NA,length(Train_Idx)),
                                testStart=rep(NA,length(Train_Idx)),testEnd=rep(NA,length(Train_Idx)))
for (i in 1:length(Train_Idx)){
  #train coords
  trainStart <- Train_Idx[[i]][1]-1
  trainEnd <- Train_Idx[[i]][length(Train_Idx[[i]])]-1
  #test coords
  testStart <- Test_Idx[[i]][1]-1
  testEnd <- Test_Idx[[i]][length(Test_Idx[[i]])]-1
  df[i,]<-c(trainStart,trainEnd,testStart,testEnd)
}

write.table(df,"coordinates.tsv",sep="\t",quote = FALSE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
}
