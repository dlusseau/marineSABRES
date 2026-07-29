
boolean_file_creation_no_write<-function(SES.mat) {

SES.bin<-sign(SES.mat)
colnames(SES.bin)<-gsub("\\s*\\([^\\)]+\\)","",colnames(SES.bin))
row.names(SES.bin)<-gsub("\\s*\\([^\\)]+\\)","",row.names(SES.bin))
colnames(SES.bin)<-gsub("[^[:alnum:]]","",colnames(SES.bin))
row.names(SES.bin)<-gsub("[^[:alnum:]]","",row.names(SES.bin))

colnames(SES.bin)<-gsub(" ","_",colnames(SES.bin))
row.names(SES.bin)<-gsub(" ","_",row.names(SES.bin))

if (all(colSums(abs(SES.bin))!=0)==FALSE) { #fixed error 7 mar 2025
SES.bin<-SES.bin[,-which(colSums(abs(SES.bin))==0)]
}
targets<-colnames(SES.bin)
factors<-rep(NA,length(targets))
boolean.df<-data.frame(targets=targets,factors=factors) # to is the columns

for (i in 1:ncol(SES.bin)) {
poss<-names(which(SES.bin[,i]==1))
negs<-names(which(SES.bin[,i]==-1))
if (length(negs)>0) {
negs<-paste0("!",negs)
}
all<-c(poss,negs)

boolean.df$factors[i]<-paste(all,collapse="|")
}
boolean.df$targets<-as.character(boolean.df$targets)

bool_file<-array(NA,nrow(boolean.df)+1)
bool_file[1]<-paste(colnames(boolean.df),collapse=",")
for (k in 1:nrow(boolean.df)) {
bool_file[k+1]<-paste((boolean.df[k,c("targets","factors")]),collapse=",")
}
return(bool_file)
#write.csv(boolean.df,file=paste0(folder,filename,".csv"),row.names = F,quote=FALSE)
}

