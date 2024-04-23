#CM.SESself is a signed stable SES
library(LoopAnalyst)
library(RConics)
library(reshape2)
library(ggplot2)

cSESself<-sign(t(CM.SESself))  #careful here row col swapped
diag(cSESself)[diag(cSESself)==0]<-(-1) ### temp fix to deal with singular matrix

det(cSESself)

c.adj<-RConics::adjoint(-cSESself)
absSESself<-abs(cSESself)

library(LoopAnalyst)
T<-make.T(cSESself)

# T<-cSESself

# for (i in 1:nrow(T)) {
	# for (j in 1:ncol(T)) {
		# T[i,j]<-permanent(absSESself[-i,-j]) 
	# }
# }

W<-abs(c.adj)/T

row.names(c.adj)<-row.names(CM.SESself)
colnames(c.adj)<-row.names(CM.SESself)

adj.melt<-melt(c.adj[8,])
adj.melt$Component=rownames(adj.melt)

WPosidonia<-melt(W[8,])
WPosidonia$Component=rownames(WPosidonia)

adj.melt$weight<-WPosidonia$value
ggplot(adj.melt,aes(x=Component,y=value,fill=weight))+
geom_bar(stat="identity")+
ylab("Posidonia response prediction")+
xlab("SES Component")+
theme_minimal()+
theme(axis.text.x = element_text(angle = 45, hjust = 1))
