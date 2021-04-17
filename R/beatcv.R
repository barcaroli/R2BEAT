beat.cv <- function(n_file,stratif,errors,des_file,psu_file,rho,epsilon=10^(-11))

{
  
colnames(stratif) <- toupper(colnames(stratif))
colnames(errors) <- toupper(colnames(errors))
colnames(des_file) <- toupper(colnames(des_file))
colnames(rho) <- toupper(colnames(rho))
colnames(psu_file) <- toupper(colnames(psu_file))
colnames(n_file) <- toupper(colnames(n_file))

#parameters
nvar=ncol(errors)-1
nstrat=nrow(stratif)
strloop <- c(1:nstrat)
ndom=nrow(errors)

val<-NULL

nom_dom <-sapply(1:ndom, function(i) paste("DOM",i,sep=""))
dom     <-as.vector(stratif[,names(stratif) %in% nom_dom])
nvalues<-sapply( nom_dom, function(vari) {val<-c(val, nlevels(as.factor(dom[,vari])))})


stratif_completo<-stratif
disegno_completo<-des_file

########

sapply(1:nvar, function(i) { stratif_completo[,paste("S_ORIG",i,sep="")] <<- stratif_completo[,paste("S",i,sep="")] }  )

disegno_completo <- disegno_completo[order(disegno_completo$STRATUM),]
disegno_completo <- cbind(disegno_completo,n_file)

disegno_completo$F=(disegno_completo$SIZE/disegno_completo$STRAT_MOS)
disegno_completo$SOGLIA=( (disegno_completo$MINIMUM/disegno_completo$F)*disegno_completo$DELTA)

anagrafico1<-merge(disegno_completo[,c("STRATUM","SOGLIA")],psu_file,by="STRATUM")
anagrafico1$AR<-0
anagrafico1$AR[ (anagrafico1$PSU_MOS>=anagrafico1$SOGLIA)]<-1
anagrafico1$POPAR=anagrafico1$PSU_MOS*anagrafico1$AR
anagrafico1$POPNAR=anagrafico1$PSU_MOS-anagrafico1$POPAR
popolaz<-aggregate(anagrafico1[,c("POPAR","POPNAR")],by=list(STRATUM=anagrafico1$STRATUM), sum)

disegno_completo<-merge(disegno_completo, popolaz, by="STRATUM")
disegno_completo<-merge(disegno_completo,rho,by="STRATUM")

disegno_completo$CAMPAR=round((disegno_completo$F)*(disegno_completo$POPAR))
disegno_completo$CAMPNAR=disegno_completo$SIZE-disegno_completo$CAMPAR

disegno_completo$BDIS_AR=1*disegno_completo$DELTA
disegno_completo$BDIS_NAR=disegno_completo$MINIMUM*disegno_completo$DELTA

sapply(1:nvar, function(i) {
 disegno_completo[,paste("DEFF",i,sep="")] <<-
(( disegno_completo$SIZE/( disegno_completo$STRAT_MOS**2 ) )*
(
( ((disegno_completo$POPAR**2)/(disegno_completo$CAMPAR+epsilon)) * ( 1 + ((disegno_completo[,paste("RHO_AR",i,sep="")])*(disegno_completo$BDIS_AR-1) )) )
   +
( ((disegno_completo$POPNAR**2)/(disegno_completo$CAMPNAR+epsilon)) * ( 1 + ((disegno_completo[,paste("RHO_NAR",i,sep="")])*(disegno_completo$BDIS_NAR-1) )) )
 ))
} )


sapply(1:nvar, function(i)
{disegno_completo[,paste("DEFF",i,sep="")][ disegno_completo[,paste("DEFF",i,sep="")] < 0 ] <<- 1  })

sapply(1:nvar, function(i)
{disegno_completo[,paste("DEFT",i,sep="")] <<- sqrt(disegno_completo[,paste("DEFF",i,sep="")]) } )


stratif_completo<- merge(stratif_completo, disegno_completo[,c("STRATUM",sapply(1:nvar, function(i) (paste("DEFT",i,sep="") )))], by="STRATUM")


sapply(1:nvar, function(i)
{ stratif_completo[,paste("S",i,sep="")] <<- (stratif_completo[,paste("S_ORIG",i,sep="")]*disegno_completo[,paste("DEFT",i,sep="")]) } )

stratif_completo<-cbind(stratif_completo,n_file)


sapply(1:nvar, function(i){stratif_completo[,paste("varfin",i,sep="")] <<- ( stratif_completo[,paste("S",i,sep="")]*stratif_completo$N )**2*((1-stratif_completo$SIZE/stratif_completo$N )/stratif_completo$SIZE)})
sapply(1:nvar, function(i){stratif_completo[,paste("TOT",i,sep="")] <<- round(stratif_completo[,paste("M",i,sep="")]*stratif_completo$N) })


sapply(1:nvar, function(i){stratif_completo[,paste("varfin",i,sep="")] <<- ( stratif_completo[,paste("S",i,sep="")]*stratif_completo$N )**2*((1-stratif_completo$SIZE/stratif_completo$N )/stratif_completo$SIZE)})
sapply(1:nvar, function(i){stratif_completo[,paste("TOT",i,sep="")] <<- round(stratif_completo[,paste("M",i,sep="")]*stratif_completo$N) })

k=0
cv<-NULL
qu=NULL
for (j in (1:ndom))
for (i in (1:nvar))
{
k=k+1
statement_aggr<-paste( "cv[[",k,"]] <-aggregate(    stratif_completo[,c('TOT",i,"','N','varfin",i,"' )] ,by=list(Domain=c(stratif_completo$DOM" ,j,") ), sum) ", sep="" )
eval(parse(text=statement_aggr))
nome1=paste ("dom",j,"var",i,sep="")
statement_cv<-paste( "cv[[",k,"]]$cv_",nome1, " <- sqrt(cv[[",k,"]]$varfin",i,")/cv[[",k,"]]$TOT",i, sep="")
eval(parse(text=statement_cv))

names(cv)[k]<-nome1

#######################################################

statement_qu<-paste( "qu[[",k,"]] <- quantile( cv[[",k,"]]$cv_",nome1," , c(0,0.25,0.5,0.75,0.9,0.99,1))", sep="" )
eval(parse(text=statement_qu))

nome2=paste ("qu_dom",j,"var",i,sep="")

names(qu)[k]<-nome2

}
cv <<- cv
#assign("cv",cv,envir=.BaseNamespaceEnv)
qu <<- qu
#assign("qu",qu,envir=.BaseNamespaceEnv)

# New output (Barcaroli)

exp_cv <- NULL

st1 <- "exp_cv$Type <- c(rep("
st2 <- "exp_cv$Dom <- c(1:"

for (i in (1:ndom)) {
  if (i < ndom) st1 <- paste0(st1,"'DOM",i,"',",nvalues[i],"),rep(")
  if (i == ndom) st1 <- paste0(st1,"'DOM",i,"',",nvalues[i],"))")
  if (i < ndom) st2 <- paste0(st2,nvalues[i],",1:")
  if (i == ndom) st2 <- paste0(st2,nvalues[i],")")
}
eval(parse(text=st1))
eval(parse(text=st2))
exp_cv <- as.data.frame(exp_cv)

for (j in (1:nvar)) {
  st <- paste0("exp_cv$V",j," <- c(")
  for (i in (1:ndom)) {
    if (i < ndom) st <- paste0(st,"cv$dom",i,"var",j,"$cv_dom",i,"var",j,",")
    if (i == ndom) st <- paste0(st,"cv$dom",i,"var",j,"$cv_dom",i,"var",j,")")
  }
  eval(parse(text=st))
}
# return(cv)
return(exp_cv)

}