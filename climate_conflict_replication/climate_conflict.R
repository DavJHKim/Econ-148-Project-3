
#generates Figures 1 and 2 and Table 2 of Burke, Miguel, Satyanath, Dykema, Lobell, PNAS 2009.

rm(list=ls())

library(maps)

#LOAD IN DATA
data=read.csv("/Documents/climate/conflict/data_for_distribution/clim_conflict_for_R.csv") #core data
load("/Documents/climate/conflict/data_for_distribution/RData_deltas2030monthly.Rdata") #climate change projections
ssa=read.csv('/Documents/climate/conflict/data_for_distribution/gseasSSA.csv') #growing season data
area=read.csv("/Documents/climate/conflict/data_for_distribution/croparea_weight.csv") #crop and land area from FAO
cltol=read.csv("/Documents/climate/conflict/data_for_distribution/precip_climatology.csv",row.names=1) #precipitation climatology
mss_country=read.csv('/Documents/climate/conflict/data_for_distribution/mss_countries_new.csv',header=T) #country names from Miguel, Satyanath, Sergenti 2004
mss_country[,1]=as.vector(mss_country[,1])
mss=as.vector(unique(mss_country[,1])) #list of unique country names


#regions are
	#1==Sahel
	#2==West Africa
	#3==Central Africa
	#4==East Africa
	#5==Southern Africa

#make area weights by region
z=which(area[,2]!="Cereals,Total +")
for (i in z) {
	area[i,5]=area[which(area[1:42,1]==area[i,1]),5]
	}
unweight=c()
for (j in 1:5) {
	unweight=c(unweight,sum(area[which(area[,2]=="land area"&area[,5]==j),4]))
	}
unweight=unweight/sum(unweight)


##climate model projections
#for temperature, averaging over monthly projected changes
#for precip, apply percentage changes to baseline climatatology, then compute total change in precip across months
clim=array(dim=c(6,3,2,20),dimnames=list(c(dimnames(altogether.monthly)[[1]][4:8],"AllSSA_unweight"),dimnames(altogether.monthly)[[3]],c("temp","precip"),dimnames(altogether.monthly)[[5]][1:20]))
k=1
for (i in 1:5) {
	for (j in 1:3) {
	clim[i,j,k,]=apply(altogether.monthly[i+3,,j,k,1:20],2,mean)
	}}
k=2
for (i in 1:5) {
	for (j in 1:3) {
	for (n in 1:20) {
		chg=altogether.monthly[i+3,,j,k,n]/100+1
		#apply chg to climatology, sum over months, calculate difference from climatology sum
		clim[i,j,k,n]=(sum(cltol[i,]*chg)/sum(cltol[i,])-1)*100
		}	
	}}	
clim[6,,,]=apply(clim[1:5,,,]*unweight,c(2,3,4),sum)


#Generate projections of income and regime change to 2030.
ystart=data[data$year_actual==1981,10]
yend=data[data$year_actual==2002,10]
ygro=exp(log(yend/ystart)/21)  #annual p.c. income growth over the historical period
incchg=c(median(yend,na.rm=T)*(1.02)^27,median(yend,na.rm=T)*median(ygro,na.rm=T)^27)-median(yend,na.rm=T)
polstart=data[data$year_actual==1981,228]
polend=data[data$year_actual==2002,228]
polchg=polend-polstart
pol2030=(polend+median(polchg))
pol2030[pol2030>10]=10
polchg=mean(pol2030-polend)


#bootstrap climate-conflict regressions, using different climate variables or controls
hold=list()
hold[[1]]=data.frame(data[,c(106,42,45,110:150,178:218)]) #levels with country-spec time trends, temp only
hold[[2]]=data.frame(data[,c(106,42,45,54,57,110:150,178:218)]) #levels with country-spec time trends, temp and prec
hold[[3]]=data.frame(data[,c(106,42,45,10,228,110:150,177)]) #levels with simple time trend, temp only, gdp and polity lagged
hold[[4]]=data.frame(data[,c(106,42,45,54,57,10,228,110:150,177)]) #levels with simple time trend, temp and prec, gdp and polity lagged
hold[[5]]=data.frame(data[,c(106,42,45,54,57,110:150,177)]) #levels with simple time trends, temp and prec
hold[[6]]=data.frame(data[,c(106,48,51,60,63,110:150,178:218)])  #diffs w/ country time trend, temp and precip
hold[[7]]=data.frame(data[,c(106,219:222,110:150,178:218)]) #unweight dev trend, temp and precip, country time trend
hold[[8]]=data.frame(data[,c(106,42,45,11,12,110:150,178:218)])  #unweight levels with gpcp precip, country time trend 
hold[[9]]=data.frame(data[,c(106,70,77,71,78,110:150,178:218)])  #NCC levels w/ country time trends
for (j in 1:9) {
aa=hold[[j]]
aa=aa[which(!is.na(aa[,2])),]  #drop NAs so we don't sample them
coef=c()
for (i in 1:10000) {
	spl=sample(1:dim(aa)[1],replace=T)  #sample observations with replacement
	bb=aa[spl,]
	model=lm(formula(bb),bb)
	coef=rbind(coef,model$coefficients[2:7])  #collect coefficient estimates
	if (substr(i,nchar(i)-2,nchar(i))=="000") {print(i)} #show progress
	}
save(coef,file=paste("/Documents/climate/conflict/boot/bootstrap_",j,".Rdata",sep=""))
}



#FIGURE 1
#combine bootstrap coefficients with clim projections for a1b, by region, plot middle 95% of observations for climate uncertainty only, regression uncertainty only, and both together
plt=function(x,y) {  #function boxplots the middle 95% of obs for a given distribution, at the assigned spot
	n=length(x)
	all95=sort(x)*100
	all95=all95[(n-0.95*n):(0.95*n)]  #5-95% of distribution
	boxplot(all95+loc[i,1]+setzero,range=0,at=loc[i,2]+y,col=color[i],horizontal=T,add=T,boxwex=4,axes=F)
	}
color=c("red","orange","green","blue","grey","white")
names=c("Sahel","West Africa","Central Africa","Eastern Africa","Southern Africa","Sub-Saharan Africa")
quartz(height=7,width=14)
par(mfrow=c(1,2),mar=c(2,2,2,2))
plot(1,type="n",xlim=c(-25,50),ylim=c(-40,30),xlab="",ylab="",axes=F)
box()
countrycols=c(5,5,2,4,1,3,3,3,3,4,4,2,3,2,2,2,4,2,5,5,5,1,1,5,1,2,2,4,5,2,2,4,1,2,4,4,2,5,5,5,5) #assign countries to regions
ssa1=as.character(as.vector(ssa[c(1:10,12:42),2]))
ssa1[12]="Gambia"  #rename so map() function will work. 
for (i in 1:length(ssa1)) {
	map(regions=ssa1[i],exact=T,add=T,col=color[countrycols[i]],fill=T,axes=F)
	}
locx=c(10,-25,-10,22,18,-25)  #gives xmin start for grey rectangles
locy=c(12,8,-12,-8,-30,-36)  #gives ymin start for grey rectangles
loc=data.frame(locx,locy)
for (i in 1:6) {
	projt=clim[i,1,1,]*10
	projp=clim[i,1,2,]
	rect(loc[i,1]-1,loc[i,2]-2,loc[i,1]+27,loc[i,2]+15,col="grey95",density=70,border=NA)  #background rectangle
	#temp change
	setzero=8  #how far to the right of the grey rectangle xmin you want 0.0 to be
	text(seq(loc[i,1]+(setzero-5),loc[i,1]+20+(setzero-5),5),loc[i,2]+1,seq(0,2,0.5),cex=0.8) #axis labels
	segments(loc[i,1]+(setzero-5)-1,loc[i,2]+2,loc[i,1]+20+(setzero-5)+1,loc[i,2]+2) #axis line
	segments(seq(loc[i,1]+(setzero-5),loc[i,1]+20+(setzero-5),5),loc[i,2]+1.5,seq(loc[i,1]+(setzero-5),loc[i,1]+20+(setzero-5),5),loc[i,2]+2) #axis ticks
	boxplot(projt+loc[i,1]+3,range=0,at=loc[i,2]+3.5,col=color[i],horizontal=T,add=T,boxwex=4,axes=F)
	#precip change
	setzero=8  #how far to the right of the grey rectangle xmin you want 0.0 to be
	text(seq(loc[i,1]+(setzero-5),loc[i,1]+20+(setzero-5),5),loc[i,2]+7,seq(-5,15,5),cex=0.8) #axis labels
	segments(loc[i,1]+(setzero-5)-1,loc[i,2]+8,loc[i,1]+20+(setzero-5)+1,loc[i,2]+8) #axis line
	segments(seq(loc[i,1]+(setzero-5),loc[i,1]+20+(setzero-5),5),loc[i,2]+7.5,seq(loc[i,1]+(setzero-5),loc[i,1]+20+(setzero-5),5),loc[i,2]+8) #axis ticks
	text(loc[i,1]+setzero+5,loc[i,2]+13.5,names[i],cex=1.3)
	boxplot(projp+loc[i,1]+8,range=0,at=loc[i,2]+10,col=color[i],horizontal=T,add=T,boxwex=4,axes=F,density=10)
	}

#conflict projections
plot(1,type="n",xlim=c(-25,50),ylim=c(-40,30),xlab="",ylab="",axes=F)
box()
countrycols=c(5,5,2,4,1,3,3,3,3,4,4,2,3,2,2,2,4,2,5,5,5,1,1,5,1,2,2,4,5,2,2,4,1,2,4,4,2,5,5,5,5) #assign countries to regions
ssa1=as.character(as.vector(ssa[c(1:10,12:42),2]))
mss=as.character(unique(data[,1]))
matchcol=countrycols[match(mss,ssa1)]
matchcol[8]=2  #fix countries with non-matched names. ivory coast
matchcol[20:21]=3  #congos
region=c()  #assign each MSS country to a region
for (n in 1:dim(data)[1]) {
	region=c(region,matchcol[which(mss==as.character(data[n,1]))])
	}
ssa1[12]="Gambia"  
for (i in 1:length(ssa1)) {
	map(regions=ssa1[i],exact=T,add=T,col=color[countrycols[i]],fill=T,axes=F)
	}
locx=c(10,-25,-10,22,18,-25)  #gives xmin start for grey rectangles
locy=c(12,8,-12,-8,-30,-36)  #gives ymin start for grey rectangles
loc=data.frame(locx,locy)
setzero=8  #how far to the right of the grey rectangle xmin you want 0.0 to be
load("/Documents/climate/conflict/boot/bootstrap_1.Rdata")
for (i in 1:6) {
	projt=clim[i,1,1,]
	if (i<6) {projp=clim[i,1,2,]*mean(data[which(region==i),56],na.rm=T)/100} #evaluate precip at regional avg precip
		else {projp=clim[i,1,2,]*mean(data[,56],na.rm=T)/100}  #use all precip for all SSA
	projt=projt[!is.na(projt)]
	projp=projp[!is.na(projp)]
	all=c()
	for (j in 1:length(projt)) {
		all=c(all,projt[j]*coef[,1]+projt[j]*coef[,2])
		}
	climr=median(coef[,1])*projt+median(coef[,2])*projt
	confr=median(projt)*coef[,1]+median(projt)*coef[,2]
	rect(loc[i,1]-1,loc[i,2]-2,loc[i,1]+27,loc[i,2]+15,col="grey95",density=70,border=NA)  #background rectangle
	text(seq(loc[i,1]+(setzero-5),loc[i,1]+20+(setzero-5),5),loc[i,2]+1,seq(-5,15,5),cex=0.8) #axis labels
	segments(loc[i,1]+(setzero-5)-1,loc[i,2]+2,loc[i,1]+20+(setzero-5)+1,loc[i,2]+2) #axis line
	segments(seq(loc[i,1]+(setzero-5),loc[i,1]+20+(setzero-5),5),loc[i,2]+1.5,seq(loc[i,1]+(setzero-5),loc[i,1]+20+(setzero-5),5),loc[i,2]+2) #axis ticks
	text(loc[i,1]+setzero+5,loc[i,2]+13.5,names[i],cex=1.3)
	plt(climr,4)
	plt(confr,7)
	plt(all,10)
	text(rep(loc[i,1]+setzero-6.5,3),loc[i,2]+c(4,7,10),c("(3)","(2)","(1)"))
	}
dev.copy2eps(file="/Documents/climate/conflict/plots/Figure1.eps")


#Table 2
hold=list()
load("/Documents/climate/conflict/boot/bootstrap_1.Rdata")
hold[[1]]=coef
load("/Documents/climate/conflict/boot/bootstrap_2.Rdata")
hold[[2]]=coef
load("/Documents/climate/conflict/boot/bootstrap_4.Rdata")
hold[[3]]=coef
load("/Documents/climate/conflict/boot/bootstrap_4.Rdata")
hold[[4]]=coef
projections=matrix(nrow=6,ncol=5)
rownames(projections)=c("A1btemp","A1btemp_prec","A2temp","A2temp_prec","B1temp","B1temp_prec")
colnames(projections)=c("median","med_inc","5th_inc","95th_inc","% obs<0")
base=mean(data$war_prio_new[data$year_actual<2003])  #baseline conflict incidence
where=matrix(c(1:6),nrow=2)
for (i in 1:2) {
	coeff=hold[[i]]
	k=6
	for (x in 1:3) {
	projt=clim[k,x,1,]
	projp=clim[k,x,2,]*mean(data[,54],na.rm=T)/100  #evaluate precip at mean continental precip
	projt=projt[!is.na(projt)]
	projp=projp[!is.na(projp)]
	all=c()
	if (i==1) {coeff[,3:6]=0}
	if (i==2) {coeff[,5:6]=0}
	for (j in 1:length(projt)) {
		all=c(all,projt[j]*(coeff[,1]+coeff[,2])+projp[j]*(coeff[,3]+coeff[,4]))
		}
	n=length(all)
	all=sort(all)*100
	all10=all[(n-0.95*n):(0.95*n)] #5th-95th of distribution
	projections[where[i,x],]=c(median(all10),median(all10)/base,min(all10)/base,max(all10)/base,sum(all<0)/length(all)*100)
	}}
write.csv(projections,"/Documents/climate/conflict/plots/Table2.csv")



#FIGURE 2
hold=list()
load("/Documents/climate/conflict/boot/bootstrap_1.Rdata")
hold[[1]]=coef
load("/Documents/climate/conflict/boot/bootstrap_2.Rdata")
hold[[2]]=coef
load("/Documents/climate/conflict/boot/bootstrap_6.Rdata")
hold[[3]]=coef
load("/Documents/climate/conflict/boot/bootstrap_7.Rdata")
hold[[4]]=coef
load("/Documents/climate/conflict/boot/bootstrap_8.Rdata")
hold[[5]]=coef
load("/Documents/climate/conflict/boot/bootstrap_9.Rdata")
hold[[6]]=coef
load("/Documents/climate/conflict/boot/bootstrap_4.Rdata")
hold[[7]]=coef
load("/Documents/climate/conflict/boot/bootstrap_4.Rdata")
hold[[8]]=coef
quartz(height=8,width=8)
plot(1,type="n",ylim=c(0,length(hold)+0.5),xlim=c(-15,18),axes=F,ylab="",xlab="percentage point change in civil war",cex.axis=1.1,cex.lab=1.1)
axis(1,at=seq(-5,15,5),labels=seq(-5,15,5),cex.axis=1.1,cex.lab=1.1)
abline(v=0,lty=3)
nms=c("Temperature only","Temperature and precipitation","Temperature and precipitation,","Temperature and precipitation,","Temperature (CRU) and precipitation (GPCP)","Temperature and precipitation (NCC)","Temperature and precipitation, linear","Temperature and precipitation,")
nms2=c("(Table 1, Model 1)","(Table 1, Model 2)"," first differences (Table S2, Model 3)","deviations from trend (Table S2, Model 5)","(Table S3, Model 1)","(Table S3, Model 3)","extrapolation of 1981-02 income and","optimistic scenario for income and")
nms3=c(rep("",6),"democracy trends (Table 1, Model 3)","democracy trends (Table 1, Model 3)")
count=length(hold)
k=6
inc=0
where=c(count:1)-0.1
where[where<2.75]=where[where<2.75]-0.7 #lowering the bottom two down a bit
for (i in 1:count) {
	coeff=hold[[i]]
	projt=clim[k,1,1,]
	projp=clim[k,1,2,]*mean(data[,56],na.rm=T)/100  #evaluate precip at mean continental precip
	projt=projt[!is.na(projt)]
	projp=projp[!is.na(projp)]
	all=c()
	if (i==1) {coeff[,3:6]=0}
	if (i%in%2:6) {coeff[,5:6]=0}
	if (i==7) {inc=incchg[2]} 
	if (i==8) {inc=incchg[1]} 
	for (j in 1:length(projt)) {
		all=c(all,projt[j]*(coeff[,1]+coeff[,2])+projp[j]*(coeff[,3]+coeff[,4])+inc*coeff[,5]+polchg*coeff[,6])
		}
	n=length(all)
	all=sort(all)*100
	all10=all[(n-0.95*n):(0.95*n)] #5th-95th of distribution
	boxplot(all10,range=0,at=where[i],horizontal=TRUE,add=TRUE,col="grey",axes=F)
	text(-16,where[i],nms[i],cex=0.8,pos=4)
	text(-15,where[i]-0.2,nms2[i],cex=0.8,pos=4)
	text(-15,where[i]-0.4,nms3[i],cex=0.8,pos=4)
	}
rect(-30,length(hold)+0.3,20,length(hold)+1,col="light grey")
text(-16,8.55,"Impact of climate change, other factors fixed",pos=4,cex=1.1)
rect(-30,1.75,20,2.25,col="light grey")
text(-16,2,"Combined impacts of changes in climate, per capita income, and democracy",pos=4,cex=1.1)
box()
dev.copy2eps(file="/Documents/climate/conflict/plots/Figure2.eps")

