##############################
####### SBM Version 1 ########
##############################
#required packages
library(biostrings)
library(seqinr)
library(stringdist)
############################################################################################
## Truncate and Concatenate Seqeunces for each H1N1 Antigenic Site (Sa, Sb, Ca1, Ca2, Cb) ##
############################################################################################
#Use known antigenic site amino acid positions to truncate and concatenate proteins seqeunces to contain only amino acids in the antigenic site
mySeqs  <- read.alignment(file = "/myLocation/AllSeqs.fasta", format = "fasta")#read in aligned sequences  (e.g. MUSCLE, MAAFT, clustalx)

for (i in 1:length(mySeqs$nam)) {
 	name<-mySeqs$nam[[i]]
 	seq<-mySeqs$seq[[i]]
 	assign(name, seq)
 	#Sa epitope; truncate and concatenate seqeunces	
	a<-substr(seq, 141, 142)
	b<-substr(seq, 170, 174)
	c<-substr(seq, 176, 181)
	Sa<-paste0(a,b,c)
	epiNam<-paste0(name,"_Sa")
	assign(epiNam, Sa)
	if(i == 1){
		SaVir<-epiNam}
		else {SaVir<-c(SaVir,epiNam)}
	#Sb epitope; truncate and concatenate seqeunces	
	Sb<-substr(seq, 201, 212)
	epiNam<-paste0(name,"_Sb")
	assign(epiNam, Sb)
	if(i == 1){
		SbVir<-epiNam}
		else {SbVir<-c(SbVir,epiNam)}
	#Ca1 epitope; truncate and concatenate seqeunces	
	d<-substr(seq, 183, 187)
	e<-substr(seq, 220, 222)
	f<-substr(seq, 252, 254)
	Ca1<-paste0(d,e,f)
	epiNam<-paste0(name,"_Ca1")
	assign(epiNam, Ca1)
	if(i == 1){
		Ca1Vir<-epiNam}
		else {Ca1Vir<-c(Ca1Vir,epiNam)}
	#Ca2 epitope; truncate and concatenate seqeunces	
	g<-substr(seq, 154, 159)
	h<-substr(seq, 238, 239)
	Ca2<-paste0(g,h)
	epiNam<-paste0(name,"_Ca2")
	assign(epiNam, Ca2)
	if(i == 1){
		Ca2Vir<-epiNam}
		else {Ca2Vir<-c(Ca2Vir,epiNam)}
	#Cb epitope; truncate and concatenate seqeunces	
	Cb<-substr(seq, 87, 92)
	epiNam<-paste0(name,"_Cb")
	assign(epiNam, Cb)
	if(i == 1){
		CbVir<-epiNam}
		else {CbVir<-c(CbVir,epiNam)}
 	}
 	
####################################################################
## Calculate Hamming Distance matrix for each H1N1 Antigenic Site ##
####################################################################
#create data frames of antigenic distances (i.e. Hamming Distance divided by antigenic site length, multiplied by 20) for each 
##Sa
mat_Sa<-matrix(data=NA, nrow=length(SaVir), ncol=length(SaVir))
df_Sa<-as.data.frame(mat_Sa)
colnames(df_Sa)<-SaVir
rownames(df_Sa)<-SaVir
for (i in 1:length(SaVir)){
	name<-paste0(SaVir[i],"_vec")
	v1<-get(SaVir[i])
	for (j in 1:length(SaVir)){
		v2<-get(SaVir[j])
		dis<-stringdist(v1, v2, method="h")
		numAA<-nchar(v1)
		perdif<-dis/numAA
		AD<-perdif*20
		ADint<-round(AD)
		if(j==1){vec<-ADint} else{vec<-c(vec, ADint)}
		assign(name,vec)}
		if(j==length(SaVir)){df_Sa[i]<-vec}
		
	}
	
##Sb
mat_Sb<-matrix(data=NA, nrow=length(SbVir), ncol=length(SbVir))
df_Sb<-as.data.frame(mat_Sb)
colnames(df_Sb)<-SbVir
rownames(df_Sb)<-SbVir
for (i in 1:length(SbVir)){
	name<-paste0(SbVir[i],"_vec")
	v1<-get(SbVir[i])
	for (j in 1:length(SbVir)){
		v2<-get(SbVir[j])
		dis<-stringdist(v1, v2, method="h")
		numAA<-nchar(v1)
		perdif<-dis/numAA
		AD<-perdif*20
		ADint<-round(AD)
		if(j==1){vec<-ADint} else{vec<-c(vec, ADint)}
		assign(name,vec)}
		if(j==length(SbVir)){df_Sb[i]<-vec}
		
	}

##Ca1
mat_Ca1<-matrix(data=NA, nrow=length(Ca1Vir), ncol=length(Ca1Vir))
df_Ca1<-as.data.frame(mat_Ca1)
colnames(df_Ca1)<-Ca1Vir
rownames(df_Ca1)<-Ca1Vir
for (i in 1:length(Ca1Vir)){
	name<-paste0(Ca1Vir[i],"_vec")
	v1<-get(Ca1Vir[i])
	for (j in 1:length(Ca1Vir)){
		v2<-get(Ca1Vir[j])
		dis<-stringdist(v1, v2, method="h")
		numAA<-nchar(v1)
		perdif<-dis/numAA
		AD<-perdif*20
		ADint<-round(AD)
		if(j==1){vec<-ADint} else{vec<-c(vec, ADint)}
		assign(name,vec)}
		if(j==length(Ca1Vir)){df_Ca1[i]<-vec}
		
	}
	
##Ca2
#create data frames with for single epitope comparisons
mat_Ca2<-matrix(data=NA, nrow=length(Ca2Vir), ncol=length(Ca2Vir))
df_Ca2<-as.data.frame(mat_Ca2)
colnames(df_Ca2)<-Ca2Vir
rownames(df_Ca2)<-Ca2Vir
for (i in 1:length(Ca2Vir)){
	name<-paste0(Ca2Vir[i],"_vec")
	v1<-get(Ca2Vir[i])
	for (j in 1:length(Ca2Vir)){
		v2<-get(Ca2Vir[j])
		dis<-stringdist(v1, v2, method="h")
		numAA<-nchar(v1)
		perdif<-dis/numAA
		AD<-perdif*20
		ADint<-round(AD)
		if(j==1){vec<-ADint} else{vec<-c(vec, ADint)}
		assign(name,vec)}
		if(j==length(Ca2Vir)){df_Ca2[i]<-vec}
		
	}

##Cb
#create data frames with for single epitope comparisons
mat_Cb<-matrix(data=NA, nrow=length(CbVir), ncol=length(CbVir))
df_Cb<-as.data.frame(mat_Cb)
colnames(df_Cb)<-CbVir
rownames(df_Cb)<-CbVir
for (i in 1:length(CbVir)){
	name<-paste0(CbVir[i],"_vec")
	v1<-get(CbVir[i])
	for (j in 1:length(CbVir)){
		v2<-get(CbVir[j])
		dis<-stringdist(v1, v2, method="h")
		numAA<-nchar(v1)
		perdif<-dis/numAA
		AD<-perdif*20
		ADint<-round(AD)
		if(j==1){vec<-ADint} else{vec<-c(vec, ADint)}
		assign(name,vec)}
		if(j==length(CbVir)){df_Cb[i]<-vec}
		
	}
	

#########################################
## Write Out Antigenic Distance Matrix ##
#########################################
write.csv(df_Sa, file="Sa_AD.csv")
write.csv(df_Sb, file="Sb_AD.csv")	
write.csv(df_Cb, file="Cb_AD.csv")
write.csv(df_Ca1, file="Ca1_AD.csv")
write.csv(df_Ca2, file="Ca2_AD.csv")	
#write out composite antigenic distance
allDF<-(df_Sa + df_Sb + df_Ca1 + df_Ca2 + df_Cb)
write.csv(allDF, file="Composite_AD.csv")


