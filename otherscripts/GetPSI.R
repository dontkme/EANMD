splicetype="SE"      #type of alternative splicing, e.g., SE, A3SS, A5SS, MXE, IR
counttype="JCEC"       #JCEC (junction count + exon body count) or JC (junction count only)

##################
#Input parameters#
##################
# inputpath="./02_PSI_value_quantification/01_Get_PSI_from_rMATS_output/example_input"     #input path for rMATS output
inputpath="/Users/kaininghu/HKN/0Postdoc/202108/HumanForebrain/combinedAll.od"     #input path for rMATS output"     #input path for rMATS output
outputpath="/Users/kaininghu/HKN/0Postdoc/202108/HumanForebrain/combinedAll.example_output"

inputcount=paste(counttype,".raw.input.",splicetype,".txt",sep="")
inputexon=paste("fromGTF.",splicetype,".txt",sep="")

command=paste("mkdir -p ",outputpath,sep="")    #create output folder if not exist
system(command)

b1path=inputpath
b1name="../Male_timecourse_brain_BAMs_list.txt2"

######################
#read in count result#
######################
setwd(inputpath)
RCtable=read.table(inputcount,sep="\t",header=T)
samplesize=length(strsplit(as.character(RCtable[1,"IJC_SAMPLE_1"]),split=",")[[1]])

#get exon name as row name
setwd(inputpath)
exoninfo=read.table(inputexon,sep="\t",header=T)
nospace_exoninfo=apply(exoninfo,2,function(x)gsub('\\s+', '',x)) #remove all the white spaces in exoninfo
exonname=apply(nospace_exoninfo,1,paste,collapse="|")

#get SRR ID as column name
setwd(b1path)
# raw_SRR=read.table(b1name,sep="\t")   
raw_SRR=read.table(b1name,sep=",",stringsAsFactors = F)   
SRRlist=rep(NA,dim(raw_SRR)[2])
for(i in 1:dim(raw_SRR)[2]){
  tmp <- strsplit(as.character(raw_SRR[1,i]),split="[.]")
  #SRRlist[i]=strsplit(as.character(raw_SRR[1,i]),split="[.]")[[1]][1]
  SRRlist[i]=paste(tmp[[1]][1],tmp[[1]][4],sep="_")
}

#generate tables for output
oriICtable=matrix(NA,dim(RCtable)[1],samplesize)    #original inclusion count
rownames(oriICtable)=exonname
colnames(oriICtable)=SRRlist      #we use SRR ID because it is unique to each sample

oriSCtable=matrix(NA,dim(RCtable)[1],samplesize)      #original skipping count
rownames(oriSCtable)=exonname
colnames(oriSCtable)=SRRlist

oritotalRCtable=matrix(NA,dim(RCtable)[1],samplesize)     #original total count
rownames(oritotalRCtable)=exonname
colnames(oritotalRCtable)=SRRlist

PSItable=matrix(NA,dim(RCtable)[1],samplesize)
rownames(PSItable)=exonname
colnames(PSItable)=SRRlist

#################
#fill the tables#
#################
get_psi<-function(x){  
  #get count
  inclulength=as.numeric(as.character(x[,"IncFormLen"]))
  skilength=as.numeric(as.character(x[,"SkipFormLen"]))
  ratio=inclulength/skilength
  inclucount=as.numeric(strsplit(as.character(x[,"IJC_SAMPLE_1"]),split=",")[[1]])   #original inclusion count
  skicount=as.numeric(strsplit(as.character(x[,"SJC_SAMPLE_1"]),split=",")[[1]])     #original skipping count
  #calculate psi
  PSI=(inclucount/inclulength)/(inclucount/inclulength+skicount/skilength)
  return(list(inclucount,   #original inclusion count
              skicount,     #original skipping count
              inclucount+skicount,     #original total count
              PSI,           #inclusion level
              inclucount,         #normalized inclusion count
              skicount*ratio,     #normalized skipping count (we use skipping count*ratio instead of inclusion count/ratio)
              c(inclucount+skicount*ratio)))   #normalized total count
}

for (i in 1:dim(RCtable)[1]){
  temp=get_psi(RCtable[i,])
  oriICtable[i,]=temp[[1]]
  oriSCtable[i,]=temp[[2]]
  oritotalRCtable[i,]=temp[[3]]
  PSItable[i,]=temp[[4]]
}


###########
#Filtering#
###########
IC=oriICtable
totalRC=oritotalRCtable
PSI=PSItable

#we filter out exons with average PSI<5% or >95% and average total read count <10
avepsi=apply(PSI,1,mean,na.rm=T)
averc=apply(totalRC,1,mean,na.rm=T)    #we use original total read count

rccutoff=10
bool=intersect(intersect(which(avepsi>0.05),which(avepsi<0.95)),which(averc>=rccutoff))
PSI_filter=PSI[bool,]   

#since PCA cannot deal with missing data, we further remove exons with more than 5% missing value and then impute the rest
sumna=function(x){
  temp=sum(is.na(x))
  return(temp)
}
sumNA=apply(PSI_filter,1,sumna)
missing_cutoff=0.05
PSI_filter=PSI_filter[which((sumNA/dim(PSI_filter)[2])<missing_cutoff),]

#we also filter out exons with max(PSI)-min(PSI)<=5%
deltapsi=function(x){
  temp=max(x,na.rm=T)-min(x,na.rm=T)   #if we use this, we will have 12132 exons left
  return(temp)
}
delta_psi=apply(PSI_filter,1,deltapsi)
PSI_filter=PSI_filter[which(delta_psi>0.05),]

#output result
setwd(outputpath)

IC_filter=IC[rownames(PSI_filter),]
totalRC_filter=totalRC[rownames(PSI_filter),]
write.table(PSI_filter,paste(splicetype,"_",counttype,"_PSI_filter.txt",sep=""),sep="\t",col.names = NA)
write.table(IC_filter,paste(splicetype,"_",counttype,"_IC_filter.txt",sep=""),sep="\t",col.names = NA)
write.table(totalRC_filter,paste(splicetype,"_",counttype,"_totalRC_filter.txt",sep=""),sep="\t",col.names = NA)



