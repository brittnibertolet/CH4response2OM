###### Analyses for Bertolet et al. 
###### Lake sediment methane responses to organic matter are related to 
###### microbial community composition in experimental microcosms

###### Updated 2021-12-08

###### Calculate CH4 and CO2 production rates from incubations ####
# Read in cleaned and quality controlled CH4 and CO2 data 
# Raw data were adjusted to account for headspace dilution with sampling 
# See adjustHS.R for function
ydata=read.csv("github/data/gcAdjData_22lakes.csv", stringsAsFactors = F )
#write.csv(ydata,"github/data/gcAdjData_22lakes.csv", quote = F, row.names = F)
# Calculate production rate using linear regression of the timeseries
# Loop through each lake, treatment, and replicate combination for both CO2 and CH4
lakes=unique(ydata$lakeID)
treat=unique(ydata$treatment)

rates=data.frame()
for(i in 1:length(lakes)){
  for(j in 1:length(treat)){
    for(k in 1:3){
      gc=ydata[ydata$lakeID==lakes[i] & ydata$treatment==treat[j] & ydata$rep==k,]
      CH4fit=lm(gc$CH4adjConc_umolL~gc$incub_days)
      CO2fit=lm(gc$CO2adjConc_umolL~gc$incub_days)
      rates=rbind(rates,data.frame(lakeID=lakes[i], treatment=treat[j], rep=k, 
                               CH4prod=CH4fit$coefficients[2], 
                               CO2prod=CO2fit$coefficients[2]))
    }
  }
}

# Set incubations in which CH4 did not accumulate over time to 0
rates$CH4prod[rates$CH4prod<0]=0

# Kruskall wallace test to test for the effect of OM additions across all lakes
kruskal.test(CH4prod~treatment, data=rates)

rates$treatment=factor(rates$treatment, levels=c("control", "algal"), labels = c("control", "OM addition"))
p1a=ggplot(rates, aes(x=treatment, y=CH4prod))+
  #geom_boxplot(aes(fill=treatment))+
  stat_summary(geom="bar", fun = "mean", aes(fill=treatment))+
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.2)+
  ylab(expression(paste(CH[4], " production (", mu,"mol ", L^-1," ", d^-1,")")))+
  xlab("Treatment")+
  ylim(0, 35)+
  theme_bw()+
  theme(panel.grid=element_blank(),
        plot.margin=unit(c(0.2,0,0.2,0.2),"cm"))+
  scale_fill_manual(values=c("grey60", "darkolivegreen3"))+
  guides(fill="none")
p1b=ggplot(rates, aes(x=lakeID, y=CH4prod))+
  stat_summary(geom="bar", fun = "mean", aes(fill=treatment), position = "dodge")+
  stat_summary(geom="errorbar", fun.data = "mean_se", aes(group=treatment), position = "dodge")+
  scale_fill_manual(name="Treatment",values=c("grey70", "darkolivegreen3"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5))+
  xlab("Lake ID")+
  ylim(0, 35)+
  guides(fill="none")+
  ylab(expression(paste(CH[4], " production (", mu,"mol ", L^-1," ", d^-1,")")))

plot_grid(p1a, p1b, ncol=2, rel_widths = c(0.5, 1), align="h", axis="l",
          labels = c("A", "B"))
ggsave("github/figures/Figure1.pdf", height=3.25, width=6.5)

# Get the difference in CH4 and CO2 between each
diff=data.frame(LakeID=lakes, CH4diff=NA, CO2diff=NA, baseCH4=NA)
for(i in 1:length(lakes)){
  controlCH4=mean(rates$CH4prod[rates$lakeID==lakes[i] & rates$treatment=="control"])  
  algalCH4=mean(rates$CH4prod[rates$lakeID==lakes[i] & rates$treatment=="OM addition"])  
  controlCO2=mean(rates$CO2prod[rates$lakeID==lakes[i] & rates$treatment=="control"])  
  algalCO2=mean(rates$CO2prod[rates$lakeID==lakes[i] & rates$treatment=="OM addition"])  
  # calculate the difference and store 
  diff[i,2]=algalCH4-controlCH4
  diff[i,3]=algalCO2-controlCO2
  diff[i,4]=controlCH4
}

# Maximum potential utilization of algae
diff$maxAutilizied=(diff$CH4diff*28*0.1/1040)*100
range(diff$maxAutilizied)


###### Generate metrics of microbial community composition ####
# read in microbial metrics 
ddPCR=read.csv("github/data/ddPCRdata_22lakes.csv", stringsAsFactors = F)
# read in meta data
meta=read.csv(file ="github/data/metaData_22lakes.csv", stringsAsFactors = F)
# read in OTU table and taxonomy 
otus=read.table("github/data/3M.opti_mcc.shared", sep="\t", header=T, stringsAsFactors = F)
tax=read.table("dgithub/data/3M.taxonomy", sep="\t", header=T, stringsAsFactors = F)
rownames(otus)=otus$Group

# make sure OTUs are in more than one sample
# throw out singletons 
csums=colSums(otus[,-(1:3)])
singleOTUs=names(csums[csums==1])
otus=otus[, !colnames(otus)%in%singleOTUs]

# rarify to lowest sample 
min(rowSums(otus[,-(1:3)]))
otus2=rrarefy(otus[,-(1:3)], sample=39853)
rowSums(otus2)

# split otu table into methanogens and non-methanogens 
methanogens=tax$OTU[grep("Methano", tax$Taxonomy)]
motus=otus2[,colnames(otus2)%in%methanogens]
botus=otus2[,!(colnames(otus2)%in%methanogens)]

#look at read count and relative abundances
counts=data.frame(LakeID=gsub("3M_", "", rownames(motus)), 
                  nonmeth_reads=rowSums(botus),
                  meth_reads=rowSums(motus),
                  meth_relAbund=rowSums(motus)/rowSums(botus))
range(counts$meth_relAbund)*100
sum(counts$meth_reads)/(sum(counts$nonmeth_reads)+sum(counts$meth_reads))*100
# merge with microbial metrics
microbes=merge(ddPCR, counts, by="LakeID")

# transform read counts to relative abundances 
botu.ra=decostand(botus, method="total", MARGIN = 1) 
rowSums(botu.ra)
motu.ra=decostand(motus, method="total", MARGIN = 1) 
rowSums(motu.ra)

#Generate distance matrix
distB=vegdist(botu.ra, method="bray")
distM=vegdist(motu.ra, method="bray")
#Generate principle coordinate scaling with eigenvalues
pcoaB=cmdscale(distB, eig=T)
pcoaM=cmdscale(distM, eig=T)
#Create new dataframe 
pcoa.DF=data.frame(LakeID=gsub("3M_", "",rownames(pcoaM$points)),
                   B.axis1=pcoaB$points[,1],
                   B.axis2=pcoaB$points[,2],
                   M.axis1=pcoaM$points[,1],
                   M.axis2=pcoaM$points[,2])
rownames(pcoa.DF)=NULL

# Statistical relationships
adonis(distM~data$sediment_pH)
adonis(distM~data$percOM)
adonis(distB~data$sediment_pH)
adonis(distB~data$percOM)

# merge with microbial metrics
microbes=merge(microbes,pcoa.DF, by="LakeID")
# merge with meta data 
DF=merge(microbes, meta, by="LakeID")
# merge with difference in CH4 production data frame
DF=merge(diff, DF, by="LakeID")

###### Linear regressions to determine best predictor of response of CH4 production to OM ####
summary(lm(log(CH4diff)~copies_mcrA_per_g, data=DF))
summary(lm(log(CH4diff)~copies_16S_per_g, data=DF))
summary(lm(log(CH4diff)~B.axis1, data=DF))
summary(lm(log(CH4diff)~M.axis1, data=DF))
summary(lm(log(CH4diff)~sediment_pH, data=DF))
summary(lm(log(CH4diff)~percOM, data=DF))
summary(lm(log(CH4diff)~chl.a, data=DF))
summary(lm(log(CH4diff)~baseCH4, data=DF))
summary(lm(log(CH4diff)~baseCH4, data=DF))


# Figure 2
p2a=ggplot(DF, aes(x=B.axis1, y=log(CH4diff)))+
  stat_smooth(method="lm", formula=y~x, color="black")+
  geom_point(aes(color=sediment_pH))+
  ylab(expression(paste("log(",Delta*CH[4]," production)")))+
  xlab("non-MCC PCoA 1")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_continuous(name="pH")
p2b=ggplot(DF, aes(x=M.axis1, y=log(CH4diff)))+
  stat_smooth(method="lm", formula=y~x, color="black")+
  geom_point(aes(color=sediment_pH))+
  ylab(expression(paste("log(",Delta*CH[4]," production)")))+
  xlab("MCC PCoA 1")+
  theme_bw()+
  theme(panel.grid = element_blank())

plot_grid(p2a+guides(color="none"),
          p2b+guides(color="none"),
          get_legend(p2a), ncol=3, rel_widths = c(1,1,0.3))
ggsave("github/figures/Figure2.pdf", height=2.75, width=6.5)

##### Supplementary Figures ####
# Figure 1 
gcSub1=c("BH",  "BO",  "BR",  "BY",  "CB",  "CR",  "CYB", "HB",  "HH",  "JS",  "MI")

ydata1=ggplot(ydata[ydata$lakeID%in%gcSub1,], aes(x=incub_days, y=CH4adjConc_umolL))+geom_point(aes(color=as.factor(rep)))+
  facet_grid(treatment~lakeID)+
  stat_smooth(method="lm", formula=y~x, se=F, aes(color=as.factor(rep)))+
  guides(color="none")+
  ylab(expression(paste(CH[4], " concentration (", mu,"mol ", L^-1,")")))+
  xlab("Incubation Day")+
  ylim(-50, 1100)
ydata2=ggplot(ydata[!(ydata$lakeID%in%gcSub1),], aes(x=incub_days, y=CH4adjConc_umolL))+geom_point(aes(color=as.factor(rep)))+
  facet_grid(treatment~lakeID)+
  stat_smooth(method="lm", formula=y~x, se=F, aes(color=as.factor(rep)))+
  guides(color="none")+
  ylab(expression(paste(CH[4], " concentration (", mu,"mol ", L^-1,")")))+
  xlab("Incubation Day")+
  ylim(-50, 1100)
all=plot_grid(ydata1, ydata2, ncol=1)
ggsave("github/figures/SuppFigure1.pdf", height=5.5, width=9)

# Figure 2
nonMCC=ggplot(DF, aes(x=B.axis1, y=B.axis2))+geom_point(aes(color=sediment_pH, size=percOM))+
  xlab("non-MCC PCoA 1 [29.2%]")+
  ylab("non-MCC PCoA 2 [11.6%]")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.box = "horizontal")
MCC=ggplot(DF, aes(x=M.axis1, y=M.axis2))+geom_point(aes(color=sediment_pH, size=percOM))+
  xlab("MCC PCoA 1 [40.6%]")+
  ylab("MCC PCoA 2 [12.6%]")+
  scale_size_continuous(name="OM %")+
  scale_color_continuous(name="pH")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.box = "horizontal")
plot_grid(MCC+guides(color="none", size="none"), nonMCC+guides(color="none", size="none"),
          get_legend(MCC), 
          rel_widths = c(1,1,0.5),
          ncol=3, align="hv", axis="t", labels=c("A", "B"))
ggsave("github/figures/SuppFigure2.pdf", height=2.75, width=6.5)

# Supplemental Figure 3
# See ddPCRcounts script 

# Supplemental Figure 4
p16S=ggplot(DF, aes(x=chl.a, y=copies_16S_per_g))+geom_point()+
  ylab("Bacterial + archaeal abundance\n(copies 16S per g sediment)")+
  xlab(expression(paste("chl a (", mu,"g per L)")))+
  stat_smooth(method="lm", formula=y~x)+
  annotate("text", x=10, y=8.5e7, label="R^2 == 0.24", parse=T)+
  theme_bw()+
  theme(panel.grid=element_blank())
pmcrA=ggplot(DF, aes(x=chl.a, y=copies_mcrA_per_g))+geom_point()+
  ylab("Methanogen abundance\n(copies mcrA per g sediment)")+
  xlab(expression(paste("chl a (", mu,"g per L)")))+
  stat_smooth(method="lm", formula=y~x)+
  annotate("text", x=10, y=1.2e7, label="R^2 == 0.19", parse=T)+
  theme_bw()+
  theme(panel.grid=element_blank())
plot_grid(p16S,pmcrA, ncol=2, labels=c("A", "B"))
ggsave("github/figures/SuppFigure4.pdf", height=3.25, width=6.5)
