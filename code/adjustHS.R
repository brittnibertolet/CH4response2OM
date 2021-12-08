#function for adjusting headspace concentration to account for dilution in sampling
#based on methods in West et al. 2012

adjustHS=function(df, lakeIDs){
  outAdjustedConc=data.frame()
  treatments=c("control", "algal")
  for(k in 1:length(lakeIDs)){
    temp=df[df$lakeID==lakeIDs[k],]
    for(i in 1:length(treatments)){
      temp2=temp[temp$depthTop==treatments[i],]
      for(j in 1:3){
        temp3=temp2[temp2$replicate==j,]
        
        # adjusting for N2 added for repeated sampling
        # for UNDERC 2012 headspace vol=219-50-50 ml and 10 ml of sample was taken and replaced with 10ml of N2 each time
        HSvol=(250-50-50)/1000
        sampVol=0.01
        N2pctadd=sampVol/HSvol
        
        if(length(unique(temp3$subsampleDateTime))>2){
          forRate=data.frame(lakeID=lakeIDs[k], treatment=treatments[i],
                             sample_times=strptime(temp3$subsampleDateTime,format="%m/%d/%y %H:%M"),
                             conc=temp3$CH4original_umolL, 
                             rep=temp3$replicate)
          forRate$incub_days=as.numeric(forRate$sample_times-forRate$sample_times[1])/(60*60*24)
          
          afterMeasure=forRate$conc*(1-N2pctadd)
          delta=forRate$conc[2:length(afterMeasure)]-afterMeasure[1:(length(afterMeasure)-1)]
          
          forRate$adj_conc=c(forRate$conc[1],forRate$conc[1:(length(forRate$conc)-1)]+delta)

          
          # Bind to outputs
          outAdjustedConc=rbind(outAdjustedConc, forRate)
        }
      }
    }
  }
  return(list(outAdjustedConc))
}
