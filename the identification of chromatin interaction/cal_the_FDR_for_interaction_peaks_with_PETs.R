#aa=read.table("interaction_peak_with_PETs_summary");
#FDR=1-phyper(aa$V1,aa$V9,(N_peak-aa$V9),aa$V5);
#write.table(mydata,file="interaction_peak_with_PETs_summary_with_FDR",sep=" ",quote=FALSE,append=FALSE,na="NA");

aa=read.table("./ear_H3K27ac_PETs_in_peak_both_end");
bb=read.table("./ear_H3K27ac_interaction_peak_with_PETs_summary");
#N_peak=sum(as.numeric(bb$V5))+sum(as.numeric(bb$V9));
N_peak=sum(as.numeric(aa$V1));
p.val=1-phyper(bb$V1-1,bb$V9,(N_peak-bb$V9),bb$V5);
p.val <- as.numeric(p.val)
p.val[p.val > 1] <- 1
p.val[p.val < 0] <- 0
bb <- cbind(bb, p.val);
cat("P-value cbind Done!\n\n");
 p_adjust=p.adjust(p.val, method = "BH")
#library(qvalue)
#q.val <- qvalue(p.val)
#cat("Q-value calculated Done!\n\n");

#mydata=cbind(aa,q.val$qvalues);
mydata=cbind(bb,p_adjust);
write.table(mydata,file="ear_H3K27ac_interaction_peak_with_PETs_summary_with_FDR",sep="\t",quote=FALSE,append=FALSE,na="NA");

