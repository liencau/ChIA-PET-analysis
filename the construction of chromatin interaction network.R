library(igraph)
pdf(file="ear.pdf")
par(mfrow=c(1,1),mai=c(0.3,0.3,0.3,0.3))
#g1=read.table("/NAS4/lien/Chia-pet/seedling_H3K27ac/only_use_Rep1/get_interaction_in_diff_peaks_and_one_peak/interaction_between_peak_within_interaction_3pets_for_network_format_lien_add_E.txt",header=T)
g1=read.table("ear_interaction_peak_3pets_intrachr_merged_anontataion_for_igraph_formated_filter.txt",header=T)
g2 = graph.data.frame(d = g1,directed = F)
com = walktrap.community(g2, steps = 1)
MEM=membership(com)
for(i in 1:max(as.vector(MEM)))
        {
        a=MEM[MEM==i]
        b=names(a)
        g3=induced_subgraph(g2, b)
	print(paste("coreset",i,sep=""))
	print (E(g3)[[]])
	#next
	for(j in 1:10)
		{
		V(g3)[degree(g3)>=j]$size=j*1
		V(g3)[degree(g3)>=j]$cex=j*0.05
		}
	for(j in 1:length(names(V(g3))))
                {
                if(gregexpr("E",names(V(g3))[j])[[1]][1]==1)
                        {
                        V(g3)[j]$color="red"
                        }
                if(gregexpr("E",names(V(g3))[j])[[1]][1]==-1)
                        {
                        V(g3)[j]$color="blue"
                        }
                }
			
	plot(g3,layout=layout_nicely,edge.width=2,vertex.color=V(g3)$color,vertex.color=V(g3)$color,edge.arrow.size=0.05,vertex.size=V(g3)$size,edge.color=E(g3)$type,vertex.label.cex=V(g3)$cex,vertex.label.color="white",vertex.frame.color=V(g3)$color,main=paste("Core",i,sep=""))

	}
dev.off()
