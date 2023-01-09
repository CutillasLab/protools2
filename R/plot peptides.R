

barplot.top.peptides <- function(df.up,df.do,graph.heading=""){

  df.up$effect <- "Increased"
  df.do$effect <- "Decreased"

  df.up.down <-  na.omit(rbind.data.frame( df.up[1:25,],df.do[1:25,]))


  a <- df.up.down$sites
  a <- ifelse(nchar(a) > 20, paste0(strtrim(a, 20), '...'), a)
  df.up.down$sites <- a


  barplot.top.up <- ggplot(df.up.down,aes(x=reorder(sites,fold),y=fold))+
    geom_bar(stat = "identity", aes(fill=effect))+
    scale_fill_manual(values = c("skyblue","red"))+
    #coord_flip()+
    #geom_text(aes(label=paste0("q-val = ",round(qvalue,digits = 3))),hjust=0)+
    theme_bw()+
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(title = paste0("Top increased or decreased phosphosites in ", graph.heading),
         y="log2(fold)",
         x="",
         subtitle = "All values q<0.05")
}
