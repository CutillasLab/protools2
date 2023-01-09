
pca.plot <- function(df,df.design,colorfactor="",shapefactor=""){

  df <- dplyr::select_if(df,is.numeric)

  res.pca <- prcomp(df, scale = F)

  df.pca <- data.frame(res.pca$rotation)
  df.pca <- data.frame(heading=rownames(df.pca),df.pca)

  df.pca <- merge.data.frame(df.pca,df.design,by="heading")
  ss <- summary(res.pca)
  df.sum <- data.frame(ss$importance)
  pc1 <- df.sum$PC1[2]
  pc2 <- df.sum$PC2[2]
  pc3 <- df.sum$PC3[2]

  f <- round(nrow(df.pca)/20)

  plot1 <- ggplot(df.pca,aes_string(x="PC1",y="PC2",color=colorfactor,shape=shapefactor))+
    geom_point()+
    scale_color_manual(values = rep( mycolors()$C23,f))+
    scale_shape_manual(values = rep(c(1:20),f))+
    theme_bw()+
    labs(title = "", x=paste("PC1,", round(pc1,digits = 2)*100,"%"),
         y=paste("PC2,", round(pc2,digits = 2)*100,"%"))

  plot2 <- ggplot(df.pca,aes_string(x="PC2",y="PC3",color=colorfactor,shape=shapefactor))+
    geom_point()+
    scale_color_manual(values =  rep( mycolors()$C23,f))+
    scale_shape_manual(values = rep(c(1:20),f))+
    theme_bw()+
    labs(title = "", x=paste("PC2,", round(pc2,digits = 2)*100,"%"),
         y=paste("PC3,", round(pc3,digits = 2)*100,"%"))

  return(cowplot::plot_grid( plot1,plot2))
}
