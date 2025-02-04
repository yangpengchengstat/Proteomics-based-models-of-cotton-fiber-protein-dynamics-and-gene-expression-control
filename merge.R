meannorm_informative <- NULL
for (frac in fraction_names){
  meannorm_informative <- rbind(meannorm_informative, mean_of_normalized_informative[[frac]])
}
Data_meannorm <- as.matrix(meannorm_informative)


merge_informative_codes <- function(som_grid_x, som_grid_y, Data,codes_informative){
  
  x_limit <- som_grid_x
  y_limit <- som_grid_y
  grid_size <- paste0(paste0("(SOM grid size of ", x_limit*y_limit), ")")
  set.seed(123)
  
  system.time(SOM_wcc <- supersom(Data, 
                                  grid=somgrid(xdim = x_limit, ydim = y_limit, topo = c("rectangular"), 
                                               neighbourhood.fct = "gaussian", toroidal = F),  
                                  dist.fcts = c("WCCd0"), 
                                  rlen=500,
                                  normalizeDataLayers = FALSE))
  prot_id <- rownames(Data)
  label_of_sample <- SOM_wcc$unit.classif
  cluster_ID <- levels(factor(label_of_sample))
  
  sample_in_each_unit <- list()
  for (i in c(1:(x_limit*y_limit))){
    postion_index <- c(label_of_sample==i)
    name <- paste0('V',as.character(i), sep="")
    if (length(label_of_sample[postion_index])>0){
      sample_in_each_unit[[name]] <- prot_id[postion_index]
    }else{
      sample_in_each_unit[[name]] <- NULL
    }
  }
  clustering_inform_code <- NULL
  for (c in names(sample_in_each_unit)){
    df_c <- codes_informative[sample_in_each_unit[[c]], ]
    df_c[["cluster_ID"]] <- rep(c, dim(df_c)[1])
    clustering_inform_code <- rbind(clustering_inform_code, df_c)
  }
  write.csv(clustering_inform_code, 
            paste(paste(paste0("Table SOM clustering of informative mean(norm) of APO S200 and P200",grid_size), Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
  codes_merge <- SOM_wcc$codes[[1]]
  write.csv(codes_merge, 
            paste(paste(paste0("Table SOM Codes of merging informative mean(norm) of APO S200 and P200",grid_size), Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
  plot_data <- list()
  for (v in names(sample_in_each_unit)){
    df <- codes_informative[sample_in_each_unit[[v]], ]
    value <- c()
    type_ <- c()
    DPA <- c()
    names <- c()
    for (r in rownames(df)){
      zx= t(df[r,])
      cc=zx[,colnames(zx)]
      value <- c(value, cc)
      type_ <- c(type_, rep(strsplit(r, split = "_")[[1]][1], 20))
      names <- c(names, rep(r, 20))
      DPA <- c(DPA, 5:24)
    }
    value <- c(value, SOM_wcc$codes[[1]][v,])
    type_ <- c(type_, rep(v, 20))
    names <- c(names, rep(v, 20))
    DPA <- c(DPA, 5:24)
    plot_data[[v]] <- data.frame(value=value, DPA = DPA, Fraction=type_, Code=names)
  }
  
  mx <- c()
  mn <- c()
  for (i in names(plot_data)){
    mx<-c(mx,max(plot_data[[i]]$value))
    mn<-c(mn,min(plot_data[[i]]$value))
  }
  max_y <- max(mx)
  z = c("APO", "P200","S200","V")
  color_set <- c("#ed5151","#149ece","#a7c636","black")
  pdf(paste(paste(paste0(paste0("Merge informative ", "mean(norm)", " plot of APO P200 S200"),grid_size), Sys.Date(), sep=' '),
            "pdf", sep='.'), width=12,height=8)
  plot(SOM_wcc, type = "codes", main="Codes")
  plot(SOM_wcc, type = "mapping",
       pch = 1, main = "Count plot", keepMargins = TRUE)
  for (V in names(plot_data)){
    df_t <- plot_data[[V]]
    colors_t <- levels(factor((df_t$Fraction)))
    if (substr(colors_t[length(colors_t)],1,1)=="V"){
      colors_t <- c(colors_t[-length(colors_t)], "V")
    }
    colors <- color_set[z %in% colors_t]
    lineso <- c(2:(length(levels(factor((df_t$Code))))), 1)
    plot_p <- ggplot(df_t, aes(x=DPA, y=value,linetype=Code,color=Fraction,
                               group=interaction(Code, Fraction)))+
      geom_line()+
      scale_linetype_manual(values = lineso)+
      scale_colour_manual(values = colors)+
      labs(title=V,subtitle = paste0(paste0("{",paste(sample_in_each_unit[[V]],collapse = ";")),"}"))+
      scale_x_continuous(breaks=seq(5,24,1),limits = c(5,24))+
      scale_y_continuous(breaks=seq(0,1,0.1),limits = c(0,1))+
      theme(plot.title=element_text(hjust=0.5, vjust=12, margin=margin(t=40,b=-30)), title =element_text(size=14),
            
            axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), 
            axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))+
      theme(legend.title = element_text(size = 8,face = 8),legend.text=element_text(size=8),legend.position = "bottom")
    print(plot_p)
  }
  dev.off()  
  
}
#################################################
x_lim <- 3
y_lim <- 3
merge_informative_codes(som_grid_x=x_lim, som_grid_y=y_lim, Data=Data_meannorm, codes_informative=meannorm_informative)

x_lim <- 5
for (y_lim in c(4,6,8)){
  merge_informative_codes(som_grid_x=x_lim, som_grid_y=y_lim, Data=Data_meannorm, codes_informative=meannorm_informative)
}
