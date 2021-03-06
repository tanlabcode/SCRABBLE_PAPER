# This is the main R file to generate the results related our manuscript
# figure 2 and the supplementary figures related to figure 2.
# Please contact Tao Peng: pengt@email.chop.edu if you have any questions 
# about the scripts or data

# load the libraries 

library(dplyr)
library(vsn)
library(Seurat)
library(scater)
library(edgeR)
library(gridExtra)
library(R.matlab)
library(cowplot)
library(biomaRt)
library(data.table)
library(lattice)
library(scImpute)
library(SCRABBLE)
library(VennDiagram)
library(Rtsne)
library(DT)
library(ggpubr)
library(ggsignif)
library(scatterplot3d)
library(ggplot2)
library(reshape2)
library(ggfortify)
library(refGenome)
library(pheatmap)
library(RColorBrewer)
library(dendsort)
library(entropy)
library(DrImpute)
library(splatter)
library(RColorBrewer)
library(mcriPalettes)
library(plotly)
library(factoextra)
library(cluster)
library(NbClust)
library(fpc)
library(class)
library(VIPER)
library(SC3)

source("analysis_library.R")


# the following script is to generate the simulation data. Here we use 
# HPC to generate the simulation which could reduce the running time
for(dropout_index in c(1:3) ){
  
  for(seed_value in c(1:100)){
    
    generate_save_data(dropout_index, seed_value)
    
  }
}

# The reader could skip the above script and we have deposit the genereate 
# data in the current folder /simulation_data/

# the following codes are used to do the imputation using different methods

# the following script is to impute data. Here we use 
# HPC to impute the data using drimpute which could reduce the running time
for(dropout_index in c(1:3) ){
  
  for(seed_value in c(1:100)){
    
    run_drimpute(dropout_index, seed_value)
    
  }
  
}

# the following script is to impute data. Here we use 
# HPC to impute the data using scimpute which could reduce the running time
for(dropout_index in c(1:3) ){
  
  for(seed_value in c(1:100)){
    
    run_scimpute(dropout_index, seed_value)
    
  }
  
}


# the following script is to impute data. Here we use 
# HPC to impute the data using scimpute which could reduce the running time
for(dropout_index in c(1:3) ){
  
  for(seed_value in c(1:100)){
    
    run_viper(dropout_index, seed_value)
    
  }
  
}

# the following script is to impute data. Here we use 
# HPC to impute the data using scrabble which could reduce the running time
for(dropout_index in c(1:3) ){
  
  for(seed_value in c(1:100)){
    
    run_scrabble(dropout_index, seed_value)
    
  }
}

# the following script is to calculate the errors. Here we use 
# HPC to impute the data using scrabble which could reduce the running time
for(dropout_index in c(1:3) ){
  
  for(seed_value in c(1:100)){
    
    result <- calculate_error_splatter(drop_index, seed_value)
    
    dir.create(file.path('error_data'), showWarnings = FALSE)
    
    saveRDS(result,
            file = paste0("error_data/error_",drop_index,"_",seed_value,".rds")
    )
    
  }
  
}


# the following script is to calculate the errors. Here we use 
# HPC to impute the data using scrabble which could reduce the running time
for(dropout_index in c(1:3) ){
  
  for(seed_value in c(1:100)){
    
    cal_cell_distribution(dropout_index, seed_value)
    
  }
  
}


# the following script is to calculate the errors. Here we use 
# HPC to impute the data using scrabble which could reduce the running time
for(dropout_index in c(1:3) ){
  
  for(seed_value in c(1:100)){
    
    cal_gene_distribution(dropout_index, seed_value)
    
  }
  
}

# Gather the errors
# -----------------
error_list <- list()

error_cell_list <- list()

error_gene_list <- list()

# gather the error from the data from different dropout rates
for(i in c(1:3)){
  error_matrix <- c()
  
  error_cell_matrix <- c()
  
  error_gene_matrix <- c()
  
  for(j in c(1:100)){
    
    tmp <- readRDS(file = paste0("error_data/error_",i,"_",j,".rds"))
    
    error_matrix <- cbind(error_matrix,as.matrix(tmp$error))
    
    error_cell_matrix <- cbind(error_cell_matrix,as.matrix(tmp$error_cell))
    
    error_gene_matrix <- cbind(error_gene_matrix,as.matrix(tmp$error_gene))
    
  }
  
  error_list[[i]] <- error_matrix
  
  error_cell_list[[i]] <- error_cell_matrix
  
  error_gene_list[[i]] <- error_gene_matrix
  
}

# save the errors
saveRDS(error_list,file = "error_all.rds")

saveRDS(error_cell_list,file = "error_all_cell.rds")

saveRDS(error_gene_list,file = "error_all_gene.rds")
# ------

# Plot the boxplots in figure 2 and the supplementary figures related to Figure 2
# ----------------------------------------------------------------------------
# load the error data
error_list <- readRDS(file = "error_all.rds")

error_cell_list <- readRDS(file = "error_all_cell.rds")

error_gene_list <- readRDS(file = "error_all_gene.rds")


# define the list for boxplot comparisons
p <- list()

# Dropout rate: 71%
p[[1]] <- plot_comparison(error_list[[1]], "Error", 1400, 140)

p[[2]] <- plot_comparison(error_cell_list[[1]], "Correlation", 1, 0.1)

p[[3]] <- plot_comparison(error_gene_list[[1]], "Correlation", 0.4, 0.04)

# Dropout rate: 83%
p[[4]] <- plot_comparison(error_list[[2]], "Error", 1800, 180)

p[[5]] <- plot_comparison(error_cell_list[[2]], "Correlation", 1, 0.1)

p[[6]] <- plot_comparison(error_gene_list[[2]], "Correlation", 0.3, 0.03)

# Dropout rate:87%
p[[7]] <- plot_comparison(error_list[[3]], "Error",1800,180)

p[[8]] <- plot_comparison(error_cell_list[[3]], "Correlation", 1, 0.1)

p[[9]] <- plot_comparison(error_gene_list[[3]], "Correlation", 0.2, 0.02)

# save the PDF files
main <- grid.arrange(grobs = p,ncol = 3)
ggsave(filename="Figure_error.pdf", 
       plot = main, 
       width = 18, 
       height = 12)
# ----------------------------------------------------------------------------


# plot the mean-variance figures in Figure 2 and the supplementary figures realted to Figure 2
# ----------------------------------------------------------------------------
for(drop_index in c(1:3)){
  
  seed_value <- 10
  
  p <- plot_meansd(drop_index, seed_value)
  
  ggsave(filename=paste0("Figures_mean_variance_",drop_index,".pdf"),
         plot = p,
         width = 24,
         height = 3)
}
# ----------------------------------------------------------------------------


# plot the mean-variance figures in Figure 2 and the supplementary figures realted to Figure 2
# ----------------------------------------------------------------------------
for(drop_index in c(1:3)){
  
  seed_value <- 10
  
  p <- plot_comparison_tsne(drop_index,
                            seed_value,
                            50,100)
  
  ggsave(filename=paste0("Figures_tsne_",drop_index,".pdf"),
         plot = p,
         width = 24,
         height = 3)
}
# -------------------------------------------------------------------------

# plot the mean-variance figures in Figure 2 and the supplementary figures realted to Figure 2
# ----------------------------------------------------------------------------
for(drop_index in c(1:3)){
  
  seed_value <- 10
  
  p <- plot_ma(drop_index, seed_value)
  
  ggsave(filename=paste0("Figure_ma_",drop_index,".pdf"), 
         plot = p, 
         width = 24, 
         height = 3)
}
# -------------------------------------------------------------------------


# plot the distribution of cell-cell correlation figures in Figure 5 and the supplementary figures realted to Figure 5
# ----------------------------------------------------------------------------
for(drop_index in c(1:3)){
  
  seed_value <- 10
  
  p <- plot_cell_distribution(drop_index, seed_value)
  
  ggsave(filename=paste0("Figure_cell_distribution_",drop_index,".pdf"), 
         plot = p, 
         width = 35, 
         height = 12)
}
# -------------------------------------------------------------------------


# plot the distribution of gene-gene correlation figures in Figure 5 and the supplementary figures realted to Figure 5
# ----------------------------------------------------------------------------
for(drop_index in c(1:3)){
  
  seed_value <- 10
  
  p <- plot_gene_distribution(drop_index, seed_value)
  
  ggsave(filename=paste0("Figure_gene_distribution_",drop_index,".pdf"), 
         plot = p, 
         width = 35, 
         height = 12)
}
# -------------------------------------------------------------------------

# plot the distribution of gene-gene correlation figures in Figure 5 and the supplementary figures realted to Figure 5
# ----------------------------------------------------------------------------

methods <- c("True Data", "Dropout Data", "DrImpute", "scImpute", "MAGIC", "VIPER", "SCRABBLE")

for(drop_index in c(1:3)){
  
  data_cor <- get_cor_data(drop_index, 10)
  
  data_cell <- data_cor[[1]]
  
  data_gene <- data_cor[[2]]
  
  p <- list()
  
  k <- 1
  
  for(i in c(1:length(methods))){
    
    data_plot <- data.frame(x = data_cor_vector(data_cell[[i]]))
    
    p[[k]] <- ggplot(data_plot, aes(x=x)) + 
              geom_histogram(color="darkblue", fill="lightblue") + 
              ggtitle(paste0("Cell: ", methods[i]))
    
    k <- k + 1
    
  }
  
  for(i in c(1:length(methods))){
    
    data_plot <- data.frame(x = data_cor_vector(data_gene[[i]]))
    
    p[[k]] <- ggplot(data_plot, aes(x=x)) + 
              geom_histogram(color="darkblue", fill="lightblue") + 
              ggtitle(paste0("Gene: ", methods[i]))
    
    k <- k + 1
    
  }
  
  main <- grid.arrange(grobs = p,ncol = 8)
  
  ggsave(filename= paste0("Figure_splatter_histogram_test_all_",drop_index,".pdf"), 
         plot = main, 
         width = 35, 
         height = 8)
}
# -------------------------------------------------------------------------

# Assemble the \delta KS values between cell clusters
# -------------------------------------------------------------------------
data_tmp <- list()

k <- 1

for(i in c(1:3)){
  
  tmp1 <- c()
  
  for(j in c(1:100)){
    
    tmp <- readRDS(file = paste0("data_cell_distribution/data_",i,"_",j,".rds"))
    
    tmp1 <- rbind(tmp1, tmp)
    
  }
  
  data_tmp[[k]] <- tmp1
  
  k <- k + 1
  
}

for(i in c(1:3)){
  
  dataV0 <- t(data_tmp[[i]])
  
  dataV1 <- data.frame(as.vector(t(dataV0)))
  
  ylim_value <- 0.85
  
  h_ylim <- 0.1
  
  # calculate the dropout rate
  N <- dim(dataV1)[1]    
  
  dataV1$group <- as.factor(rep(c(1:6), each = N/6))
  
  colnames(dataV1) <- c('y','group')
  
  my_comparisons <- list( c("1", "6"), c("2", "6"), c("3", "6"), c("4", "6"), c("5", "6"))
  
  pval <- compare_means(y ~ group,data = dataV1, method = "t.test", ref.group = "6", paired = TRUE)
  
  pp <- ggboxplot(dataV1, x = "group", y = "y", fill = "group",
                  palette = c("#00AFBB","#0000CD", "#E7B800", "#FC4E07",  "#ef8a62", "#6ebb00")) +
    stat_boxplot(geom = "errorbar", width = 0.3) +
    ylim(c(0,ylim_value + 4.5*h_ylim)) + 
    theme_bw() +
    geom_signif(comparisons = my_comparisons, 
                annotations = formatC(pval$p, format = "e", digits = 2),
                tip_length = 0.01,
                y_position = c(ylim_value + 4*h_ylim, ylim_value + 3*h_ylim, ylim_value + 2*h_ylim, ylim_value + h_ylim, ylim_value)) +
    theme(text=element_text(size=14)) +
    xlab("Method") + 
    ylab("Error") + 
    ggtitle(paste0("Dropout Rate: ",i,"%")) +
    scale_fill_manual( values = c("#00AFBB","#0000CD", "#E7B800", "#FC4E07",  "#ef8a62", "#6ebb00"),
                       name="Method",
                       breaks=c("1", "2", "3", "4", "5", "6"),
                       labels=c("Dropout", "DrImpute", "scImpute", "MAGIC", "VIPER", "SCRABBLE")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  ggsave(filename= paste0("Figure_splatter_boxplot_gene_",i,".pdf"), 
         plot = pp, 
         width = 6, 
         height = 4)
  
}
# -------------------------------------------------------------------------


# Assemble the \delta KS values between marker genes
# -------------------------------------------------------------------------
data_tmp <- list()

k <- 1

for(i in c(1:3)){
  
  tmp1 <- c()
  
  for(j in c(1:100)){
    
    tmp <- readRDS(file = paste0("data_cell_distribution/data_",i,"_",j,".rds"))
    
    tmp1 <- rbind(tmp1, tmp)
    
  }
  
  data_tmp[[k]] <- tmp1
  
  k <- k + 1
  
}

for(i in c(1:3)){
  
  dataV0 <- t(data_tmp[[i]])
  
  dataV1 <- data.frame(as.vector(t(dataV0)))
  
  ylim_value <- 0.85
  
  h_ylim <- 0.1
  
  # calculate the dropout rate
  N <- dim(dataV1)[1]   
  
  dataV1$group <- as.factor(rep(c(1:6), each = N/6))
  
  colnames(dataV1) <- c('y','group')
  
  my_comparisons <- list( c("1", "6"), c("2", "6"), c("3", "6"), c("4", "6"), c("5", "6"))
  
  pval <- compare_means(y ~ group,data = dataV1, method = "t.test", ref.group = "6", paired = TRUE)
  
  pp <- ggboxplot(dataV1, x = "group", y = "y", fill = "group",
                  palette = c("#00AFBB","#0000CD", "#E7B800", "#FC4E07",  "#ef8a62", "#6ebb00")) +
    stat_boxplot(geom = "errorbar", width = 0.3) +
    ylim(c(0,ylim_value + 4.5*h_ylim)) + 
    theme_bw() +
    geom_signif(comparisons = my_comparisons, 
                annotations = formatC(pval$p, format = "e", digits = 2),
                tip_length = 0.01,
                y_position = c(ylim_value + 4*h_ylim, ylim_value + 3*h_ylim, ylim_value + 2*h_ylim, ylim_value + h_ylim, ylim_value)) +
    theme(text=element_text(size=14)) +
    xlab("Method") + 
    ylab("Error") + 
    ggtitle(paste0("Dropout Rate: ",i,"%")) +
    scale_fill_manual( values = c("#00AFBB","#0000CD", "#E7B800", "#FC4E07",  "#ef8a62", "#6ebb00"),
                       name="Method",
                       breaks=c("1", "2", "3", "4", "5", "6"),
                       labels=c("Dropout", "DrImpute", "scImpute", "MAGIC", "VIPER", "SCRABBLE")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  ggsave(filename= paste0("Figure_splatter_boxplot_gene_",i,".pdf"), 
         plot = pp, 
         width = 6, 
         height = 4)
  
}
# -------------------------------------------------------------------------


