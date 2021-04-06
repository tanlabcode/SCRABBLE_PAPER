######## Data Generation ##############
# the function generating the simulation data using bioconductor package Splatter
generate_simulation_splatter <- function(dropout_index, seed_value, nGenes = 800){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  #    seed_value: the random seed
  #        nGenes: the number of genes in the simulation data. The default is 800
  
  # Set up the parameters
  params = newSplatParams()
  
  params = setParams(params, list(batchCells = 1000,
                                  nGenes = nGenes,
                                  group.prob = c(0.20, 0.35, 0.45),
                                  de.prob = c(0.045, 0.045, 0.045),
                                  de.facLoc = 0.1,
                                  de.facScale = 0.4)
  )
  
  # Set up the vector of dropout.mid
  dropout_mid = c(4, 5, 5.5)
  
  # determine if it is a good parameter
  if(dropout_index > length(dropout_mid)){
    
    stop(
      paste0('The dropout_index shold not be greater than ', 
             length(dropout_mid), 
             ' . Please input a proper one.\n')
    )
    
  }
  
  # Generate the simulation data using Splatter package
  sim = splatSimulateGroups(params,
                            dropout.type = "experiment",
                            dropout.shape = -1,
                            dropout.mid = dropout_mid[dropout_index],
                            seed = seed_value)
  
  # genereate the cpm levels of the true simulation data
  data_true = cpm(sim@assays@data$TrueCounts)
  
  data_dropout = Matrix::Matrix(data_true)
  
  # generate the dropout data based on the counts in sim
  data_dropout[counts(sim) == 0] = 0
  
  # calculate the dropout rate
  percentage_zeros = round(Matrix::nnzero(data_dropout == 0, na.counted = NA)/
                             (dim(data_dropout)[1]*dim(data_dropout)[2])*100)
  
  
  # generate the bulk RNAseq data
  data_bulk = data.frame(val = rowMeans(data_true))
  
  # define the data list for the simulation data
  # indcluding: data_true: true data
  #          data_dropout: dropout data
  #             data_bluk: bulk data
  #      percentage_zeros: dropout rate
  #                 group: the group label
  
  data = list()
  
  data$data_bulk = data_bulk
  
  data$data_dropout = data_dropout
  
  data$data_true = data_true
  
  data$percentage_zeros = percentage_zeros
  
  data$group = colData(sim)@listData$Group
  
  return(data)
}

# generate the simulation data and save the data
generate_save_data <- function(dropout_index, seed_value){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # generate the simulation data
  data_simulation = generate_simulation_splatter(dropout_index, seed_value)
  
  # generate the folder saving the simulation data
  dir.create(file.path('simulation_data'), showWarnings = FALSE)
  
  # save the data as RDS format
  saveRDS(data_simulation, 
          file = paste0('simulation_data/simulation_data_dropout_index_',
                        dropout_index, 
                        '_seed_', 
                        seed_value,
                        '.rds')
  )
  
}

######## DrImpute ##############
run_drimpute <- function(dropout_index, seed_value){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # load the raw data
  data <- readRDS(file = paste0('simulation_data/simulation_data_dropout_index_',
                                dropout_index,
                                '_seed_',
                                seed_value,
                                '.rds')
  )
  
  # build the folder saving the imputed data using DrImpute
  path <- "imputation_drimpute_data/"
  
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using DrImpute
  data_dropout <- as.matrix(data$data_dropout)
  
  exdata <- DrImpute(data_dropout)
  
  # write the data
  write.table(exdata,
              paste0(path, "drimpute_",dropout_index,"_",seed_value,".csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
}

######## scImpute ##############
run_scimpute <- function(dropout_index, seed_value){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  
  # load the data
  data = readRDS(file = paste0('simulation_data/simulation_data_dropout_index_',
                               dropout_index,
                               '_seed_',
                               seed_value,
                               '.rds')
  )
  
  # build the folder saving the imputed data using scimpute method
  path = "imputation_scimpute_data/"
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using scImpute
  data_dropout = data$data_dropout
  
  write.table(data_dropout,paste0(path, 
                                  "dropout_scimpute_",
                                  dropout_index,"_",
                                  seed_value,".csv"),
              sep=',',
              row.names = TRUE,
              col.names = TRUE
  )
  
  file.remove(paste0(path, "scimpute_", dropout_index, "_", seed_value,"_*"))
  
  # run scImpute
  scimpute(
    paste0(path, "dropout_scimpute_",dropout_index,"_",seed_value,".csv"),
    infile = "csv",           
    outfile = "csv",         
    out_dir = paste0(path, "scimpute_", dropout_index, "_", seed_value,"_"),
    drop_thre = 0.5,         
    Kcluster = 2,
    ncores = 2)       
  # 
  # clean the data
  data_dropout = read.table( file = paste0(path, "scimpute_",
                                           dropout_index, "_",
                                           seed_value,
                                           "_scimpute_count.csv") ,
                             header = TRUE, sep=",")
  
  data_dropout$X = NULL
  
  # save the data
  write.table(data_dropout,
              paste0(path, "data_imputation_scimpute_",dropout_index,"_",seed_value,".csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
}

######## MAGIC ##############
# here we use python script to run magic and we post the python script here.
# -----------------------------------------------------------------------------
# import sys
# import os
# import magic 
# import pandas as pd
# 
# cwd = os.getcwd()
# 
# if not os.path.exists(cwd+"/imputation_magic_data"):
#   os.makedirs(cwd+"/imputation_magic_data")
#
# X =pd.read_csv("simulation_data/simulation_data_dropout_index_"+str(dropout_value)+"_seed_"+str(seed_value)+".txt",sep = ' ', header=None)
#
# magic_operator = magic.MAGIC()
# X_magic = magic_operator.fit_transform(X.T)
#
# out_magic = X_magic.T
# out_magic.to_csv(cwd+"/imputation_magic_data/magic_"+str(dropout_value)+"_"+str(seed_value)+".csv", sep = '\t', header= None)
# -----------------------------------------------------------------------------


######## VIPER ##############
run_viper <- function(dropout_index, seed_value){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # load the raw data
  data <- readRDS(file = paste0('simulation_data/simulation_data_dropout_index_',
                                dropout_index,
                                '_seed_',
                                seed_value,
                                '.rds')
  )
  
  # build the folder saving the imputed data using DrImpute
  path <- "imputation_viper_data/"
  
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using VIPER
  data_dropout <- as.matrix(data$data_dropout)
  
  exdata <- VIPER(data_dropout, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1, 
                  report = FALSE, outdir = NULL, prefix = NULL)
  
  # write the data
  write.table(exdata$imputed,
              paste0(path, "viper_",dropout_index,"_",seed_value,".csv"),
              sep=',',
              row.names = F,
              col.names = F
  )
  
}

######## SCRABBLE ##############
run_scrabble <- function(dropout_index, seed_value){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  
  # load the data
  data = readRDS(file = paste0('simulation_data/simulation_data_dropout_index_',
                               dropout_index,
                               '_seed_',
                               seed_value,
                               '.rds')
  )
  
  
  
  path = "imputation_scrabble_data/"
  
  dir.create(file.path(path), showWarnings = FALSE)
  
  # impute the data using DrImpute
  data1 = list()
  
  data1[[1]] = data$data_dropout
  
  data1[[2]] = data$data_bulk
  
  # set up the parameters
  parameter = c(1, 1e-06, 1e-04)
  
  # run scrabble
  result = scrabble(data1,
                    parameter = parameter, 
                    nIter = 60,
                    error_out_threshold = 1e-7, 
                    nIter_inner = 100,
                    error_inner_threshold = 1e-5)
  
  # write the data
  write.table(result,
              paste0(path, "scrabble_",dropout_index,"_",seed_value,"_",j,".csv"),
              sep=',',
              row.names = F,
              col.names = F)
  
}


######## Data Analysis ##############
# Calculate the similarity of two datasets
calculate_similarity <- function(data1,data2){
  
  d = cor(c(data1[lower.tri(data1)]),c(data2[lower.tri(data2)]))
  
  return(d)
  
}


# Calculate the error between true data and imputed data
run_error <- function(dropout_index, seed_value){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  options( warn = -1 )
  # load the simulationd data
  data_simulation = readRDS(file = paste0('simulation_data/simulation_data_dropout_index_',
                                          dropout_index,
                                          '_seed_',
                                          seed_value,
                                          '.rds')
  )
  
  # get the index of the genes with nonzero means
  index = rowMeans(data_simulation$data_dropout) > 0
  
  # obtain the true data
  data_true = data_simulation$data_true
  
  data_true = data_true[index,]
  
  # cell-cell correlation of the true data
  data_true_cell = cor(as.matrix((data_true)))
  
  # gene-gene correlation of the true data
  data_true_gene = cor(t((data_true)), method = "pearson")
  
  # obtain the dropout data
  data_dropout = data_simulation$data_dropout
  
  data_dropout = data_dropout[index,]
  
  # cell-cell correlation of the dropout data
  data_dropout_cell = cor((data_dropout), method = "pearson")
  
  # gene-gene correlation of the dropout data
  data_dropout_gene = cor(t((data_dropout)), method = "pearson")
  
  # load imputed data from DrImpute
  data_drimpute = read.table( file = paste0("imputation_drimpute_data/drimpute_",
                                            dropout_index, "_",
                                            seed_value,
                                            ".csv") ,
                              header = FALSE, sep=","
  )
  
  data_drimpute = data_drimpute[index,]
  
  # cell-cell correlation of the DrImpute imputed data
  data_drimpute_cell = cor((data_drimpute), method = "pearson")
  
  # gene-gene correlation of the DrImpute imputed data
  data_drimpute_gene = cor(t((data_drimpute)), method = "pearson")
  
  # load imputed data from scImpute
  data_scimpute = read.table( file = paste0("imputation_scimpute_data/scimpute_",
                                            dropout_index, "_",
                                            seed_value,
                                            ".csv") ,
                              header = FALSE, sep=","
  )
  
  data_scimpute = data_scimpute[index,]
  
  # cell-cell correlation of the scImpute imputed data
  data_scimpute_cell = cor((data_scimpute), method = "pearson")
  
  # gene-gene correlation of the scImpute imputed data
  data_scimpute_gene = cor(t((data_scimpute)), method = "pearson")
  
  # load the MAGIC imputed data 
  data = read.csv(paste0("imputation_magic_data/magic_",dropout_index,"_",seed_value,".csv"),
                  header = FALSE,
                  sep = "\t")
  
  data$V1 = NULL
  
  data_magic = as.matrix(data)
  
  data_magic[data_magic < 0] = 0
  
  data_magic[is.nan(data_magic)] = 0
  
  data_magic = data_magic[index,]
  
  # cell-cell correlation of the MAGIC imputed data
  data_magic_cell = cor((data_magic), method = "pearson")
  
  # gene-gene correlation of the MAGIC imputed data
  data_magic_gene = cor(t((data_magic)), method = "pearson")
  
  
  # load imputed data from VIPER
  data_viper = read.table( file = paste0("imputation_viper_data/viper_",
                                            dropout_index, "_",
                                            seed_value,
                                            ".csv") ,
                              header = FALSE, sep=","
  )
  
  data_viper = data_viper[index,]
  
  # cell-cell correlation of the VIPER imputed data
  data_viper_cell = cor((data_viper), method = "pearson")
  
  # gene-gene correlation of the VIPER imputed data
  data_viper_gene = cor(t((data_viper)), method = "pearson")
  
  # load imputed data from scrabble
  data_scrabble = read.table( file = paste0("imputation_scrabble_data/scrabble_",
                                            dropout_index, "_",
                                            seed_value,
                                            ".csv") ,
                              header = FALSE, sep=","
  )
  
  data_scrabble = data_scrabble[index,]
  
  # cell-cell correlation
  data_scrabble_cell = cor(as.matrix(data_scrabble), method = "pearson")
  
  # gene-gene correlation
  data_scrabble_gene = cor(t((data_scrabble)), method = "pearson")
  
  # calulate the error between the imputed data and true data
  
  error = matrix(0, nrow = 7, ncol = 1)
  
  error[1] = norm(log10(data_dropout + 1) - log10(data_true + 1), type = "2")
  
  error[2] = norm(log10(data_drimpute + 1) - log10(data_true + 1), type = "2")
  
  error[3] = norm(log10(data_scimpute + 1) - log10(data_true + 1), type = "2")
  
  error[4] = norm(log10(data_magic + 1) - log10(data_true + 1), type = "2")
  
  error[5] = norm(log10(data_viper + 1) - log10(data_true + 1), type = "2")
  
  error[6] = norm(log10(data_scrabble + 1) - log10(data_true + 1), type = "2")
  
  error[7] = data_simulation$percentage_zeros
  
  # calulate the similarity between the cell-cell
  # correlation of the imputed data and the one of the true data
  
  error_cell = matrix(0, nrow = 7, ncol = 1)
  
  error_cell[1] = calculate_similarity(data_true_cell, data_dropout_cell)
  
  error_cell[2] = calculate_similarity(data_true_cell, data_drimpute_cell)
  
  error_cell[3] = calculate_similarity(data_true_cell, data_scimpute_cell)
  
  error_cell[4] = calculate_similarity(data_true_cell, data_magic_cell)
  
  error_cell[5] = calculate_similarity(data_true_cell, data_viper_cell)
  
  error_cell[6] = calculate_similarity(data_true_cell, data_scrabble_cell)
  
  error_cell[7] = data_simulation$percentage_zeros
  
  
  # calulate the similarity between the gene-gene
  # correlation of the imputed data and the one of the true da
  error_gene = matrix(0, nrow = 7, ncol = 1)
  
  error_gene[1] = calculate_similarity(data_true_gene, data_dropout_gene)
  
  error_gene[2] = calculate_similarity(data_true_gene, data_drimpute_gene)
  
  error_gene[3] = calculate_similarity(data_true_gene, data_scimpute_gene)
  
  error_gene[4] = calculate_similarity(data_true_gene, data_magic_gene)
  
  error_gene[5] = calculate_similarity(data_true_gene, data_viper_gene)
  
  error_gene[6] = calculate_similarity(data_true_gene, data_scrabble_gene)
  
  error_gene[7] = data_simulation$percentage_zeros
  
  # gather the errors as a list
  result = list()
  
  result$error = error
  
  result$error_cell = error_cell
  
  result$error_gene = error_gene
  
  return(result)
  
}

# Define the plots for errors

plot_comparison <- function(data, ylabels = "Error", ylim_value = 100, h_ylim = 10){
  
  # this function is used to plot the boxplot with the p-values
  # Parameter in the function
  # data : data is a matrix with 6 columns, the first five colums are the errors
  # the last column is the dropout rates
  # ylabels: the y labels shown on the graph
  # ylim_value : the smallest y value where the pvalue is shown
  # h_ylim: the difference between two pvalues shown in the graph. the best ones
  # is the 10%*ylim_vlaue
  
  # extract the data with the first five columns
  dataV0 = data[c(1:6),]
  
  dataV1 = data.frame(as.vector(t(dataV0)))
  
  # calculate the dropout rate
  dropout_rate = round(mean(data[7,]))
  
  # the number of data using for plotting the boxplot
  # dataV1 is built with two columns: y values and group labels
  N = dim(dataV1)[1]     
  
  dataV1$group = rep(c(1:6), each = N/6)
  
  colnames(dataV1) = c('y','group')
  
  # define the comparison lists
  my_comparisons <- list( c("1", "6"), c("2", "6"), c("3", "6"), c("4", "6"), c("5", "6"))
  
  # compare the errors using t-test
  pval <- compare_means(y ~ group,data = dataV1, method = "t.test", ref.group = "6", paired = TRUE)
  
  # plot the boxplot with pvalues for the comparisons
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
    ylab(ylabels) + 
    ggtitle(paste0("Dropout Rate: ",dropout_rate,"%")) +
    scale_fill_manual( values = c("#00AFBB","#0000CD", "#E7B800", "#FC4E07",  "#ef8a62", "#6ebb00"),
                       name="Method",
                       breaks=c("1", "2", "3", "4", "5", "6"),
                       labels=c("Dropout", "DrImpute", "scImpute", "MAGIC", "VIPER", "SCRABBLE")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  return(pp)
  
}


# plot the mean-variance function
plot_meansd <- function(dropout_index, seed_value){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  options( warn = -1 )
  # load the simulationd data
  data_simulation = readRDS(file = paste0('simulation_data/simulation_data_dropout_index_',
                                          dropout_index,
                                          '_seed_',
                                          seed_value,
                                          '.rds')
  )
  
  # true data
  data_true = data_simulation$data_true
  
  data_true = as.matrix(data_true)
  
  # raw data (dropout data)
  data_dropout = data_simulation$data_dropout
  
  data_dropout = as.matrix(data_dropout)
  
  # load imputed data from Drimpute
  data_drimpute = read.table( file = paste0("imputation_drimpute_data/drimpute_",
                                            dropout_index, "_",
                                            seed_value,
                                            ".csv") ,
                              header = FALSE, sep=","
  )
  
  data_drimpute = as.matrix(data_drimpute)
  
  # load imputed data from scImpute
  data_scimpute <- read.table( file = paste0("imputation_scimpute_data/scimpute_",
                                             dropout_index, "_",
                                             seed_value,
                                             ".csv") ,
                               header = FALSE, sep=","
  )
  
  data_scimpute = as.matrix(data_scimpute)
  
  # load the magic results 
  data = read.csv(paste0("imputation_magic_data/magic_",
                         dropout_index,"_",
                         seed_value,".csv"),
                  header = FALSE,
                  sep = "\t")
  
  data$V1 = NULL
  
  data_magic = as.matrix(data)
  
  # load imputed data from VIPER
  data_viper = read.table( file = paste0("imputation_viper_data/viper_",
                                            dropout_index, "_",
                                            seed_value,
                                            ".csv") ,
                              header = FALSE, sep=","
  )
  
  data_viper = as.matrix(data_viper)
  
  # load imputed data from scrabble
  data_scrabble = read.table( file = paste0("imputation_scrabble_data/scrabble_",
                                            dropout_index, "_",
                                            seed_value,
                                            ".csv") ,
                              header = FALSE, sep=","
  )
  
  data_scrabble = as.matrix(data_scrabble)
  
  ratio_zeros_raw <- 1 - Matrix::nnzero(data_dropout)/length(data_dropout)
  
  pl = list()
  
  p = meanSdPlot(log10(data_true + 1), ranks = FALSE)
  
  pl[[1]] = p$gg + ggtitle(paste0("True Data: ",round(100*ratio_zeros_true),"%")) 
  
  p = meanSdPlot(log10(data_dropout + 1), ranks = FALSE)
  
  pl[[2]] = p$gg + ggtitle(paste0("Raw Data: ",round(100*ratio_zeros_raw),"%"))
  
  p = meanSdPlot(log10(data_drimpute + 1), ranks = FALSE) 
  
  pl[[3]] = p$gg + ggtitle("DrImpute") 
  
  p = meanSdPlot(log10(data_scimpute + 1), ranks = FALSE)
  
  pl[[4]] = p$gg + ggtitle("scImpute")
  
  p = meanSdPlot(log10(data_magic + 1), ranks = FALSE)
  
  pl[[5]] = p$gg + ggtitle("MAGIC")
  
  p = meanSdPlot(log10(data_viper + 1), ranks = FALSE)
  
  pl[[6]] = p$gg + ggtitle("VIPER")
  
  p = meanSdPlot(log10(data_scrabble + 1), ranks = FALSE)
  
  pl[[7]] = p$gg + ggtitle("SCRABBLE")
  
  # combine the six plots as a whole one
  main = grid.arrange(grobs = pl,ncol = 7)
  
  return(main)
  
}

# Plot tsne function
plot_comparison_tsne <- function(dropout_index,
                                 seed_value, 
                                 initial_dims_value,
                                 perplexity_value){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  # initial_dims_value: the initial dimensions used in the tsne
  # perplexity_value: the perplexity used in the tnse
  
  data_simulation = readRDS(file = paste0('simulation_data/simulation_data_dropout_index_',
                                          dropout_index,
                                          '_seed_',
                                          seed_value,
                                          '.rds')
  )
  
  # true data
  data_true = data_simulation$data_true
  
  # raw data (dropout data)
  data_dropout = data_simulation$data_dropout
  
  # load imputed data from Drimpute
  data_drimpute = read.table( file = paste0("imputation_drimpute_data/drimpute_",
                                            dropout_index, "_",
                                            seed_value,
                                            ".csv") ,
                              header = FALSE, sep=","
  )
  
  # load imputed data from scImpute
  data_scimpute = read.table( file = paste0("imputation_scimpute_data/scimpute_",
                                            dropout_index, "_",
                                            seed_value,
                                            ".csv") ,
                              header = FALSE, sep=","
  )
  
  # load the magic results 
  data = read.csv(paste0("imputation_magic_data/magic_",
                         dropout_index,"_",
                         seed_value,".csv"),
                  header = FALSE,
                  sep = "\t")
  
  data$V1 = NULL
  
  data_magic = as.matrix(data)
  
  # load imputed data from VIPER
  data_viper = read.table( file = paste0("imputation_viper_data/viper_",
                                         dropout_index, "_",
                                         seed_value,
                                         ".csv") ,
                           header = FALSE, sep=","
  )
  
  
  # load imputed data from scrabble
  data_scrabble = read.table( file = paste0("imputation_scrabble_data/scrabble_",
                                            dropout_index, "_",
                                            seed_value,
                                            ".csv") ,
                              header = FALSE, sep=","
  )
  
  # define the plot list
  pl = list()
  
  set.seed(1) # Set a seed if you want reproducible results
  
  # calculate the tsne of true data
  true_tsne = Rtsne(t(as.matrix(data_true)), 
                    initial_dims = initial_dims_value, 
                    perplexity = perplexity_value)
  
  
  pl[[1]] = plot_pca_singlecell(true_tsne$Y,data_simulation$group)
  
  set.seed(1)
  
  # calculate the tsne of dropout data
  dropout_tsne = Rtsne(t(as.matrix(cpm(data_dropout))), 
                       initial_dims = initial_dims_value, 
                       perplexity = perplexity_value)
  
  pl[[2]] = plot_pca_singlecell(dropout_tsne$Y,data_simulation$group)
  
  set.seed(1)
  
  # calculate the tsne of drimpute imputed data
  drimpute_tsne = Rtsne(t(as.matrix(data_drimpute)), 
                        initial_dims = initial_dims_value, 
                        perplexity = perplexity_value)
  
  pl[[3]] = plot_pca_singlecell(drimpute_tsne$Y,data_simulation$group)
  
  set.seed(1)
  
  # calculate the tsne of scimpute imputed data
  scimpute_tsne = Rtsne(t(as.matrix(data_scimpute)), 
                        initial_dims = initial_dims_value, 
                        perplexity = perplexity_value)
  
  pl[[4]] = plot_pca_singlecell(scimpute_tsne$Y,data_simulation$group)
  
  set.seed(1)
  
  # calculate the tsne of magic imputed data
  magic_tsne = Rtsne(t(as.matrix(data_magic)), 
                     initial_dims = initial_dims_value, 
                     perplexity = perplexity_value)
  
  pl[[5]] = plot_pca_singlecell(magic_tsne$Y,data_simulation$group)
  
  set.seed(1)
  
  # calculate the tsne of viper imputed data
  viper_tsne = Rtsne(t(as.matrix(data_viper)), 
                        initial_dims = initial_dims_value, 
                        perplexity = perplexity_value)
  
  pl[[6]] = plot_pca_singlecell(viper_tsne$Y,data_simulation$group)
  
  set.seed(1)
  
  # calculate the tsne of scrabble imputed data
  scrabble_tsne = Rtsne(t(as.matrix(data_scrabble)), 
                        initial_dims = initial_dims_value, 
                        perplexity = perplexity_value)
  
  pl[[7]] = plot_pca_singlecell(scrabble_tsne$Y,data_simulation$group)
  
  main = grid.arrange(grobs = pl,ncol = 7)
  
  return(main)
  
}


# plot MA function
plot_ma <- function(dropout_index, seed_value){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  data_simulation = readRDS(file = paste0('simulation_data/simulation_data_dropout_index_',
                                          dropout_index,
                                          '_seed_',
                                          seed_value,
                                          '.rds')
  )
  
  # true data
  data_true = data_simulation$data_true
  
  # raw data (dropout data)
  data_dropout = data_simulation$data_dropout
  
  
  # load the magic results 
  data = read.csv(paste0("imputation_magic_data/magic_",
                         dropout_index,"_",
                         seed_value,".csv"),
                  header = FALSE,
                  sep = "\t")
  data$V1 = NULL
  data_magic = as.matrix(data)
  
  # load imputed data from scImpute
  data_scimpute = read.table( file = paste0("imputation_scimpute_data/scimpute_",
                                            dropout_index, "_",
                                            seed_value,
                                            ".csv") ,
                              header = FALSE, sep=","
  )
  
  
  # load imputed data from Drimpute
  data_drimpute = read.table( file = paste0("imputation_drimpute_data/drimpute_",
                                            dropout_index, "_",
                                            seed_value,
                                            ".csv") ,
                              header = FALSE, sep=","
  )
  
  
  # load imputed data from VIPER
  data_viper = read.table( file = paste0("imputation_viper_data/viper_",
                                         dropout_index, "_",
                                         seed_value,
                                         ".csv") ,
                           header = FALSE, sep=","
  )
  
  # load imputed data from scrabble
  data_scrabble = read.table( file = paste0("imputation_scrabble_data/scrabble_",
                                            dropout_index, "_",
                                            seed_value,
                                            "_", k,
                                            ".csv") ,
                              header = FALSE, sep=","
  )
  
  pl = list()
  # dropout
  pl[[1]] = prepare_ma_data_plot(data_true, data_dropout,
                                 paste0("Dropout rate:", 
                                        round(data_simulation$percentage_zeros),"%")
  ) 
  
  # drimpute
  pl[[2]] = prepare_ma_data_plot(data_true, as.matrix(data_drimpute), "DrImpute") 
  
  # scimpute
  pl[[3]] = prepare_ma_data_plot(data_true, as.matrix(data_scimpute),"scImpute") 
  
  # magic
  pl[[4]] = prepare_ma_data_plot(data_true, as.matrix(data_magic), "MAGIC") 
  
  # viper
  pl[[5]] = prepare_ma_data_plot(data_true, as.matrix(data_viper), "VIPER") 
  
  # scrabble
  pl[[6]] = prepare_ma_data_plot(data_true, as.matrix(data_scrabble), "SCRABBLE") 
  
  main = grid.arrange(grobs = pl,ncol = 6)
  
  return(main)
  
}

# plot MA data function
prepare_ma_data_plot <- function(data_base,data,name_plot){
  
  genes = paste0("gene",c(1:dim(data_base)[1]))
  
  data_ma = data.frame(name = genes)
  
  data_ma$baseMean = rowMeans((log2(data_base + 1) + log2(data + 1))/2)
  
  data_ma$log2FoldChange = rowMeans(log2((data+1)/(data_base+1)))
  
  data_ma$padj = t(calculate_pvalue_two_matrices(data_base,data))
  
  p = ggmaplot1(data_ma,
                 size = 1.5,
                 top = 0,
                 main = name_plot,
                 font.main = c("bold",18),
                 font.legend = c("bold", 15),
                 legend = c(0.8, 0.2),
                 genenames = as.vector(data_ma$name),
                 ylim = c(-8,5)) +
    theme(axis.title = element_text(family = "Helvetica", size = (18)),
          axis.text.x = element_text(family = "Helvetica", size = (15)),
          axis.text.y = element_text(family = "Helvetica", size = (15)))
  
  
  return(p)
  
}

ggmaplot1 <- function (data, fdr = 0.05, fc = 1.5, genenames = NULL,
                       detection_call = NULL, size = NULL,
                       font.label = c(12, "plain", "black"), label.rectangle = FALSE,
                       palette = c("#B31B21", "#1465AC", "darkgray"),
                       top = 15, select.top.method = c("padj", "fc"),
                       main = NULL, xlab = "Log2 mean expression",  ylab = "Log2 fold change",
                       ggtheme = theme_classic(),...)
{
  
  if(!base::inherits(data, c("matrix", "data.frame", "DataFrame", "DE_Results", "DESeqResults")))
    
    stop("data must be an object of class matrix, data.frame, DataFrame, DE_Results or DESeqResults")
  
  if(!is.null(detection_call)){
    if(nrow(data)!=length(detection_call))
      
      stop("detection_call must be a numeric vector of length = nrow(data)")
    
  }
  else if("detection_call" %in% colnames(data)){
    
    detection_call = as.vector(data$detection_call)
    
  }
  else detection_call = rep(1, nrow(data))
  
  # Legend position
  if(is.null(list(...)$legend)) legend = c(0.12, 0.9)
  
  # Check data format
  ss = base::setdiff(c("baseMean", "log2FoldChange", "padj"), colnames(data))
  if(length(ss)>0) stop("The colnames of data must contain: ",
                        paste(ss, collapse = ", "))
  
  if(is.null(genenames)) genenames <- rownames(data)
  else if(length(genenames)!=nrow(data))
    stop("genenames should be of length nrow(data).")
  
  sig <- rep(3, nrow(data))
  
  sig[which(data$padj <= fdr & data$log2FoldChange < 0 & abs(data$log2FoldChange) >= log2(fc) & detection_call ==1)] = 2
  
  sig[which(data$padj <= fdr & data$log2FoldChange > 0 & abs(data$log2FoldChange) >= log2(fc) & detection_call ==1)] = 1
  
  data <- data.frame(name = genenames, mean = data$baseMean, lfc = data$log2FoldChange,
                     padj = data$padj, sig = sig)
  
  # Change level labels
  . <- NULL
  data$sig <- as.factor(data$sig)
  .lev <- .levels(data$sig) %>% as.numeric()
  palette <- palette[.lev]
  new.levels <- c(
    paste0("Up: ", sum(sig == 1)),
    paste0("Down: ", sum(sig == 2)),
    "NS"
  ) %>% .[.lev]
  
  data$sig <- factor(data$sig, labels = new.levels)
  
  
  # Ordering for selecting top gene
  select.top.method <- match.arg(select.top.method)
  if(select.top.method == "padj") data <- data[order(data$padj), ]
  else if(select.top.method == "fc") data <- data[order(abs(data$lfc), decreasing = TRUE), ]
  # select data for top genes
  labs_data <- stats::na.omit(data)
  labs_data <- subset(labs_data, padj <= fdr & name!="" & abs(lfc) >= log2(fc))
  labs_data <- utils::head(labs_data, top)
  
  font.label <- .parse_font(font.label)
  font.label$size <- ifelse(is.null(font.label$size), 12, font.label$size)
  font.label$color <- ifelse(is.null(font.label$color), "black", font.label$color)
  font.label$face <- ifelse(is.null(font.label$face), "plain", font.label$face)
  
  # Plot
  set.seed(42)
  mean <- lfc <- sig <- name <- padj <-  NULL
  p <- ggplot(data, aes(x = mean, y = lfc)) +
    geom_point(aes(color = sig), size = size)
  
  if(label.rectangle){
    p <- p + ggrepel::geom_label_repel(data = labs_data, mapping = aes(label = name),
                                       box.padding = unit(0.35, "lines"),
                                       point.padding = unit(0.3, "lines"),
                                       force = 1, fontface = font.label$face,
                                       size = font.label$size/3, color = font.label$color)
  }
  else{
    p <- p + ggrepel::geom_text_repel(data = labs_data, mapping = aes(label = name),
                                      box.padding = unit(0.35, "lines"),
                                      point.padding = unit(0.3, "lines"),
                                      force = 1, fontface = font.label$face,
                                      size = font.label$size/3, color = font.label$color)
  }
  
  p <- p + scale_x_continuous(breaks=seq(0, max(data$mean), 2))+
    labs(x = xlab, y = ylab, title = main, color = "")+ # to remove legend title use color = ""
    geom_hline(yintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2),
               color = c("black", "black", "black"))
  
  p <- ggpar(p, palette = palette, ggtheme = ggtheme, ...)
  p
}

.levels <- function(x){
  if(!is.factor(x)) x <- as.factor(x)
  levels(x)
}

.parse_font <- function(font){
  if(is.null(font)) res <- NULL
  else if(inherits(font, "list")) res <- font
  else{
    # matching size and face
    size <- grep("^[0-9]+$", font, perl = TRUE)
    face <- grep("plain|bold|italic|bold.italic", font, perl = TRUE)
    if(length(size) == 0) size <- NULL else size <- as.numeric(font[size])
    if(length(face) == 0) face <- NULL else face <- font[face]
    color <- setdiff(font, c(size, face))
    if(length(color) == 0) color <- NULL
    res <- list(size=size, face = face, color = color)
  }
  res
}

get_cor_data <- function(dropout_index, seed_value){
  
  options( warn = -1 )
  
  # load the simulationd data
  
  data_simulation = readRDS(file = paste0('simulation_data/simulation_data_drop_index_',
                                           dropout_index_index,
                                           '_seed_',
                                           seed_value,
                                           '.rds')
  )
  
  index = rowMeans(data_simulation$data_dropout) > 0
  
  data_true = data_simulation$data_true
  
  data_true = data_true[index,]
  
  # cell-cell correlation
  
  data_true_cell = cor(as.matrix((data_true)))
  
  data_true_cell[is.na(data_true_cell)] = 0
  
  # gene-gene correlation
  
  data_true_gene = cor(t((data_true)), method = "pearson")
  
  data_true_gene[is.na(data_true_gene)] = 0
  
  data_dropout = data_simulation$data_dropout
  
  data_dropout = data_dropout[index,]
  
  # cell-cell correlation
  
  data_dropout_cell = cor((data_dropout), method = "pearson")
  
  data_dropout_cell[is.na(data_dropout_cell)] = 0
  
  # gene-gene correlation
  
  data_dropout_gene = cor(t((data_dropout)), method = "pearson")
  
  data_dropout_gene[is.na(data_dropout_gene)] = 0
  
  # load the magic results 
  
  data = read.csv(paste0("imputation_magic_data/magic_",dropout_index,"_",seed_value,".csv"),
                   header = FALSE,
                   sep = "\t")
  
  data$V1 = NULL
  
  data_magic = as.matrix(data)
  
  data_magic[data_magic < 0] = 0
  
  data_magic[is.nan(data_magic)] = 0
  
  data_magic = data_magic[index,]
  
  # cell-cell correlation
  
  data_magic_cell = cor((data_magic), method = "pearson")
  
  data_magic_cell[is.na(data_magic_cell)] = 0
  
  # gene-gene correlation
  
  data_magic_gene = cor(t((data_magic)), method = "pearson")
  
  data_magic_gene[is.na(data_magic_gene)] = 0
  
  # load imputed data from scImpute
  
  data_scimpute = read.table( file = paste0("imputation_scimpute_data/scimpute_",
                                             dropout_index, "_",
                                             seed_value,
                                             ".csv") ,
                               header = FALSE, sep=","
  )
  
  data_scimpute = data_scimpute[index,]
  
  # cell-cell correlation
  
  data_scimpute_cell = cor((data_scimpute), method = "pearson")
  
  data_scimpute_cell[is.na(data_scimpute_cell)] = 0
  
  # gene-gene correlation
  
  data_scimpute_gene = cor(t((data_scimpute)), method = "pearson")
  
  data_scimpute_gene[is.na(data_scimpute_gene)] = 0
  
  # load imputed data from Drimpute
  
  data_drimpute = read.table( file = paste0("imputation_drimpute_data/drimpute_",
                                             dropout_index, "_",
                                             seed_value,
                                             ".csv") ,
                               header = FALSE, sep=","
  )
  
  data_drimpute = data_drimpute[index,]
  
  # cell-cell correlation
  
  data_drimpute_cell = cor((data_drimpute), method = "pearson")
  
  data_drimpute_cell[is.na(data_drimpute_cell)] = 0
  
  # gene-gene correlation
  
  data_drimpute_gene = cor(t((data_drimpute)), method = "pearson")
  
  data_drimpute_gene[is.na(data_drimpute_gene)] = 0
  

  # load imputed data from Saver
  data_viper = read.table( file = paste0("imputation_viper_data/viper_",
                                          dropout_index, "_",
                                          seed_value,
                                          ".csv") ,
                            header = FALSE, sep=","
  )
  
  data_viper = data_viper[index,]
  
  # cell-cell correlation
  
  data_viper_cell = cor((data_viper), method = "pearson")
  
  data_viper_cell[is.na(data_viper_cell)] = 0
  
  # gene-gene correlation
  
  data_viper_gene = cor(t((data_viper)), method = "pearson")
  
  data_viper_gene[is.na(data_viper_gene)] = 0
  
  # load imputed data from scrabble
  
  data_scrabble = read.table( file = paste0("imputation_scrabble_data/scrabble_",
                                             dropout_index, "_",
                                             seed_value,
                                             ".csv") ,
                               header = FALSE, sep=","
  )
  
  data_scrabble = data_scrabble[index,]
  
  
  # cell-cell correlatio
  
  data_scrabble_cell = cor(as.matrix(data_scrabble), method = "pearson")
  
  data_scrabble_cell[is.na(data_scrabble_cell)] = 0
  
  # gene-gene correlation
  
  data_scrabble_gene = cor(t((data_scrabble)), method = "pearson")
  
  data_scrabble_gene[is.na(data_scrabble_gene)] = 0
  
  data_cell = list()
  
  data_cell[[1]] = data_true_cell
  
  data_cell[[2]] = data_dropout_cell
  
  data_cell[[3]] = data_drimpute_cell
  
  data_cell[[4]] = data_scimpute_cell
  
  data_cell[[5]] = data_magic_cell
  
  data_cell[[6]] = data_viper_cell
  
  data_cell[[7]] = data_scrabble_cell
  
  data_gene = list()
  
  data_gene[[1]] = data_true_gene
  
  data_gene[[2]] = data_dropout_gene
  
  data_gene[[3]] = data_drimpute_gene
  
  data_gene[[4]] = data_scimpute_gene
  
  data_gene[[5]] = data_magic_gene
  
  data_gene[[6]] = data_viper_gene
  
  data_gene[[7]] = data_scrabble_gene
  
  data_cor = list()
  
  data_cor[[1]] = data_cell
  
  data_cor[[2]] = data_gene
  
  return(data_cor)
  
}


plot_cell_distribution <- function(dropout_index, seed_value){
  
  methods = c("True Data", "Dropout Data", "DrImpute", "scImpute", "MAGIC", "VIPER", "SCRABBLE")
  # load the simulationd data
  data_simulation = readRDS(file = paste0('simulation_data/simulation_data_drop_index_',
                                           dropout_index,
                                           '_seed_',
                                           seed_value,
                                           '.rds')
  )
  
  group_info = unique(data_simulation$group)
  
  data_cor = get_cor_data(dropout_index, seed_value)
  
  data_cell = data_cor[[1]]
  
  data_gene = data_cor[[2]]
  
  k = 1
  
  p = list()
  
  for( j in group_info){
    
    index1 = data_simulation$group == j
    
    ratio = c()
    
    for(i in c(1:length(methods))){
      
      tmp = data_cell[[i]]
      
      low_index = lower.tri(tmp)*1
      
      tmp_index = matrix(0, nrow = dim(tmp)[1], ncol = dim(tmp)[2])
      
      tmp_index[index1,index1] = 1
      
      low_index1 = low_index*tmp_index > 0 
      
      tmp_index = matrix(0, nrow = dim(tmp)[1], ncol = dim(tmp)[2])
      
      tmp_index[index1,!index1] = 1
      
      tmp_index[!index1,index1] = 1
      
      low_index2 = low_index*tmp_index > 0
      
      data_plot1 = data.frame(x = data_cor_vector(tmp[low_index1]))
      
      data_plot1$y = j
      
      data_plot2 = data.frame(x = data_cor_vector(tmp[low_index2]))
      
      data_plot2$y = "Group0"
      
      data_plot = rbind(data_plot1,data_plot2)
      
      a = ks.test(data_plot1$x, data_plot2$x)
      
      ratio = cbind(ratio, a[[1]])
      
      p[[k]] = ggplot(data_plot, aes(x=x, fill=y, color=y)) +
        geom_density(alpha = .3) +
        ggtitle(paste0("Cell: ", methods[i])) +
        theme(legend.position = c(0.2, 0.9)) 
      
      k = k + 1
      
    }
    
    data_bar = data.frame(y = t(ratio), x = methods)
    
    colnames(data_bar) = c("y","x")
    
    p[[k]] = ggplot(data_bar, aes(x=x, y=y)) +
      geom_bar(stat="identity", fill="steelblue") +
      scale_x_discrete(limits=data_bar$x) +
      xlab("Method") +
      ylab("KS Statistics") +
      ggtitle(paste0(j, "_", round(data_simulation$percentage_zeros))) +
      theme_cowplot() +
      theme(plot.title = element_text(size = 18, hjust = 0.4),
            axis.text = element_text(size = 12),
            legend.title=element_blank())
    
    k = k + 1
    
  }
  
  main = grid.arrange(grobs = p,ncol = (length(methods) + 1))
  
  return(main)
  
}


plot_gene_distribution <- function(dropout_index, seed_value){
  
  methods = c("True Data", "Dropout Data", "DrImpute", "scImpute", "MAGIC", "VIPER", "SCRABBLE")
  # load the simulationd data
  data_simulation = readRDS(file = paste0('simulation_data/simulation_data_drop_index_',
                                           dropout_index,
                                           '_seed_',
                                           seed_value,
                                           '.rds')
  )
  
  group_info = unique(data_simulation$group)
  
  data_cor = get_cor_data(dropout_index, seed_value)
  
  data_cell = data_cor[[1]]
  
  data_gene = data_cor[[2]]
  
  tmp = as.numeric(as.numeric(gsub("Group", "", data_simulation$group)))
  
  de = get_marker_genes(data_simulation$data_true, tmp)
  
  group_unique = unique(tmp)
  
  N_groups = length(group_unique)
  
  de_gene = list()
  
  for(i in c(1:N_groups)){
  
    de_gene[[i]] = paste0("Gene",which((de$auroc > 0.85) & (de$clusts == group_unique[i]) & (de$pvalue < 0.01)))
  
  }
  
  names_gene = rownames(data_gene[[1]])
  
  k = 1
  
  p = list()
  
  for( j in c(1:N_groups)){
    
    index1 = names_gene %in% de_gene[[j]]
    
    ratio = c()
    
    for(i in c(1:length(methods))){
      
      tmp = data_gene[[i]]
      
      low_index = lower.tri(tmp)*1
      
      tmp_index = matrix(0, nrow = dim(tmp)[1], ncol = dim(tmp)[2])
      
      tmp_index[index1,index1] = 1
      
      low_index1 = low_index*tmp_index > 0 
      
      tmp_index = matrix(0, nrow = dim(tmp)[1], ncol = dim(tmp)[2])
      
      tmp_index[index1,!index1] = 1
      
      tmp_index[!index1,index1] = 1
      
      low_index2 = low_index*tmp_index > 0
      
      data_plot1 = data.frame(x = data_cor_vector(tmp[low_index1]))
      
      data_plot1$y = paste0("Marker_",j)
      
      data_plot2 = data.frame(x = data_cor_vector(tmp[low_index2]))
      
      data_plot2$y = "NonMarker"
      
      data_plot = rbind(data_plot1,data_plot2)
      
      a = ks.test(data_plot1$x, data_plot2$x)
      
      ratio = cbind(ratio, a[[1]])
      
      p[[k]] = ggplot(data_plot, aes(x=x, fill=y, color=y)) +
        geom_density(alpha = .3) +
        ggtitle(paste0("Gene: ", methods[i])) +
        theme(legend.position = c(0.2, 0.9)) 
      
      k = k + 1
      
    }
    
    data_bar = data.frame(y = t(ratio), x = methods)
    
    colnames(data_bar) = c("y","x")
    
    p[[k]] = ggplot(data_bar, aes(x=x, y=y)) +
      geom_bar(stat="identity", fill="steelblue") +
      scale_x_discrete(limits=data_bar$x) +
      xlab("Method") +
      ylab("KS Statistics") +
      ggtitle(paste0(j, "_", round(data_simulation$percentage_zeros))) +
      theme_cowplot() +
      theme(plot.title = element_text(size = 18, hjust = 0.4),
            axis.text = element_text(size = 12),
            legend.title=element_blank())
    
    k = k + 1
    
  }
  
  main = grid.arrange(grobs = p,ncol = 9)
  
  return(main)
  
}


data_cor_vector <- function(data){
  
  return(data[lower.tri(data)])
  
}

cal_cell_distribution <- function(dropout_index, seed_value){
  
  methods = c("True Data", "Dropout Data", "DrImpute", "scImpute", "MAGIC", "VIPER", "SCRABBLE")
  
  # load the simulationd data
  
  data_simulation = readRDS(file = paste0('simulation_data/simulation_data_drop_index_',
                                           dropout_index,
                                           '_seed_',
                                           seed_value,
                                           '.rds')
  )
  
  group_info = unique(data_simulation$group)
  
  data_cor = get_cor_data(dropout_index, seed_value)
  
  data_cell = data_cor[[1]]
  
  data_gene = data_cor[[2]]
  
  k = 1
  
  p = list()
  
  ratio1 = c()
  
  for( j in group_info){
    
    index1 = data_simulation$group == j
    
    ratio = c()
    
    for(i in c(1:length(methods))){
      
      tmp = data_cell[[i]]
      
      low_index = lower.tri(tmp)*1
      
      tmp_index = matrix(0, nrow = dim(tmp)[1], ncol = dim(tmp)[2])
      
      tmp_index[index1,index1] = 1
      
      low_index1 = low_index*tmp_index > 0 
      
      tmp_index = matrix(0, nrow = dim(tmp)[1], ncol = dim(tmp)[2])
      
      tmp_index[index1,!index1] = 1
      
      tmp_index[!index1,index1] = 1
      
      low_index2 = low_index*tmp_index > 0
      
      data_plot1 = data.frame(x = data_cor_vector(tmp[low_index1]))
      
      data_plot1$y = j
      
      data_plot2 = data.frame(x = data_cor_vector(tmp[low_index2]))
      
      data_plot2$y = "Group0"
      
      data_plot = rbind(data_plot1,data_plot2)
      
      a = ks.test(data_plot1$x, data_plot2$x)
      
      ratio = cbind(ratio, a[[1]])
      
    }
    
    ratio1 = rbind(ratio1,abs(ratio[2:length(methods)] - ratio[1]))
  }
  
  saveRDS(colMeans(ratio1),file = paste0("data_cell_distribution/data_",dropout_index,"_",seed_value,".rds"))
  
}

cal_gene_distribution <- function(dropout_index, seed_value){
  
  methods = c("True Data", "Dropout Data", "DrImpute", "scImpute", "MAGIC", "VIPER", "SCRABBLE")
  
  # load the simulationd data
  data_simulation = readRDS(file = paste0('simulation_data/simulation_data_drop_index_',
                                           dropout_index,
                                           '_seed_',
                                           seed_value,
                                           '.rds')
  )
  
  group_info = unique(data_simulation$group)
  
  data_cor = get_cor_data(drop_index, seed_value)
  
  data_cell = data_cor[[1]]
  
  data_gene = data_cor[[2]]
  
  tmp = as.numeric(as.numeric(gsub("Group", "", data_simulation$group)))
  
  de = get_marker_genes(data_simulation$data_true, tmp)
  
  group_unique = unique(tmp)
  
  N_groups = length(group_unique)
  
  de_gene = list()
  
  for(i in c(1:N_groups)){
    
    de_gene[[i]] = paste0("Gene",which((de$auroc > 0.85) & (de$clusts == group_unique[i]) & (de$pvalue < 0.01)))
    
  }
  
  names_gene = rownames(data_gene[[1]])
  
  k = 1
  
  p = list()
  
  ratio1 = c()
  
  for( j in c(1:N_groups)){
    
    if(length(de_gene[[j]]) > 2){
      
      index1 = names_gene %in% de_gene[[j]]
      
      ratio = c()
      
      for(i in c(1:length(methods))){
        
        tmp = data_gene[[i]]
        
        low_index = lower.tri(tmp)*1
        
        tmp_index = matrix(0, nrow = dim(tmp)[1], ncol = dim(tmp)[2])
        
        tmp_index[index1,index1] = 1
        
        low_index1 = low_index*tmp_index > 0 
        
        tmp_index = matrix(0, nrow = dim(tmp)[1], ncol = dim(tmp)[2])
        
        tmp_index[index1,!index1] = 1
        
        tmp_index[!index1,index1] = 1
        
        low_index2 = low_index*tmp_index > 0
        
        data_plot1 = data.frame(x = data_cor_vector(tmp[low_index1]))
        
        data_plot1$y = paste0("Marker_",j)
        
        data_plot2 = data.frame(x = data_cor_vector(tmp[low_index2]))
        
        data_plot2$y = "NonMarker"
        
        data_plot = rbind(data_plot1,data_plot2)
        
        a = ks.test(data_plot1$x, data_plot2$x)
        
        ratio = cbind(ratio, a[[1]])
        
      }
      
      ratio1 = rbind(ratio1,abs(ratio[2:length(methods)] - ratio[1]))
      
    }
    
  }
  
  # return(ratio1)  
  saveRDS(colMeans(ratio1),file = paste0("data_gene_distribution/data_",drop_index,"_",seed_value,".rds"))
  
}







