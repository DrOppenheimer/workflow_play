sigtest <- function(data_file="filtered_counts.txt", 
                    metadata_file="filtered_counts.metadata.txt", #filtered_counts.metadata.txt", 
                    metadata_column="env_package.data.env_package", #env_package.data.env_package", #  env_package.data.body_site
                    stat_test="ANOVA-one-way", # c("Kruskal-Wallis","t-test-paired","Wilcoxon-paired","ANOVA-one-way","t-test-unpaired" = ,"Mann-Whitney-unpaired-Wilcoxon") 
                    p_adjust_method = "BH" # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")
                    ){
  
  # load some data and metadata to play with
  my_data <- import_data(data_file)
  my_metadata <- import_metadata(metadata_file)
  my_metadata <<- my_metadata
  
  # Here, make sure that the data are sorted COLUMNWISE by id
  my_data <-  my_data[,order(colnames(my_data))]
  
  # make sure that the metadata matrix is sorted (ROWWISE) by id
  my_metadata <-  my_metadata[order(rownames(my_metadata)),]
  
  # check to make sure that the two listings of ids are identical
  if( identical(colnames(my_data), rownames(my_metadata)) ){
    print("IDs are EQUAL")
  }else{
    stop("IDs are NOT EQUAL, check your data and metadata files")
  }

  # use a switch to allow the user to select the test --- yes, your grandmama really could write better code ^_^
  switch(stat_test, 
         "Kruskal-Wallis" = perform_kw(data_file,metadata_file,metadata_column,stat_test,p_adjust_method, my_data, my_metadata),
         "t-test-paired" = perform_ptt(data_file,metadata_file,metadata_column,stat_test,p_adjust_method, my_data, my_metadata),
         "Wilcoxon-paired" = perform_pw(data_file,metadata_file,metadata_column,stat_test,p_adjust_method, my_data, my_metadata),
         "t-test-unpaired" = perform_uptt(data_file,metadata_file,metadata_column,stat_test,p_adjust_method, my_data, my_metadata),
         "Mann-Whitney-unpaired-Wilcoxon" = perform_upw(data_file,metadata_file,metadata_column,stat_test,p_adjust_method, my_data, my_metadata), 
         "ANOVA-one-way" = perform_anova(data_file,metadata_file,metadata_column,stat_test,p_adjust_method,my_data, my_metadata, debug="TRUE")
         )
}

import_metadata <- function(group_table){ #, group_column, sample_names){
  metadata_matrix <- as.matrix( # Load the metadata table (same if you use one or all columns)
    read.table(
      file=group_table,row.names=1,header=TRUE,sep="\t",
      colClasses = "character", check.names=FALSE,
      comment.char = "",quote="",fill=TRUE,blank.lines.skip=FALSE
    )
  )
}

import_data <- function(file_name)
{
  data.matrix(read.table(file_name, row.names=1, header=TRUE, sep="\t", comment.char="", quote="", check.names=FALSE))
}       

export_data <- function(data_object, file_name){
  write.table(data_object, file=file_name, sep="\t", col.names = NA, row.names = TRUE, quote = FALSE, eol="\n")
} 
    
perform_anova <- function(
    data_file, #="filtered_counts.txt", 
    metadata_file, #="filtered_counts.metadata.txt", 
    metadata_column, #="env_package.data.body_site", 
    stat_test, #="ANOVA-one-way",
    p_adjust_method, # = "BH",
    my_data,
    my_metadata,
    debug = FALSE
  ){

  # prep someplace to store the stats for all of the rows in my_data
  my_stats <- matrix(nrow = nrow(my_data), ncol=7)
  rownames(my_stats) <- rownames(my_data)
  
  # parts that need to be part of the switch
  colnames(my_stats) <- c("median", "mean", "sd", "F_stat", "p", "bonferroni_p", paste(p_adjust_method,"_p",sep=""))

  # iterate through each row 
  for (i in 1:nrow(my_data)){
  
    # first calculate some simple state
    my_stats[i,"median"] <- median(my_data[i,])
    my_stats[i,"mean"] <- mean(my_data[i,])
    my_stats[i,"sd"] <- sd(my_data[i,])
  
    # prep data for anova
    #stat_input <- matrix(nrow=ncol(my_data), ncol=2) # a
    stat_input <- matrix(nrow=ncol(my_data), ncol=2) # a
    colnames(stat_input) <- c("values","groups")
    stat_input[,"values"] <- my_data[i,]
    #if(identical(debug,"TRUE")){my_data<<-my_data;stat_input<<-stat_input; my_metadata <<- my_metadata ;stop("HOLD ON THERE")}
    stat_input[,"groups"] <- my_metadata[,metadata_column]
    
    stat_input <- as.data.frame(stat_input)
  
    # Perform ANOVA using the ~ operator
    aov_result <- aov(values ~ groups, data = stat_input)
  
    # Summary of ANOVA results
    aov_result_summary <- summary(aov_result)
  
    # Get ANOVA results into my_stats (note FDR and adjusted p have to be calculated later)
    my_stats[i,"F_stat"] <- aov_result_summary[[1]]$`F value`[1]
    my_stats[i,"p"]      <- aov_result_summary[[1]]$`Pr(>F)`[1]
  
  }
  
  # Calculate the Bonferroni adjusted p
  my_stats[,"bonferroni_p"] <- p.adjust(p=my_stats[,"p"], method = "bonferroni")

  # Calculate the Benjamini & Hochberg adjusted p
  my_stats[,"BH_p"] <- p.adjust(p=my_stats[,"p"], method = p_adjust_method)

  # combine my_data and my_stats to create a single output object
  my_output_data <- cbind(my_data,my_stats)

  # sort the data by p value  
  my_output_data <- my_output_data[order(my_output_data[, "p"]), ]

  # output the object
  export_data(
    data_object = my_output_data, 
    file_name = paste(tools::file_path_sans_ext(data_file),".",stat_test,".", metadata_column, ".STAT_RESULTS.txt", sep="")
  )

}





perform_kw <- function(
    data_file, #="filtered_counts.txt", 
    metadata_file, #="filtered_counts.metadata.txt", 
    metadata_column, #="env_package.data.body_site", 
    stat_test, #="Kruskal-Wallis",
    p_adjust_method, # = "BH",
    my_data,
    my_metadata
){
  
  # prep someplace to store the stats for all of the rows in my_data
  my_stats <- matrix(nrow = nrow(my_data), ncol=7)
  rownames(my_stats) <- rownames(my_data)
  
  # parts that need to be part of the switch
  colnames(my_stats) <- c("median", "mean", "sd", "chi-squared_stat", "p", "bonferroni_p", paste(p_adjust_method,"_p",sep=""))
  
  # iterate through each row 
  for (i in 1:nrow(my_data)){
    
    # first calculate some simple state
    my_stats[i,"median"] <- median(my_data[i,])
    my_stats[i,"mean"] <- mean(my_data[i,])
    my_stats[i,"sd"] <- sd(my_data[i,])
    
    # prep data for KW
    stat_input <- matrix(nrow=ncol(my_data), ncol=2) # a 
    colnames(stat_input) <- c("values","groups")
    stat_input[,"values"] <- my_data[i,]
    stat_input[,"groups"] <- my_metadata[,metadata_column]
    stat_input <- as.data.frame(stat_input)
    
    # Perform KW using the ~ operator
    kw_result <- kruskal.test(values ~ groups, data = stat_input)
    
    # Get KW results into my_stats (note FDR and adjusted p have to be calculated later)
    my_stats[i,"chi-squared_stat"] <- as.numeric(kw_result["statistic"])
    my_stats[i,"p"]      <- as.numeric(kw_result["p.value"])
    
  }
  
  # Calculate the Bonferroni adjusted p
  my_stats[,"bonferroni_p"] <- p.adjust(p=my_stats[,"p"], method = "bonferroni")
  
  # Calculate the Benjamini & Hochberg adjusted p
  my_stats[,"BH_p"] <- p.adjust(p=my_stats[,"p"], method = p_adjust_method)
  
  # combine my_data and my_stats to create a single output object
  my_output_data <- cbind(my_data,my_stats)
  
  # sort the data by p value  
  my_output_data <- my_output_data[order(my_output_data[, "p"]), ]
  
  # output the object
  export_data(
    data_object = my_output_data, 
    file_name = paste(tools::file_path_sans_ext(data_file),".",stat_test,".", metadata_column, ".STAT_RESULTS.txt", sep="")
  )
  
}


perform_ptt <- function(
    data_file, #="filtered_counts.txt", 
    metadata_file, #="filtered_counts.metadata.txt", 
    metadata_column, #="env_package.data.env_package", 
    stat_test, #="t-test-paired",
    p_adjust_method, # = "BH",
    my_data,
    my_metadata
){
  
  # prep someplace to store the stats for all of the rows in my_data
  my_stats <- matrix(nrow = nrow(my_data), ncol=7)
  rownames(my_stats) <- rownames(my_data)
  
  # parts that need to be part of the switch
  colnames(my_stats) <- c("median", "mean", "sd", "t_stat", "p", "bonferroni_p", paste(p_adjust_method,"_p",sep=""))
  
  # iterate through each row 
  for (i in 1:nrow(my_data)){
    
    # first calculate some simple state
    my_stats[i,"median"] <- median(my_data[i,])
    my_stats[i,"mean"] <- mean(my_data[i,])
    my_stats[i,"sd"] <- sd(my_data[i,])
    
    # prep data for paired t-test
    stat_input <- matrix(nrow=ncol(my_data), ncol=2) # a 
    colnames(stat_input) <- c("values","groups")
    stat_input[,"values"] <- my_data[i,]
    stat_input[,"groups"] <- my_metadata[,metadata_column]
    
    #check to make sure that there are only two groups
    quick_test <- length(unique(stat_input[,"groups"]))
    if( quick_test != 2 ){
      stop(paste("You have", quick_test,"groups, not 2. Check your metadata"))
    }
    ##stat_input <- as.data.frame(stat_input)
    # Separate the first column based on the values in the second column
    split_input <- split(stat_input[, "values"], stat_input[, "groups"])
    # check to make sure that the two groups are the same size
    if( length(split_input[[1]]) != length(split_input[[2]]) ){
      stop("Groups are not the same size. Check your metadata")
    }
    
    # Perform paired t-test 
    ptt_result <- t.test(as.numeric(split_input[[1]]), as.numeric(split_input[[2]]), paired = TRUE)
    
    # Get ptt results into my_stats (note FDR and adjusted p have to be calculated later)
    my_stats[i,"t_stat"] <- as.numeric(ptt_result["statistic"])
    my_stats[i,"p"]      <- as.numeric(ptt_result["p.value"])
    
  }
  
  # Calculate the Bonferroni adjusted p
  my_stats[,"bonferroni_p"] <- p.adjust(p=my_stats[,"p"], method = "bonferroni")
  
  # Calculate the Benjamini & Hochberg adjusted p
  my_stats[,"BH_p"] <- p.adjust(p=my_stats[,"p"], method = p_adjust_method)
  
  # combine my_data and my_stats to create a single output object
  my_output_data <- cbind(my_data,my_stats)
  
  # sort the data by p value  
  my_output_data <- my_output_data[order(my_output_data[, "p"]), ]
  
  # output the object
  export_data(
    data_object = my_output_data, 
    file_name = paste(tools::file_path_sans_ext(data_file),".",stat_test,".", metadata_column, ".STAT_RESULTS.txt", sep="")
  )
  
}



perform_uptt <- function(
    data_file, #="filtered_counts.txt", 
    metadata_file, #="filtered_counts.metadata.txt", 
    metadata_column, #="env_package.data.env_package", 
    stat_test, #="t-test-paired",
    p_adjust_method,
    my_data,
    my_metadata # = "BH"
){
  
  # prep someplace to store the stats for all of the rows in my_data
  my_stats <- matrix(nrow = nrow(my_data), ncol=7)
  rownames(my_stats) <- rownames(my_data)
  
  # parts that need to be part of the switch
  colnames(my_stats) <- c("median", "mean", "sd", "t_stat", "p", "bonferroni_p", paste(p_adjust_method,"_p",sep=""))
  
  # iterate through each row 
  for (i in 1:nrow(my_data)){
    
    # first calculate some simple state
    my_stats[i,"median"] <- median(my_data[i,])
    my_stats[i,"mean"] <- mean(my_data[i,])
    my_stats[i,"sd"] <- sd(my_data[i,])
    
    # prep data for unpaired t-test
    stat_input <- matrix(nrow=ncol(my_data), ncol=2) # a 
    colnames(stat_input) <- c("values","groups")
    stat_input[,"values"] <- my_data[i,]
    stat_input[,"groups"] <- my_metadata[,metadata_column]
    
    #check to make sure that there are only two groups
    quick_test <- length(unique(stat_input[,"groups"]))
    if( quick_test != 2 ){
      stop(paste("You have", quick_test,"groups, not 2. Check your metadata"))
    }
    ##stat_input <- as.data.frame(stat_input)
    # Separate the first column based on the values in the second column
    split_input <- split(stat_input[, "values"], stat_input[, "groups"])
    
    # Perform unpaired t-test 
    uptt_result <- t.test(as.numeric(split_input[[1]]), as.numeric(split_input[[2]]), paired = FALSE)
    
    # Get ptt results into my_stats (note FDR and adjusted p have to be calculated later)
    my_stats[i,"t_stat"] <- as.numeric(uptt_result["statistic"])
    my_stats[i,"p"]      <- as.numeric(uptt_result["p.value"])
    
  }
  
  # Calculate the Bonferroni adjusted p
  my_stats[,"bonferroni_p"] <- p.adjust(p=my_stats[,"p"], method = "bonferroni")
  
  # Calculate the Benjamini & Hochberg adjusted p
  my_stats[,"BH_p"] <- p.adjust(p=my_stats[,"p"], method = p_adjust_method)
  
  # combine my_data and my_stats to create a single output object
  my_output_data <- cbind(my_data,my_stats)
  
  # sort the data by p value  
  my_output_data <- my_output_data[order(my_output_data[, "p"]), ]
  
  # output the object
  export_data(
    data_object = my_output_data, 
    file_name = paste(tools::file_path_sans_ext(data_file),".",stat_test,".", metadata_column, ".STAT_RESULTS.txt", sep="")
  )
  
}


perform_pw <- function(
    data_file, #="filtered_counts.txt", 
    metadata_file, #="filtered_counts.metadata.txt", 
    metadata_column, #="env_package.data.env_package", 
    stat_test, #="t-test-paired",
    p_adjust_method,
    my_data,
    my_metadata# = "BH"
){
  
  # prep someplace to store the stats for all of the rows in my_data
  my_stats <- matrix(nrow = nrow(my_data), ncol=7)
  rownames(my_stats) <- rownames(my_data)
  
  # parts that need to be part of the switch
  colnames(my_stats) <- c("median", "mean", "sd", "v_stat", "p", "bonferroni_p", paste(p_adjust_method,"_p",sep=""))
  
  # iterate through each row 
  for (i in 1:nrow(my_data)){
    
    # first calculate some simple state
    my_stats[i,"median"] <- median(my_data[i,])
    my_stats[i,"mean"] <- mean(my_data[i,])
    my_stats[i,"sd"] <- sd(my_data[i,])
    
    # prep data for paired wilcoxon test
    stat_input <- matrix(nrow=ncol(my_data), ncol=2) # a 
    colnames(stat_input) <- c("values","groups")
    stat_input[,"values"] <- my_data[i,]
    stat_input[,"groups"] <- my_metadata[,metadata_column]
    
    #check to make sure that there are only two groups
    quick_test <- length(unique(stat_input[,"groups"]))
    if( quick_test != 2 ){
      stop(paste("You have", quick_test,"groups, not 2. Check your metadata"))
    }
    ##stat_input <- as.data.frame(stat_input)
    # Separate the first column based on the values in the second column
    split_input <- split(stat_input[, "values"], stat_input[, "groups"])
    # check to make sure that the two groups are the same size
    if( length(split_input[[1]]) != length(split_input[[2]]) ){
      stop("Groups are not the same size. Check your metadata")
    }
    
    # Perform paired pw-test 
    pw_result <- t.test(as.numeric(split_input[[1]]), as.numeric(split_input[[2]]), paired = TRUE)
    
    # Get ptt results into my_stats (note FDR and adjusted p have to be calculated later)
    my_stats[i,"v_stat"] <- as.numeric(pw_result["statistic"])
    my_stats[i,"p"]      <- as.numeric(pw_result["p.value"])
    
  }
  
  # Calculate the Bonferroni adjusted p
  my_stats[,"bonferroni_p"] <- p.adjust(p=my_stats[,"p"], method = "bonferroni")
  
  # Calculate the Benjamini & Hochberg adjusted p
  my_stats[,"BH_p"] <- p.adjust(p=my_stats[,"p"], method = p_adjust_method)
  
  # combine my_data and my_stats to create a single output object
  my_output_data <- cbind(my_data,my_stats)
  
  # sort the data by p value  
  my_output_data <- my_output_data[order(my_output_data[, "p"]), ]
  
  # output the object
  export_data(
    data_object = my_output_data, 
    file_name = paste(tools::file_path_sans_ext(data_file),".",stat_test,".", metadata_column, ".STAT_RESULTS.txt", sep="")
  )
  
}


perform_upw <- function(
    data_file, #="filtered_counts.txt", 
    metadata_file, #="filtered_counts.metadata.txt", 
    metadata_column, #="env_package.data.env_package", 
    stat_test, #="t-test-paired",
    p_adjust_method,
    my_data,
    my_metadata # = "BH"
){
  
  # prep someplace to store the stats for all of the rows in my_data
  my_stats <- matrix(nrow = nrow(my_data), ncol=7)
  rownames(my_stats) <- rownames(my_data)
  
  # parts that need to be part of the switch
  colnames(my_stats) <- c("median", "mean", "sd", "w_stat", "p", "bonferroni_p", paste(p_adjust_method,"_p",sep=""))
  
  # iterate through each row 
  for (i in 1:nrow(my_data)){
    
    # first calculate some simple state
    my_stats[i,"median"] <- median(my_data[i,])
    my_stats[i,"mean"] <- mean(my_data[i,])
    my_stats[i,"sd"] <- sd(my_data[i,])
    
    # prep data for unpaired wilcoxon
    stat_input <- matrix(nrow=ncol(my_data), ncol=2) # a 
    colnames(stat_input) <- c("values","groups")
    stat_input[,"values"] <- my_data[i,]
    stat_input[,"groups"] <- my_metadata[,metadata_column]
    
    #check to make sure that there are only two groups
    quick_test <- length(unique(stat_input[,"groups"]))
    if( quick_test != 2 ){
      stop(paste("You have", quick_test,"groups, not 2. Check your metadata"))
    }
    ##stat_input <- as.data.frame(stat_input)
    # Separate the first column based on the values in the second column
    split_input <- split(stat_input[, "values"], stat_input[, "groups"])
    
    # Perform unpaired t-test 
    upw_result <- t.test(as.numeric(split_input[[1]]), as.numeric(split_input[[2]]), paired = FALSE)
    
    # Get ptt results into my_stats (note FDR and adjusted p have to be calculated later)
    my_stats[i,"w_stat"] <- as.numeric(upw_result["statistic"])
    my_stats[i,"p"]      <- as.numeric(upw_result["p.value"])
    
  }
  
  # Calculate the Bonferroni adjusted p
  my_stats[,"bonferroni_p"] <- p.adjust(p=my_stats[,"p"], method = "bonferroni")
  
  # Calculate the Benjamini & Hochberg adjusted p
  my_stats[,"BH_p"] <- p.adjust(p=my_stats[,"p"], method = p_adjust_method)
  
  # combine my_data and my_stats to create a single output object
  my_output_data <- cbind(my_data,my_stats)
  
  # sort the data by p value  
  my_output_data <- my_output_data[order(my_output_data[, "p"]), ]
  
  # output the object
  export_data(
    data_object = my_output_data, 
    file_name = paste(tools::file_path_sans_ext(data_file),".",stat_test,".", metadata_column, ".STAT_RESULTS.txt", sep="")
  )
  
}

