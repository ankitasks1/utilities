Upsetout <- function(input_list){
  if(missing(input_list)){
    stop("No input list is provided. Please provide an appropriate input list")
  }
  if(is.null(input_list)){
    stop("Input list is empty. Please add data your list")
  }else{
    message("Input list is correctly provided. performing interactions")
    # Combine all elements in all the the vectors
    all_elements <- unique(as.character(unlist(input_list)))
    # Create the binary list
    binary_list <- lapply(input_list, function(x) as.numeric(all_elements%in%x))
    #Bind the binary matrix
    binary_df <- data.frame(do.call(cbind.data.frame, binary_list))
    rownames(binary_df) <- all_elements
    # make a dataframe which adress intersection of all elemnts and respective samples
    for (names_i in colnames(binary_df)){
      binary_df_names_i <- binary_df[,names_i]
      #Insert name of column into the table elements
      binary_df[names_i] <- ifelse(binary_df_names_i >= 1, names_i, 0)
      
    }
    binary_df$combos <- apply(binary_df[ , colnames(binary_df)] , 1 , paste , collapse = "_" )
    
    # Prepare the max combination required
    combos_for_venn <- expand.grid(rep(list(0:1),length(names(input_list))))
    colnames(combos_for_venn) <- names(input_list)
    for (names_i in names(combos_for_venn)){
      combos_for_venn_names_i <- combos_for_venn[,names_i]
      #print(head(temp50_names_i))
      combos_for_venn[names_i] <- ifelse(combos_for_venn_names_i >= 1, names_i, 0)
    }
    combos_for_venn$combos <- apply(combos_for_venn[ , colnames(combos_for_venn)] , 1 , paste , collapse = "_" )
    # Extract combinations as output list
    finalout_upsetlist <- list()
    for (combos in combos_for_venn$combos){
      finalout_upsetlist[[combos]] <- binary_df[which(binary_df$combos == combos),]
    }
    # print(finalout_upsetlist)
    message("Combinations are ready. Explore the output list assigned by you")
    return(finalout_upsetlist)
  }
}