# R function to convert data size when they are in different format and extarcted by ls -ltrh command 
# input is a column which has data size with example 108K, 92G, 67M etc.
size_converter <- function(string){
  converted_size <- c()
  for (size in string){
    # print(size)
    temp_unit <- gsub("[0-9 . ]","",size)
    temp_value <- as.numeric(gsub("[A-Za-z ]","",size))
    # print(temp)
    bytes_size <- case_when(
      temp_unit == "" ~ 0,
      temp_unit == "B" ~ 1024,
      temp_unit == "K" ~ 1024^2,
      temp_unit == "M" ~ 1024^2,
      temp_unit == "G" ~ 1024^3,
      temp_unit == "T" ~ 1024^4,
      temp_unit == "P" ~ 1024^5,
      .default = 0 
    )
    calc_size <- temp_value * bytes_size
    converted_size <- c(converted_size, calc_size)
    # print(converted_size)
  }
  return(converted_size)
}

# execute
# df["converted_size"] <- size_converter(df$V1.x)
