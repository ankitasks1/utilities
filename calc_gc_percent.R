# fasta <- readLines("~/Downloads/random_fasta.fa")
# read files
args<-commandArgs(TRUE)
fasta <- readLines(args[1])

# ignore >. And collapse all fasta seqeunce in one string
sequencef <- c()
for (f in fasta){
  if  (!grepl(">", f)){
    sequencef <- paste0(f, sequencef)
  }
}
var <- toupper(sequencef)
vars <- unlist(strsplit(var, ""))
len_var <- length(vars)

# calculate GC percent
gc_count <- 0
tot_dublets <- c()
for (i in seq(len_var)){
  dublets <- paste0(vars[i], vars[i+1])
  tot_dublets <- c(tot_dublets, dublets)
  if (dublets == "GC"){
    # print(dublets)
    gc_count <- gc_count + 1
  }
}

gc_per <- (gc_count / len_var) * 100
print(paste0("GC percentage : ",gc_per))


