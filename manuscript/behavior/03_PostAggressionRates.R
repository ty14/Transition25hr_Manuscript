# get post aggression information. 

library(tidyverse)


### Start/End Times for Cage Aggression
vals <-c("Start","End")

# Get the full path to the data_raw subfolder
data_raw_path <- file.path(getwd(), "data_raw/")


# Get a list of all csv files in the data_raw subfolder
files <- list.files(path = data_raw_path, pattern = "_Post.*\\.csv", full.names = TRUE)

# Initialize empty list to store csv data
csv_list <- list()

# Loop over csv files and import them
for (file in files) {
  # Import csv data
  csv_data <- read.csv(file)
  
  # Add csv data to list
  csv_list[[file]] <- csv_data
}

## List of files
sends <- lapply(csv_list, function(x) x[x$Start.End %in% vals,c(1,10)])

names(sends) <- stringr::str_extract(names(sends), "Batch\\d+")


sends

# batch12, add start sends[[4]]
sends[[4]]<- rbind(data.frame(Timestamp="1/16/2020 9:36:12",Start.End="Start"),
                 sends[[4]][,])

# batch14, add start sends[[6]]

sends[[6]] <- rbind(data.frame(Timestamp="2/1/2020 9:36:12",Start.End="Start"),
                 sends[[6]][,])
# batch2, add start sends[[7]]

sends[[7]] <-  rbind(data.frame(Timestamp="10/24/2019 10:41:12",Start.End="Start"),
                     sends[[7]][,])

# batch3, add move end [[8]]
sends[[8]] <- rbind(sends[[8]][c(1:3,6,4:5),])

# batch6, add start sends[[11]]
sends[[11]] <-  rbind(data.frame(Timestamp="11/21/2019 9:36:12",Start.End="Start"),
                     sends[[11]][,])

# batch9, add start sends[[14]]
sends[[14]] <- rbind(sends[[14]][1:2,])




## Function to get differences in starts/ends
get_time_dif<-function(x){
  # Convert timestamp to proper date-time format
  timestamp <- as.POSIXct(x$Timestamp, format="%m/%d/%Y %H:%M:%S")
  # Calculate number of seconds since midnight
  seconds_since_midnight <- difftime(timestamp, trunc(timestamp, "day"), units="secs")
  # Compute differences between consecutive elements
  diffs <- diff(seconds_since_midnight)
  # Extract every second difference starting from the first
  even_diffs <- diffs[seq(1, length(diffs), by=2)]
  return(even_diffs)
}

#Apply to all batches
lapply(sends,get_time_dif)

lapply(lapply(sends,get_time_dif),function(x) as.numeric(sum(x)))

ts <- lapply(lapply(sends,get_time_dif),function(x) as.numeric(sum(x))) 
tsx <- do.call(rbind, ts) %>% as.data.frame()
tsx$V1 <- as.numeric(tsx$V1)
mean(tsx$V1) #36720.21 secs 
sec <- 36720.21/10
min <- sec/60
min/60


# from pre_startend_times
obs.secs <- lapply(lapply(sends,get_time_dif),function(x) as.numeric(sum(x)))


# get post reorganization behavior data
df1 <- read_csv("behavior/Post_WL.csv")

# create list of dataframes for each individual cage
l1 <- split(df1, df1$post_batchcage)

# just keep winner loser columns
l1.wl <- l1 %>% map(~ select(., winner,loser) %>% as.data.frame) 

# convert list to win-loss matrices organized by David's Scores
l1.mat <- l1.wl %>% map(~ compete::get_wl_matrix(.)) %>% map(~ compete::org_matrix(., method="ds"))

# melt matrices.

library(reshape2)


# Define function to melt a matrix and calculate sum of row values and column values
melt_matrix <- function(mat) {
  mat_melted <- melt(mat)   # Melt matrix to long format
  row_sum <- aggregate(value ~ winner, data = mat_melted, sum)   # Sum row values
  col_sum <- aggregate(value ~ loser, data = mat_melted, sum)   # Sum column values
  sum_df <- merge(row_sum, col_sum, by.x = "winner", by.y = "loser", all = TRUE)   # Combine row sum and column sum
  names(sum_df) <- c("postid", "given", "received")  # Rename columns
  return(sum_df)
}

l.aggr <- lapply(l1.mat, melt_matrix)
l.aggr.df <- data.table::rbindlist(Map(cbind, l.aggr, name = names(l.aggr)))
l.aggr.df$batch <- stringr::str_extract(l.aggr.df$name, "\\d+")
l.aggr.df$post_cage <- stringr::str_extract(l.aggr.df$name, "\\d+$")
l.aggr.df <- l.aggr.df %>% left_join(
  data.frame(time=unlist(obs.secs), batch = stringr::str_extract(names(obs.secs), "\\d+$"))
)
l.aggr.df$post_batchcage <- gsub("\\s+", "", l.aggr.df$name)

# generate per hour value 
l.aggr.df$post.given1 <- 3600*(l.aggr.df$given / l.aggr.df$time)
l.aggr.df$post.received1 <- 3600*(l.aggr.df$received / l.aggr.df$time)

l.aggr.df$name <- NULL
l.aggr.df$post_cage <- as.numeric(l.aggr.df$post_cage)


l.aggr.df$postid <- gsub("Cage", "", l.aggr.df$postid)
l.aggr.df$post_idbatch <- paste0(l.aggr.df$postid, l.aggr.df$post_batchcage)

l.aggr.df$post_idbatch <- substr(l.aggr.df$post_idbatch,2,11) 
l.aggr.df$post_idbatch <- gsub("C", "", l.aggr.df$post_idbatch)

df <- read_csv("manuscript/brain/results_tables/coldata_ALL.csv")
head(df)
dfx <- df %>% filter(time == 25)
dfx$SampleID

ad <- l.aggr.df %>% select(8:10)
adx <- df %>% full_join(ad)
data <- adx[1:52,]

datax  <- data %>% filter(time == 25)
datax$SampleID
write.csv(data, "manuscript/brain/results_tables/coldata_ALLxAGG.csv", row.names =F)
