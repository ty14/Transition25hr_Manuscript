
### total aggression

# import start end times of pre-reorganization behavior.

source("manuscript/behavior/01_pre_startend_times.R")

# from pre_startend_times
obs.secs <- lapply(lapply(sends,get_time_dif),function(x) as.numeric(sum(x)))


# get preorganization behavior data
df1 <- read_csv("behavior/Pre_WL.csv")

# create list of dataframes for each individual cage
l1 <- split(df1, df1$pre_batchcage)

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
  names(sum_df) <- c("preid", "given", "received")  # Rename columns
  return(sum_df)
}

l.aggr <- lapply(l1.mat, melt_matrix)
l.aggr.df <- data.table::rbindlist(Map(cbind, l.aggr, name = names(l.aggr)))
l.aggr.df$batch <- stringr::str_extract(l.aggr.df$name, "\\d+")
l.aggr.df$pre_cage <- stringr::str_extract(l.aggr.df$name, "\\d+$")
l.aggr.df <- l.aggr.df %>% left_join(
  data.frame(time=unlist(obs.secs), batch = stringr::str_extract(names(obs.secs), "\\d+$"))
)
l.aggr.df$pre_idbatchcage <- gsub("\\s+", "", l.aggr.df$name)

# generate per hour value 
l.aggr.df$given1 <- 3600*(l.aggr.df$given / l.aggr.df$time)
l.aggr.df$received1 <- 3600*(l.aggr.df$received / l.aggr.df$time)

l.aggr.df$batch <- NULL
l.aggr.df$name <- NULL
l.aggr.df$pre_cage <- as.numeric(l.aggr.df$pre_cage)
l.aggr.df$pre_idbatchcage <- paste0(l.aggr.df$preid, l.aggr.df$pre_idbatchcage)



#join with coldata to use in WGCNA clustering. 

coldata <- read.csv("manuscript/brain/results_tables/coldata.csv", row.names = 1)
# get rid of controls
coldata <- coldata %>%  select(-AggGiven70min, -AggRec70min)
str(coldata)
str(l.aggr.df)

agg_df <- l.aggr.df %>% select(pre_idbatchcage, given, received, given1,received1)
coldata$pre_idbatchcage
post.df <- coldata %>% full_join(agg_df) %>% na.omit(.)

write.csv(post.df, "manuscript/brain/results_tables/coldata.csv", row.names = F)
