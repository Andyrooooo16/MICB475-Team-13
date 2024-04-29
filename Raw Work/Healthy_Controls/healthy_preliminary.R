file_path <- "../halfvarson_metadata.tsv"
data <- read.table(file_path, sep = "\t", header = TRUE)



# Filter rows where the column "timepoint" contains the value 1
# Subset the data to include only rows where timepoint equals 1
timepoint1_data <- data[data$timepoint == 1, ]

# View the unique values of timepoint in timepoint1_data
print(unique(timepoint1_data$timepoint))

# Subset timepoint1_data based on the condition cd_resection != "not applicable"
timepoint1_and_cd_resection_data <- timepoint1_data[timepoint1_data$cd_resection != "not applicable", ]

# View the unique values of timepoint in timepoint1_and_cd_resection_data
print(unique(timepoint1_and_cd_resection_data$cd_resection))

# Export the filtered dataset to a TSV file
write.table(timepoint1_and_cd_resection_data, file = "C:/Users/andre/Downloads/timepoint1_and_cd_resection_data.tsv", sep = "\t", row.names = FALSE)




#calprotectin datset & stats

calprotectin_data <- subset(cd_resection_data, !(calprotectin %in% c("not applicable", "not collected")))


# Calculate statistics for the "calprotectin" column
calprotectin_stats <- c(
  mean = mean(as.numeric(calprotectin_data$calprotectin)),
  median = median(as.numeric(calprotectin_data$calprotectin)),
  sd = sd(as.numeric(calprotectin_data$calprotectin)),
  min = min(as.numeric(calprotectin_data$calprotectin)),
  max = max(as.numeric(calprotectin_data$calprotectin))
)

# Print the results with labels
print(calprotectin_stats)




