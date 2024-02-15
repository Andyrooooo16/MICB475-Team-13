### Importing Packages ###
library(tidyverse)

# Loading in the dataset
file_path <- "/Users/Max/Desktop/MICB475_group_project_personal/halfvarson_metadata.tsv"
data <- read.delim(file_path, sep = "\t", header = TRUE)



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
# write.table(timepoint1_and_cd_resection_data, file = "timepoint1_and_cd_resection_data.tsv", sep = "\t", row.names = FALSE)




#calprotectin datset & stats

calprotectin_data <- subset(timepoint1_and_cd_resection_data, !(calprotectin %in% c("not applicable", "not collected")))


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


# Determine the number of rows that are below our inflammation cutoff (calprotectin <= 150)
no_inflammation <- calprotectin_data %>%
    filter(as.integer(calprotectin) <= 150) 

count_no_inflammation <- no_inflammation %>%
    nrow()

# Determine the number of rows that are above our inflammation cutoff (calprotect > 150)
yes_inflammation <- calprotectin_data %>%
  filter(as.integer(calprotectin) > 150) 

count_yes_inflammation <- yes_inflammation %>%
  nrow()

# Display number of samples w/in the 48 filtered rows that don't have inflammation
count_no_inflammation
count_yes_inflammation

# Any samples that have a calprotectin level over 150 will be considered to have inflammation.
# Any samples with a calprotectin level below or equal to 150 will not be considered to have inflammation.
# We will add a column called 'inflammation' which will be a logical.
# We will also add a column called disease severity where samples will be assigned a value of 1, 2 or 3
# depending on whether their cd_behaviour is classified as either B1, B2, B3

mutated_calprotectin_data <- calprotectin_data %>%
      mutate(inflammation = ifelse(as.integer(calprotectin) > 150, TRUE, FALSE)) %>% # Adds 'inflammation' column based on calprotectin level
      mutate(disease_severity = case_when(
          cd_behavior == "Non-stricturing, non-penetrating (B1)" ~ "Low",
          cd_behavior == "Stricturing (B2)" ~ "Medium",
          cd_behavior == "Penetrating (B3)" ~ "High"
      )) # Adds 'disease_severity' column depending on score in 'cd_behavior'

mutated_calprotectin_data


