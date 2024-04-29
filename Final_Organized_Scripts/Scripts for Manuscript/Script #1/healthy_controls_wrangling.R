### Importing Packages ###
library(tidyverse)

# Loading in the dataset
file_path <- "halfvarson_metadata.tsv"
data <- read.delim(file_path, sep = "\t", header = TRUE)

# Filter rows where the column "timepoint" contains the value 1
# Subset the data to include only rows where timepoint equals 1
wrangled <- filter(data, timepoint == 1)

# View the unique values of timepoint in timepoint1_data
print(unique(wrangled$timepoint))

# Subset timepoint1_data based on the condition cd_resection != "not applicable"
wrangled2 <- filter(wrangled, !(cd_resection == "not applicable" & diagnosis_full != "HC"))

# View the unique values of timepoint in timepoint1_and_cd_resection_data
print(unique(wrangled2$cd_resection))

# converts remaining calprotectin values to numerics
wrangled2$calprotectin <- as.numeric(wrangled2$calprotectin)

# Determine the number of rows that are below our inflammation cutoff (calprotectin <= 150)
no_inflammation <- wrangled2 %>%
    filter(calprotectin<= 150) 

count_no_inflammation <- no_inflammation %>%
    nrow()

# Determine the number of rows that are above our inflammation cutoff (calprotect > 150)
yes_inflammation <- wrangled2 %>%
  filter(calprotectin > 150) 

count_yes_inflammation <- yes_inflammation %>%
  nrow()

# Display number of samples w/in the 57 filtered rows that don't have inflammation
count_no_inflammation
count_yes_inflammation

# Any samples that have a calprotectin level over 150 will be considered to have inflammation.
# Any samples with a calprotectin level below or equal to 150 will not be considered to have inflammation.
# We will add a column called 'inflammation' which will be a logical.
# We will also add a column called disease severity where samples will be assigned a value of 1, 2 or 3
# depending on whether their cd_behaviour is classified as either B1, B2, B3

mutated_calprotectin_data <- wrangled2 %>%
      mutate(inflammation = ifelse(calprotectin > 150, TRUE, FALSE)) %>% # Adds 'inflammation' column based on calprotectin level
      mutate(disease_severity = case_when(
          cd_behavior == "Non-stricturing, non-penetrating (B1)" ~ "Low",
          cd_behavior == "Stricturing (B2)" ~ "Medium",
          cd_behavior == "Penetrating (B3)" ~ "High"
      )) # Adds 'disease_severity' column depending on score in 'cd_behavior') 
# mutated_calprotectin_data

# Adds a column that determines whether the sample has inflammation and had a resection
data_final <- mutated_calprotectin_data %>%
      mutate(inflammation_with_surgery = case_when(
        inflammation == TRUE & cd_resection == "yes" ~ "inflammation_with_surgery",
        inflammation == TRUE & cd_resection == "no" ~ "inflammation_no_surgery",
        inflammation == FALSE & cd_resection == "yes" ~ "no_inflammation_with_surgery",
        TRUE ~ "no_inflammation_no_surgery"
      ))
# data_final 



# Export the filtered dataset to a TSV file
write.table(data_final, file = "hc_halfvarson_metadata_wrangled.tsv", sep = "\t", row.names = FALSE)

