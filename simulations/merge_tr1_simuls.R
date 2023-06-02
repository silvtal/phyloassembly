# Load required libraries
library(tidyverse)


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Three positional arguments are required: output_filename, directory, and the core+sample name")
}
output_filename <- args[1]
directory <- args[2]
i <- args[3]

# Get a list of all CSV files in the directory that end with .csv and contain "i_sa"
csv_files <- list.files(directory, pattern = paste0(i, "_sa.*\\.csv$"), full.names = TRUE)
print(csv_files)
combined_data <- data.frame()
# Loop
for (file in csv_files) {
  data <- read.csv(file, stringsAsFactors = FALSE, check.names = F, row.names = 1)
  headers <- unique(names(data))
  # Add any new headers to the combined data frame
  if (nrow(combined_data) > 0) {
    new_headers <- setdiff(headers, names(combined_data))
    for (header in new_headers) {
      combined_data[[header]] <- NA
    }
  } 
  
  combined_data <- bind_rows(combined_data, data, .id = NULL) %>% replace(is.na(.), 0)
  
}
if ("dummy" %in% colnames(combined_data)) {
  combined_data <- subset(combined_data, select = -c(dummy)) # remove dummy column !
}
write.csv(combined_data, output_filename, row.names = TRUE)

