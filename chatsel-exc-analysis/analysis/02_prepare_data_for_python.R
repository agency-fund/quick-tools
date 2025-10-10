library(tidyverse)
source("load_data.R")

# Load messages
messages <- load_messages()

# Save as CSV for Python scripts to use
write_csv(messages, "data/messages_for_python.csv")
cat("Data exported to data/messages_for_python.csv\n")
