library(tidyverse)
source("load_data.R")

# Load data from saved RDS files
data <- load_chat_data()
messages <- data$messages
conversations <- data$conversations
users <- data$users

# Explore the data structure
cat("=== Data Overview ===\n\n")

cat("Messages:", nrow(messages), "rows\n")
cat("Columns:", paste(names(messages), collapse = ", "), "\n\n")

cat("Conversations:", nrow(conversations), "rows\n")
cat("Columns:", paste(names(conversations), collapse = ", "), "\n\n")

cat("Users:", nrow(users), "rows\n")
cat("Columns:", paste(names(users), collapse = ", "), "\n\n")

# Preview data
cat("=== Messages Preview ===\n")
print(head(messages, 5))

cat("\n=== Conversations Preview ===\n")
print(head(conversations, 5))

cat("\n=== Users Preview ===\n")
print(head(users, 10))