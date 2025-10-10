library(tidyverse)

# Load raw data from RDS files
load_chat_data <- function() {
  if (!file.exists("data/messages_raw.rds")) {
    stop("Data files not found. Please run data_prep.R first to fetch data from database.")
  }

  list(
    messages = readRDS("data/messages_raw.rds"),
    conversations = readRDS("data/conversations_raw.rds"),
    users = readRDS("data/users_raw.rds")
  )
}

# Convenience function to load individual tables
load_messages <- function() {
  readRDS("data/messages_raw.rds")
}

load_conversations <- function() {
  readRDS("data/conversations_raw.rds")
}

load_users <- function() {
  readRDS("data/users_raw.rds")
}
