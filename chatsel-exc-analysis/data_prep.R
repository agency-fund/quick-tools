library(RPostgres)
library(tidyverse)

# Load credentials from .secrets file
load_secrets <- function(file = ".secrets") {
  if (!file.exists(file)) {
    stop("Secrets file not found. Please create .secrets file with your credentials.")
  }
  lines <- readLines(file)
  lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
  secrets <- list()
  for (line in lines) {
    parts <- strsplit(line, "=")[[1]]
    if (length(parts) == 2) {
      secrets[[parts[1]]] <- parts[2]
    }
  }
  secrets
}

secrets <- load_secrets()

# Connect to PostgreSQL
con <- dbConnect(
  Postgres(),
  dbname   = secrets$PGDB,
  host     = secrets$PGHOST,
  port     = as.integer(secrets$PGPORT),
  user     = secrets$PGUSER,
  password = secrets$PGPASSWORD
)

# Test connection
cat("Connected to database. Available tables:\n")
print(dbListTables(con))

# Load data into tibbles
messages_raw <- dbGetQuery(con, "SELECT * FROM messages") %>% as_tibble()
conversations_raw <- dbGetQuery(con, "SELECT * FROM conversations") %>% as_tibble()
users_raw <- dbGetQuery(con, "SELECT * FROM users") %>% as_tibble()

# Create data directory if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Save raw data as RDS files
saveRDS(messages_raw, "data/messages_raw.rds")
saveRDS(conversations_raw, "data/conversations_raw.rds")
saveRDS(users_raw, "data/users_raw.rds")

cat("\nData saved to data/ directory:\n")
cat("- messages_raw.rds:", nrow(messages_raw), "rows\n")
cat("- conversations_raw.rds:", nrow(conversations_raw), "rows\n")
cat("- users_raw.rds:", nrow(users_raw), "rows\n")

# Close connection
dbDisconnect(con)
