library(tidyverse)
library(tidytext)

source("load_data.R")

# Load data
messages <- load_messages()

# Filter user messages only
user_messages <- messages %>%
  filter(role == "user") %>%
  mutate(
    msg_length = nchar(as.character(content)),
    created_date = as.Date(created_at)
  )

# Basic summary statistics
summary_stats <- list(
  total_messages = nrow(user_messages),
  avg_length = mean(user_messages$msg_length, na.rm = TRUE),
  median_length = median(user_messages$msg_length, na.rm = TRUE),
  unique_conversations = n_distinct(user_messages$conversation_id),
  date_range = range(user_messages$created_date)
)

cat("=== User Message Summary ===\n")
cat("Total user messages:", summary_stats$total_messages, "\n")
cat("Average message length:", round(summary_stats$avg_length, 1), "characters\n")
cat("Median message length:", summary_stats$median_length, "characters\n")
cat("Unique conversations:", summary_stats$unique_conversations, "\n")
cat("Date range:", format(summary_stats$date_range[1]), "to", format(summary_stats$date_range[2]), "\n\n")

# Tokenize and get most common words
# For Spanish text, we'll use Spanish stop words
spanish_stopwords <- data.frame(
  word = c("el", "la", "de", "que", "y", "a", "en", "un", "ser", "se", "no",
           "haber", "por", "con", "su", "para", "como", "estar", "tener", "le",
           "lo", "todo", "pero", "más", "hacer", "o", "poder", "decir", "este",
           "ir", "otro", "ese", "si", "me", "ya", "ver", "porque", "dar", "cuando",
           "él", "muy", "sin", "vez", "mucho", "saber", "qué", "sobre", "mi", "alguno",
           "mismo", "yo", "también", "hasta", "año", "dos", "querer", "entre", "así",
           "primero", "desde", "grande", "eso", "ni", "nos", "llegar", "pasar", "tiempo",
           "ella", "sí", "día", "uno", "bien", "poco", "deber", "entonces", "poner",
           "cosa", "tanto", "hombre", "parecer", "nuestro", "tan", "donde", "ahora",
           "parte", "después", "vida", "quedar", "siempre", "creer", "hablar", "llevar",
           "dejar", "nada", "cada", "seguir", "menos", "nuevo", "encontrar", "algo",
           "solo", "decir", "pueden", "puedo", "puede", "estoy", "está", "son", "hay",
           "he", "has", "ha", "hemos", "han", "eres", "es", "somos", "soy", "fue",
           "sido", "siendo", "era", "iba", "tengo", "tienes", "tiene", "tenemos", "tienen")
)

user_words <- user_messages %>%
  select(id, content) %>%
  unnest_tokens(word, content) %>%
  anti_join(spanish_stopwords, by = "word") %>%
  anti_join(stop_words, by = "word") %>%  # English stopwords
  filter(!str_detect(word, "^\\d+$")) %>%  # Remove pure numbers
  count(word, sort = TRUE)

cat("=== Top 30 Most Common Words in User Messages ===\n")
print(head(user_words, 30))

# Messages per conversation
msgs_per_conv <- user_messages %>%
  count(conversation_id, name = "user_msgs_count") %>%
  summarise(
    mean_msgs = mean(user_msgs_count),
    median_msgs = median(user_msgs_count),
    max_msgs = max(user_msgs_count)
  )

cat("\n=== User Messages per Conversation ===\n")
cat("Mean:", round(msgs_per_conv$mean_msgs, 1), "\n")
cat("Median:", msgs_per_conv$median_msgs, "\n")
cat("Max:", msgs_per_conv$max_msgs, "\n")

# Save results
results <- list(
  summary_stats = summary_stats,
  top_words = user_words,
  msgs_per_conv = msgs_per_conv,
  user_messages = user_messages
)

saveRDS(results, "analysis/user_message_summary.rds")
cat("\nResults saved to analysis/user_message_summary.rds\n")
