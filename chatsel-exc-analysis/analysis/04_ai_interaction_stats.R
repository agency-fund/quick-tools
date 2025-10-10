library(tidyverse)
library(lubridate)

source("load_data.R")

# Load data
messages <- load_messages()
conversations <- load_conversations()

cat("=== AI Interaction Statistics ===\n\n")

# Overall statistics
total_messages <- nrow(messages)
total_conversations <- n_distinct(messages$conversation_id)
ai_messages <- messages %>% filter(role == "assistant")
user_messages <- messages %>% filter(role == "user")

cat("## Overall Statistics\n")
cat("Total messages:", total_messages, "\n")
cat("Total conversations:", total_conversations, "\n")
cat("AI messages:", nrow(ai_messages), "\n")
cat("User messages:", nrow(user_messages), "\n")
cat("Messages per conversation (avg):", round(total_messages / total_conversations, 1), "\n\n")

# Conversation turns analysis
conv_turns <- messages %>%
  arrange(conversation_id, created_at) %>%
  group_by(conversation_id) %>%
  summarise(
    total_messages = n(),
    user_messages = sum(role == "user"),
    ai_messages = sum(role == "assistant"),
    turns = ceiling(total_messages / 2),  # Approximate turn count
    duration_minutes = as.numeric(difftime(max(created_at), min(created_at), units = "mins")),
    .groups = "drop"
  )

cat("## Conversation Turns\n")
cat("Average turns per conversation:", round(mean(conv_turns$turns), 1), "\n")
cat("Median turns per conversation:", median(conv_turns$turns), "\n")
cat("Max turns in a conversation:", max(conv_turns$turns), "\n")
cat("Average conversation duration:", round(mean(conv_turns$duration_minutes, na.rm = TRUE), 1), "minutes\n\n")

# AI response characteristics
ai_stats <- ai_messages %>%
  mutate(
    response_length = nchar(as.character(content)),
    has_tool_call = !is.na(was_tool_call) & was_tool_call == TRUE,
    created_date = as.Date(created_at),
    created_hour = hour(created_at)
  )

cat("## AI Response Characteristics\n")
cat("Average AI response length:", round(mean(ai_stats$response_length, na.rm = TRUE), 0), "characters\n")
cat("Median AI response length:", median(ai_stats$response_length, na.rm = TRUE), "characters\n")
cat("Responses with tool calls:", sum(ai_stats$has_tool_call, na.rm = TRUE),
    sprintf("(%.1f%%)", 100 * mean(ai_stats$has_tool_call, na.rm = TRUE)), "\n\n")

# Response time analysis (if sequential messages)
response_times <- messages %>%
  arrange(conversation_id, created_at) %>%
  group_by(conversation_id) %>%
  mutate(
    prev_role = lag(role),
    time_diff = as.numeric(difftime(created_at, lag(created_at), units = "secs"))
  ) %>%
  filter(role == "assistant" & prev_role == "user") %>%
  ungroup()

if (nrow(response_times) > 0) {
  cat("## AI Response Time\n")
  cat("Average response time:", round(mean(response_times$time_diff, na.rm = TRUE), 1), "seconds\n")
  cat("Median response time:", round(median(response_times$time_diff, na.rm = TRUE), 1), "seconds\n\n")
}

# Temporal patterns
temporal_stats <- messages %>%
  mutate(
    date = as.Date(created_at),
    hour = hour(created_at),
    weekday = wday(created_at, label = TRUE)
  ) %>%
  group_by(date) %>%
  summarise(
    messages_per_day = n(),
    conversations_per_day = n_distinct(conversation_id),
    .groups = "drop"
  )

cat("## Temporal Patterns\n")
cat("Active days:", nrow(temporal_stats), "\n")
cat("Average messages per active day:", round(mean(temporal_stats$messages_per_day), 1), "\n")
cat("Peak day:", temporal_stats$date[which.max(temporal_stats$messages_per_day)],
    "with", max(temporal_stats$messages_per_day), "messages\n\n")

# Most active hours
hour_stats <- messages %>%
  mutate(hour = hour(created_at)) %>%
  count(hour) %>%
  arrange(desc(n))

cat("## Most Active Hours\n")
print(head(hour_stats, 5))
cat("\n")

# Conversation length distribution
conv_length_dist <- conv_turns %>%
  mutate(
    length_category = case_when(
      total_messages <= 2 ~ "1-2 messages",
      total_messages <= 5 ~ "3-5 messages",
      total_messages <= 10 ~ "6-10 messages",
      total_messages <= 20 ~ "11-20 messages",
      TRUE ~ "20+ messages"
    )
  ) %>%
  count(length_category)

cat("## Conversation Length Distribution\n")
print(conv_length_dist)

# Save results
results <- list(
  overall_stats = list(
    total_messages = total_messages,
    total_conversations = total_conversations,
    ai_messages_count = nrow(ai_messages),
    user_messages_count = nrow(user_messages)
  ),
  conv_turns = conv_turns,
  ai_stats = ai_stats,
  response_times = response_times,
  temporal_stats = temporal_stats,
  hour_stats = hour_stats,
  conv_length_dist = conv_length_dist
)

saveRDS(results, "analysis/ai_interaction_stats.rds")
cat("\n✓ Results saved to analysis/ai_interaction_stats.rds\n")
