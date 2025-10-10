#!/usr/bin/env python3
"""
Sentiment analysis for user and AI messages using transformers.
Uses a multilingual sentiment model that works for both Spanish and English.
Analyzes user and AI messages separately to compare sentiment patterns.
"""

import pandas as pd
import numpy as np
from transformers import pipeline
import warnings
warnings.filterwarnings('ignore')

def analyze_sentiment(messages_df, role_name):
    """Analyze sentiment for messages of a specific role."""
    print(f"\nAnalyzing {len(messages_df)} {role_name} messages...")

    # Load sentiment analysis model (multilingual)
    # Using nlptown/bert-base-multilingual-uncased-sentiment for multilingual support
    sentiment_analyzer = pipeline(
        "sentiment-analysis",
        model="nlptown/bert-base-multilingual-uncased-sentiment",
        max_length=512,
        truncation=True
    )

    sentiments = []
    for idx, row in messages_df.iterrows():
        if idx % 50 == 0:
            print(f"  Processed {idx}/{len(messages_df)} {role_name} messages...")

        try:
            result = sentiment_analyzer(str(row['content'])[:512])[0]
            sentiments.append({
                'id': row['id'],
                'role': role_name,
                'sentiment_label': result['label'],
                'sentiment_score': result['score']
            })
        except Exception as e:
            print(f"Error processing message {idx}: {e}")
            sentiments.append({
                'id': row['id'],
                'role': role_name,
                'sentiment_label': 'neutral',
                'sentiment_score': 0.33
            })

    return pd.DataFrame(sentiments)

# Load the data
print("Loading data...")
messages = pd.read_csv("data/messages_for_python.csv")

# Separate user and AI messages
user_messages = messages[messages['role'] == 'user'].copy()
ai_messages = messages[messages['role'] == 'assistant'].copy()

print(f"\nTotal messages: {len(messages)}")
print(f"User messages: {len(user_messages)}")
print(f"AI messages: {len(ai_messages)}")

# Load sentiment analysis model once
print("\nLoading sentiment analysis model...")

# Analyze user messages
user_sentiment = analyze_sentiment(user_messages, 'user')

# Analyze AI messages
ai_sentiment = analyze_sentiment(ai_messages, 'assistant')

# Combine results
all_sentiment = pd.concat([user_sentiment, ai_sentiment], ignore_index=True)

# Summary statistics
print("\n" + "="*60)
print("SENTIMENT ANALYSIS RESULTS (Spanish/English text)")
print("="*60)

print("\n--- USER MESSAGES ---")
print(f"Total analyzed: {len(user_sentiment)}")
print("\nSentiment distribution (Usuario):")
print(user_sentiment['sentiment_label'].value_counts())
print("\nAverage confidence scores by sentiment:")
print(user_sentiment.groupby('sentiment_label')['sentiment_score'].agg(['mean', 'median', 'count']))

print("\n--- AI MESSAGES ---")
print(f"Total analyzed: {len(ai_sentiment)}")
print("\nSentiment distribution (Asistente):")
print(ai_sentiment['sentiment_label'].value_counts())
print("\nAverage confidence scores by sentiment:")
print(ai_sentiment.groupby('sentiment_label')['sentiment_score'].agg(['mean', 'median', 'count']))

# Save results
user_sentiment.to_csv("analysis/sentiment_user.csv", index=False)
ai_sentiment.to_csv("analysis/sentiment_ai.csv", index=False)
all_sentiment.to_csv("analysis/sentiment_all.csv", index=False)

print("\n✓ Results saved:")
print("  - analysis/sentiment_user.csv")
print("  - analysis/sentiment_ai.csv")
print("  - analysis/sentiment_all.csv")
