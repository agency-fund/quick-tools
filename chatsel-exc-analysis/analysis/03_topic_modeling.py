#!/usr/bin/env python3
"""
Topic modeling for user and AI messages using BERTopic.
Works with multilingual text (Spanish and English).
Analyzes user and AI messages separately to understand different themes.
"""

import pandas as pd
import numpy as np
from bertopic import BERTopic
from sklearn.feature_extraction.text import CountVectorizer
import warnings
warnings.filterwarnings('ignore')

# Load the data
print("Loading data...")
messages = pd.read_csv("data/messages_for_python.csv")

# Separate user and AI messages
user_messages = messages[messages['role'] == 'user'].copy()
user_messages = user_messages[user_messages['content'].notna()].reset_index(drop=True)

ai_messages = messages[messages['role'] == 'assistant'].copy()
ai_messages = ai_messages[ai_messages['content'].notna()].reset_index(drop=True)

print(f"\nUser messages: {len(user_messages)}")
print(f"AI messages: {len(ai_messages)}")

# Prepare documents
user_docs = user_messages['content'].astype(str).tolist()
ai_docs = ai_messages['content'].astype(str).tolist()

# Create custom stop words for Spanish and English
spanish_stop = [
    "el", "la", "de", "que", "y", "a", "en", "un", "ser", "se", "no",
    "haber", "por", "con", "su", "para", "como", "estar", "tener", "le",
    "lo", "todo", "pero", "más", "hacer", "o", "poder", "decir", "este",
    "ir", "otro", "ese", "si", "me", "ya", "ver", "porque", "dar", "cuando",
    "muy", "sin", "vez", "mucho", "saber", "sobre", "mi", "alguno",
    "mismo", "yo", "también", "hasta", "dos", "querer", "entre", "así",
    "desde", "grande", "eso", "ni", "nos", "llegar", "pasar", "tiempo",
    "ella", "día", "uno", "bien", "poco", "deber", "entonces", "poner",
    "cosa", "tanto", "parecer", "nuestro", "tan", "donde", "ahora",
    "parte", "después", "vida", "quedar", "siempre", "creer", "hablar",
    "llevar", "dejar", "nada", "cada", "seguir", "menos", "nuevo",
    "encontrar", "algo", "solo", "pueden", "puedo", "puede", "estoy",
    "está", "son", "hay", "he", "has", "ha", "hemos", "han", "eres",
    "es", "somos", "soy", "fue", "sido", "siendo", "era", "iba",
    "tengo", "tienes", "tiene", "tenemos", "tienen"
]

def create_topic_model():
    """Create a BERTopic model with Spanish/English support."""
    vectorizer_model = CountVectorizer(
        stop_words=spanish_stop,
        ngram_range=(1, 2),
        min_df=2
    )

    return BERTopic(
        language="multilingual",
        calculate_probabilities=False,
        vectorizer_model=vectorizer_model,
        verbose=False,
        min_topic_size=5
    )

def display_topics(topic_model, messages_df, role_name, n_topics=5):
    """Display topic information with examples."""
    topic_info = topic_model.get_topic_info()
    print(f"\n{'='*60}")
    print(f"{role_name.upper()} TOPICS (Temas en Español/Inglés)")
    print(f"{'='*60}")
    print(f"Number of topics found: {len(topic_info) - 1}")  # -1 for outliers

    print(f"\nTop topics for {role_name}:")
    for topic_id in range(min(n_topics, len(topic_info) - 1)):
        topic_words = topic_model.get_topic(topic_id)
        if topic_words:
            keywords = ', '.join([word for word, _ in topic_words[:8]])
            count = len(messages_df[messages_df['topic'] == topic_id])
            print(f"\nTopic {topic_id} ({count} messages):")
            print(f"  Keywords: {keywords}")

            # Show 2 example messages
            topic_docs = messages_df[messages_df['topic'] == topic_id]['content'].head(2)
            for i, doc in enumerate(topic_docs, 1):
                preview = doc[:120] + "..." if len(doc) > 120 else doc
                print(f"  Example {i}: {preview}")

# ============ USER MESSAGE TOPICS ============
print("\n" + "="*60)
print("ANALYZING USER MESSAGE TOPICS...")
print("="*60)

user_topic_model = create_topic_model()
print(f"Fitting topic model for {len(user_docs)} user messages...")
user_topics, _ = user_topic_model.fit_transform(user_docs)
user_messages['topic'] = user_topics

user_topic_info = user_topic_model.get_topic_info()
display_topics(user_topic_model, user_messages, "User", n_topics=8)

# Save user topic results
user_topic_info.to_csv("analysis/topic_info_user.csv", index=False)
user_messages[['id', 'topic']].to_csv("analysis/message_topics_user.csv", index=False)
user_topic_model.save("analysis/topic_model_user")

# ============ AI MESSAGE TOPICS ============
print("\n" + "="*60)
print("ANALYZING AI MESSAGE TOPICS...")
print("="*60)

ai_topic_model = create_topic_model()
print(f"Fitting topic model for {len(ai_docs)} AI messages...")
ai_topics, _ = ai_topic_model.fit_transform(ai_docs)
ai_messages['topic'] = ai_topics

ai_topic_info = ai_topic_model.get_topic_info()
display_topics(ai_topic_model, ai_messages, "AI Assistant", n_topics=8)

# Save AI topic results
ai_topic_info.to_csv("analysis/topic_info_ai.csv", index=False)
ai_messages[['id', 'topic']].to_csv("analysis/message_topics_ai.csv", index=False)
ai_topic_model.save("analysis/topic_model_ai")

print("\n" + "="*60)
print("✓ Results saved:")
print("="*60)
print("\nUser topics:")
print("  - analysis/topic_info_user.csv")
print("  - analysis/message_topics_user.csv")
print("  - analysis/topic_model_user/")
print("\nAI topics:")
print("  - analysis/topic_info_ai.csv")
print("  - analysis/message_topics_ai.csv")
print("  - analysis/topic_model_ai/")
