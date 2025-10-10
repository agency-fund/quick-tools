# Chat Selection Exception Analysis

Analysis of chat data from AI bot stored in PostgreSQL. This repo performs comprehensive text analysis including sentiment analysis, topic modeling, and interaction statistics.

## Setup

### 1. Install R Packages

```r
install.packages(c("tidyverse", "tidytext", "lubridate", "scales", "knitr"))
```

### 2. Install Python Dependencies

```bash
pip install -r requirements.txt
```

### 3. Install Quarto

Download and install Quarto from [https://quarto.org](https://quarto.org)

## Quick Start

### Fetch Data from Database

```bash
Rscript data_prep.R
```

This connects to PostgreSQL (credentials in [.secrets](.secrets), gitignored) and saves data to `data/`.

### Run All Analyses and Generate Report

```bash
./run_analysis.sh
```

This will:
1. Analyze user messages (word frequency, patterns)
2. Perform sentiment analysis (using multilingual transformers)
3. Run topic modeling (BERTopic)
4. Calculate AI interaction statistics
5. Generate HTML report ([chat_analysis_report.html](chat_analysis_report.html))

## Files

### Data Preparation
- **[data_prep.R](data_prep.R)**: Connects to PostgreSQL and saves data locally
- **[load_data.R](load_data.R)**: Helper functions to load saved data
- **[.secrets](.secrets)**: Database credentials (gitignored)

### Analysis Scripts (in [analysis/](analysis/))
1. **[01_user_message_summary.R](analysis/01_user_message_summary.R)**: User message analysis with word frequency
2. **[02_prepare_data_for_python.R](analysis/02_prepare_data_for_python.R)**: Export data for Python scripts
3. **[02_sentiment_analysis.py](analysis/02_sentiment_analysis.py)**: Multilingual sentiment analysis
4. **[03_topic_modeling.py](analysis/03_topic_modeling.py)**: BERTopic for thematic analysis
5. **[04_ai_interaction_stats.R](analysis/04_ai_interaction_stats.R)**: Conversation turns and AI response stats

### Reports
- **[chat_analysis_report.qmd](chat_analysis_report.qmd)**: Quarto document combining all analyses
- **[analysis.R](analysis.R)**: Quick data exploration script

## Analysis Overview

### User Messages
- Word frequency analysis (Spanish + English)
- Message length distribution
- Temporal patterns
- Messages per conversation

### Sentiment Analysis
- Multilingual sentiment detection
- Positive/Negative/Neutral classification
- Sentiment score distributions

### Topic Modeling
- Automatic topic discovery using BERTopic
- Theme identification across messages
- Topic distribution and keywords

### AI Interaction Statistics
- Conversation turns and duration
- AI response characteristics
- Response time analysis
- Hourly and daily usage patterns

## Data Structure

The database contains three main tables:
- **messages**: Individual chat messages (role, content, timestamps)
- **conversations**: Conversation metadata
- **users**: User information
