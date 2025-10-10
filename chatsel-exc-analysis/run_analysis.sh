#!/bin/bash
# Script to run all analyses and generate the report

set -e  # Exit on error

echo "=========================================="
echo "Chat Analysis Pipeline"
echo "=========================================="

# 1. User Message Summary (R)
echo ""
echo "[1/5] Running user message summary analysis..."
Rscript analysis/01_user_message_summary.R

# 2. Prepare data for Python
echo ""
echo "[2/5] Preparing data for Python analyses..."
Rscript analysis/02_prepare_data_for_python.R

# 3. Sentiment Analysis (Python)
echo ""
echo "[3/5] Running sentiment analysis..."
echo "  (This may take a few minutes to download models on first run)"
python3 analysis/02_sentiment_analysis.py

# 4. Topic Modeling (Python)
echo ""
echo "[4/5] Running topic modeling..."
echo "  (This may take several minutes)"
python3 analysis/03_topic_modeling.py

# 5. AI Interaction Stats (R)
echo ""
echo "[5/5] Running AI interaction statistics..."
Rscript analysis/04_ai_interaction_stats.R

# 6. Generate Quarto report
echo ""
echo "=========================================="
echo "Rendering Quarto report..."
echo "=========================================="
quarto render chat_analysis_report.qmd

echo ""
echo "=========================================="
echo "✓ Analysis complete!"
echo "=========================================="
echo ""
echo "Results:"
echo "  - HTML report: chat_analysis_report.html"
echo "  - Analysis outputs: analysis/"
echo ""
