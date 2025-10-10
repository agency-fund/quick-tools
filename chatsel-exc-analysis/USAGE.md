# Usage Guide / Guía de Uso

## Initial Setup (Only Once) / Configuración Inicial (Solo Una Vez)

### 1. Fetch Data / Obtener Datos

```bash
Rscript data_prep.R
```

This pulls data from PostgreSQL and saves it locally in `data/`.

*Esto extrae datos de PostgreSQL y los guarda localmente en `data/`.*

### 2. Run All Analyses / Ejecutar Todos los Análisis

```bash
./run_analysis.sh
```

This runs the complete analysis pipeline (takes 5-10 minutes).

*Esto ejecuta el pipeline completo de análisis (toma 5-10 minutos).*

## Quick Report Generation / Generación Rápida de Reporte

If you've already run the analyses once and just want to regenerate the report:

*Si ya ejecutaste los análisis y solo quieres regenerar el reporte:*

```bash
quarto render chat_analysis_report.qmd
```

## Individual Analysis Scripts / Scripts Individuales de Análisis

You can run analyses separately:

*Puedes ejecutar análisis por separado:*

### User Message Summary / Resumen de Mensajes de Usuario
```bash
Rscript analysis/01_user_message_summary.R
```

### Sentiment Analysis / Análisis de Sentimiento
```bash
Rscript analysis/02_prepare_data_for_python.R
python3 analysis/02_sentiment_analysis.py
```

### Topic Modeling / Modelado de Tópicos
```bash
python3 analysis/03_topic_modeling.py
```

### AI Interaction Stats / Estadísticas de Interacción con IA
```bash
Rscript analysis/04_ai_interaction_stats.R
```

## Report Sections / Secciones del Reporte

The final HTML report includes:

*El reporte HTML final incluye:*

1. **Part 1: User Message Analysis** / *Análisis de Mensajes de Usuario*
   - Word frequency / Frecuencia de palabras
   - Message length distribution / Distribución de longitud de mensajes
   - Temporal patterns / Patrones temporales

2. **Part 2: Sentiment Analysis** / *Análisis de Sentimiento*
   - User sentiment / Sentimiento de usuario
   - AI sentiment / Sentimiento de IA
   - Comparison / Comparación

3. **Part 3: Topic Modeling** / *Modelado de Tópicos*
   - User topics / Tópicos de usuario
   - AI topics / Tópicos de IA
   - Topic distribution / Distribución de tópicos

4. **Part 4: Interaction Statistics** / *Estadísticas de Interacción*
   - Conversation metrics / Métricas de conversación
   - Response characteristics / Características de respuesta
   - Usage patterns / Patrones de uso

## Notes / Notas

- **Spanish Language**: The chatbot data is primarily in Spanish. The report uses Spanish stopwords and handles Spanish text appropriately.

  *El lenguaje de los datos del chatbot es principalmente español. El reporte usa palabras vacías en español y maneja el texto en español apropiadamente.*

- **Bilingual Output**: All report titles, labels, and summaries are in both Spanish and English.

  *Todas las titles, etiquetas y resúmenes del reporte están en español e inglés.*

- **User vs AI Analysis**: Sentiment and topic modeling are performed separately for user messages and AI responses to understand communication patterns from both sides.

  *El análisis de sentimiento y el modelado de tópicos se realizan por separado para mensajes de usuario y respuestas de IA para entender los patrones de comunicación de ambos lados.*
