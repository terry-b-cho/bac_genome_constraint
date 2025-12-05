#!/usr/bin/env python3
"""
Script 05: Environment Prediction from GO Terms and Genome Size

Predict environment (one of 8 categories) using:
- 99% prevalence GO terms (208 terms)
- Genome size (genes_total)
- Optional normalized features

Uses supervised classification (Random Forest, Gradient Boosting, Logistic Regression).
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import argparse
# Add path to import prevalence_utils
sys.path.insert(0, str(Path(__file__).parent.parent / "4_statistical_analyses"))
from prevalence_utils import get_prevalence_prefix, filter_go_columns_by_prevalence

# Try to import sklearn
try:
    from sklearn.model_selection import train_test_split, StratifiedKFold
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler, MinMaxScaler, LabelEncoder
    from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                                classification_report, confusion_matrix,
                                precision_recall_fscore_support, roc_auc_score,
                                roc_curve, auc)
    from sklearn.multiclass import OneVsRestClassifier
    from sklearn.dummy import DummyClassifier
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    print("WARNING: sklearn not available. Install with: pip install scikit-learn")

# Try to import XGBoost with GPU support (optional)
try:
    import xgboost as xgb
    HAS_XGBOOST = True
    # Check if GPU is available
    try:
        import os
        if os.environ.get('CUDA_VISIBLE_DEVICES') is not None:
            GPU_AVAILABLE = True
        else:
            GPU_AVAILABLE = False
    except:
        GPU_AVAILABLE = False
except ImportError:
    HAS_XGBOOST = False
    GPU_AVAILABLE = False

# Try to import plotting libraries
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend for O2
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("WARNING: matplotlib/seaborn not available. Install with: pip install matplotlib seaborn")

# Try to import joblib for model saving
try:
    import joblib
    HAS_JOBLIB = True
except ImportError:
    try:
        import pickle
        HAS_JOBLIB = False
        HAS_PICKLE = True
    except ImportError:
        HAS_PICKLE = False
        print("WARNING: Neither joblib nor pickle available for model saving")

# Parse arguments
parser = argparse.ArgumentParser(description='Predict environment from GO terms and genome size')
parser.add_argument('--prevalence-threshold', type=float, default=99,
                    help='Prevalence threshold (default: 99)')
parser.add_argument('--test-mode', action='store_true',
                    help='Run on small test subset')
parser.add_argument('--normalize', type=str, choices=['none', 'per_gene', 'log', 'both'],
                    default='none', help='Feature normalization method')
parser.add_argument('--model', type=str, choices=['all', 'rf', 'gb', 'lr', 'xgb', 'baseline'],
                    default='all', help='Which model(s) to train')
parser.add_argument('--output-dir', type=str, default=None,
                    help='Output directory (default: results/5_environment_prediction/)')
parser.add_argument('--use-gpu', action='store_true',
                    help='Use GPU for XGBoost (requires CUDA and XGBoost with GPU support)')
parser.add_argument('--plot', action='store_true', default=True,
                    help='Generate visualization plots (default: True)')
args = parser.parse_args()

# Define paths
BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")
prefix = get_prevalence_prefix(args.prevalence_threshold)

# Input files
MASTER_TABLE_FILE = BASE_DIR / "results/4_statistical_analyses/02_env_cohorts" / f"{prefix}master_table_env_filtered.parquet"
MASTER_TABLE_TSV = BASE_DIR / "results/4_statistical_analyses/02_env_cohorts" / f"{prefix}master_table_env_filtered.tsv"
VALID_ENV_FILE = BASE_DIR / "results/4_statistical_analyses/02_env_cohorts" / f"{prefix}valid_environments_min20.tsv"
GO_LABELS_FILE = BASE_DIR / "results/4_statistical_analyses/05_go_labels/go_term_labels_for_plots.tsv"

# Output directory
if args.output_dir is None:
    OUTPUT_DIR = BASE_DIR / "results/5_environment_prediction"
else:
    OUTPUT_DIR = Path(args.output_dir)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

QC_LOG_FILE = OUTPUT_DIR / f"{prefix}qc_05_env_prediction.log"

# Initialize QC log
qc_log = []

def log_message(msg):
    """Add message to QC log and print to stdout."""
    qc_log.append(msg)
    print(msg)

log_message("=" * 80)
log_message("Script 05: Environment Prediction from GO Terms and Genome Size")
log_message("")
log_message("O2 Module Requirements:")
log_message("  module load conda/miniforge3/24.11.3-0")
log_message("  module load gcc/14.2.0")
if args.use_gpu or args.model in ['all', 'xgb']:
    log_message("  module load cuda/12.8  # For GPU support")
log_message("")
if args.test_mode:
    log_message("  TEST MODE: Using small subset")
log_message(f"  PREVALENCE THRESHOLD: {args.prevalence_threshold}%")
log_message(f"  NORMALIZATION: {args.normalize}")
log_message(f"  MODELS: {args.model}")
if args.use_gpu:
    log_message(f"  GPU SUPPORT: {'Enabled' if GPU_AVAILABLE else 'Requested but not available'}")
if not HAS_SKLEARN:
    log_message("  ✗ ERROR: sklearn required")
    sys.exit(1)
if not HAS_MATPLOTLIB and args.plot:
    log_message("  ⚠ WARNING: matplotlib/seaborn not available, plotting disabled")
    args.plot = False
log_message("=" * 80)
log_message("")

# ============================================================================
# 1. Data Loading and Preparation
# ============================================================================

log_message("1. Data Loading and Preparation")
log_message("-" * 80)

# Load master table
log_message("Loading master table...")
try:
    if MASTER_TABLE_FILE.exists():
        df = pd.read_parquet(MASTER_TABLE_FILE, engine='pyarrow')
        log_message(f"  ✓ Loaded {len(df)} genomes from Parquet")
    elif MASTER_TABLE_TSV.exists():
        df = pd.read_csv(MASTER_TABLE_TSV, sep='\t')
        log_message(f"  ✓ Loaded {len(df)} genomes from TSV")
    else:
        log_message(f"  ✗ ERROR: Master table not found")
        log_message(f"    Expected: {MASTER_TABLE_FILE} or {MASTER_TABLE_TSV}")
        sys.exit(1)
except Exception as e:
    log_message(f"  ✗ ERROR: Failed to load master table: {e}")
    sys.exit(1)

# Load valid environments
log_message("Loading valid environments...")
try:
    valid_envs = pd.read_csv(VALID_ENV_FILE, sep='\t')
    valid_env_list = valid_envs['environment'].tolist()
    log_message(f"  ✓ Loaded {len(valid_env_list)} valid environments")
    log_message(f"    Environments: {', '.join(valid_env_list)}")
except Exception as e:
    log_message(f"  ✗ ERROR: Failed to load valid environments: {e}")
    sys.exit(1)

# Filter to valid environments
n_before = len(df)
df = df[df['environment'].isin(valid_env_list)].copy()
n_after = len(df)
log_message(f"  Filtered to valid environments: {n_before} → {n_after} genomes")

# Verify class balance
log_message("  Class distribution:")
env_counts = df['environment'].value_counts().sort_values(ascending=False)
for env, count in env_counts.items():
    pct = 100 * count / len(df)
    log_message(f"    {env}: {count} genomes ({pct:.1f}%)")

# Test mode: subsample
if args.test_mode:
    log_message("  TEST MODE: Subsampling to 50 genomes per environment...")
    df = df.groupby('environment', group_keys=False).apply(
        lambda x: x.sample(min(50, len(x)), random_state=42)
    ).reset_index(drop=True)
    log_message(f"  After subsampling: {len(df)} genomes")

# Check for missing values in key columns
missing_env = df['environment'].isna().sum()
missing_genes = df['genes_total'].isna().sum()
if missing_env > 0:
    log_message(f"  ✗ WARNING: {missing_env} genomes with missing environment")
if missing_genes > 0:
    log_message(f"  ✗ WARNING: {missing_genes} genomes with missing genes_total")
    df = df[df['genes_total'].notna()].copy()

log_message("")

# ============================================================================
# 2. Feature Engineering
# ============================================================================

log_message("2. Feature Engineering")
log_message("-" * 80)

# Identify GO columns
go_cols = filter_go_columns_by_prevalence(df, args.prevalence_threshold)
log_message(f"  Found {len(go_cols)} GO term columns")

# Check for missing values in GO columns
go_missing = df[go_cols].isna().sum().sum()
if go_missing > 0:
    log_message(f"  ⚠ WARNING: {go_missing} missing values in GO columns (filling with 0)")
    df[go_cols] = df[go_cols].fillna(0)

# Base features: GO terms + genes_total
feature_cols = go_cols + ['genes_total']
log_message(f"  Base features: {len(feature_cols)} ({len(go_cols)} GO terms + genes_total)")

# Optional normalizations
if args.normalize in ['per_gene', 'both']:
    log_message("  Adding per-gene normalized features...")
    for col in go_cols:
        df[f"{col}_per_gene"] = df[col] / (df['genes_total'] + 1)  # +1 to avoid division by zero
    feature_cols += [f"{col}_per_gene" for col in go_cols]
    log_message(f"    Added {len(go_cols)} per-gene normalized features")

if args.normalize in ['log', 'both']:
    log_message("  Adding log-transformed features...")
    df['genes_total_log'] = np.log1p(df['genes_total'])
    for col in go_cols:
        df[f"{col}_log"] = np.log1p(df[col])
    feature_cols += ['genes_total_log'] + [f"{col}_log" for col in go_cols]
    log_message(f"    Added {len(go_cols) + 1} log-transformed features")

# Total annotated domains
df['go_total'] = df[go_cols].sum(axis=1)
feature_cols.append('go_total')
log_message(f"  Added go_total (sum of all GO counts)")

# Remove features with zero variance
log_message("  Checking feature variance...")
feature_df = df[feature_cols].copy()
variances = feature_df.var()
zero_var_cols = variances[variances == 0].index.tolist()
if len(zero_var_cols) > 0:
    log_message(f"    Removing {len(zero_var_cols)} features with zero variance")
    feature_cols = [col for col in feature_cols if col not in zero_var_cols]

log_message(f"  Final feature count: {len(feature_cols)}")
log_message(f"  Feature matrix shape: {len(df)} × {len(feature_cols)}")

# Extract feature matrix and labels
X = df[feature_cols].copy()
y = df['environment'].copy()

# Check for any remaining missing values
if X.isna().sum().sum() > 0:
    log_message(f"  ✗ ERROR: {X.isna().sum().sum()} missing values in feature matrix")
    X = X.fillna(0)
    log_message("    Filled with 0")

# Check for infinite values
if np.isinf(X).sum().sum() > 0:
    log_message(f"  ✗ WARNING: {np.isinf(X).sum().sum()} infinite values in feature matrix")
    X = X.replace([np.inf, -np.inf], 0)

log_message("")

# ============================================================================
# 3. Data Splitting
# ============================================================================

log_message("3. Data Splitting")
log_message("-" * 80)

# Stratified split: 80% train, 10% val, 10% test
log_message("  Splitting data (stratified): 80% train / 10% val / 10% test")

# First split: 80% train, 20% temp
X_train, X_temp, y_train, y_temp = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)

# Second split: 50% val, 50% test from temp (10% each of total)
X_val, X_test, y_val, y_test = train_test_split(
    X_temp, y_temp, test_size=0.5, stratify=y_temp, random_state=42
)

log_message(f"  Training set: {len(X_train)} genomes")
log_message(f"  Validation set: {len(X_val)} genomes")
log_message(f"  Test set: {len(X_test)} genomes")

# Log class distribution in each split
for split_name, y_split in [("Training", y_train), ("Validation", y_val), ("Test", y_test)]:
    log_message(f"  {split_name} set class distribution:")
    split_counts = y_split.value_counts().sort_index()
    for env in valid_env_list:
        count = split_counts.get(env, 0)
        pct = 100 * count / len(y_split) if len(y_split) > 0 else 0
        log_message(f"    {env}: {count} ({pct:.1f}%)")

log_message("")

# ============================================================================
# 4. Model Training
# ============================================================================

log_message("4. Model Training")
log_message("-" * 80)

# Preprocessing: scale features for models that need it
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_val_scaled = scaler.transform(X_val)
X_test_scaled = scaler.transform(X_test)

# Save scaler
if HAS_JOBLIB:
    scaler_file = OUTPUT_DIR / f"{prefix}env_prediction_scaler.pkl"
    joblib.dump(scaler, scaler_file)
    log_message(f"  ✓ Saved scaler to {scaler_file.name}")
elif HAS_PICKLE:
    scaler_file = OUTPUT_DIR / f"{prefix}env_prediction_scaler.pkl"
    with open(scaler_file, 'wb') as f:
        pickle.dump(scaler, f)
    log_message(f"  ✓ Saved scaler to {scaler_file.name}")

models = {}
results = {}

# Baseline: Majority class
if args.model in ['all', 'baseline']:
    log_message("  Training baseline (majority class)...")
    baseline = DummyClassifier(strategy='most_frequent', random_state=42)
    baseline.fit(X_train, y_train)
    models['baseline'] = baseline
    log_message("    ✓ Baseline trained")

# Random Forest
if args.model in ['all', 'rf']:
    log_message("  Training Random Forest...")
    rf = RandomForestClassifier(
        n_estimators=100,
        max_depth=20,
        min_samples_split=5,
        min_samples_leaf=2,
        class_weight='balanced',
        random_state=42,
        n_jobs=-1
    )
    rf.fit(X_train, y_train)
    models['rf'] = rf
    log_message("    ✓ Random Forest trained")

# Gradient Boosting
if args.model in ['all', 'gb']:
    log_message("  Training Gradient Boosting...")
    gb = GradientBoostingClassifier(
        n_estimators=100,
        learning_rate=0.1,
        max_depth=5,
        random_state=42
    )
    gb.fit(X_train, y_train)
    models['gb'] = gb
    log_message("    ✓ Gradient Boosting trained")

# XGBoost (with optional GPU support)
if args.model in ['all', 'xgb'] and HAS_XGBOOST:
    log_message("  Training XGBoost...")
    # XGBoost needs numeric labels
    label_encoder_xgb = LabelEncoder()
    y_train_xgb = label_encoder_xgb.fit_transform(y_train)
    y_val_xgb = label_encoder_xgb.transform(y_val)
    y_test_xgb = label_encoder_xgb.transform(y_test)
    
    xgb_params = {
        'n_estimators': 100,
        'learning_rate': 0.1,
        'max_depth': 5,
        'random_state': 42,
        'tree_method': 'hist',
        'objective': 'multi:softprob',
        'num_class': len(valid_env_list),
        'eval_metric': 'mlogloss'
    }
    if args.use_gpu and GPU_AVAILABLE:
        xgb_params['tree_method'] = 'gpu_hist'
        xgb_params['gpu_id'] = 0
        log_message("    Using GPU acceleration")
    xgb_model = xgb.XGBClassifier(**xgb_params)
    xgb_model.fit(X_train, y_train_xgb)
    # Store encoder for predictions
    xgb_model.label_encoder_ = label_encoder_xgb
    models['xgb'] = xgb_model
    log_message("    ✓ XGBoost trained")
elif args.model in ['all', 'xgb'] and not HAS_XGBOOST:
    log_message("  ⚠ XGBoost requested but not available (install with: pip install xgboost)")

# Logistic Regression
if args.model in ['all', 'lr']:
    log_message("  Training Logistic Regression...")
    lr = LogisticRegression(
        max_iter=1000,
        C=1.0,
        class_weight='balanced',
        random_state=42,
        n_jobs=-1
    )
    lr.fit(X_train_scaled, y_train)
    models['lr'] = lr
    log_message("    ✓ Logistic Regression trained")

log_message("")

# ============================================================================
# 5. Evaluation Metrics
# ============================================================================

log_message("5. Evaluation Metrics")
log_message("-" * 80)

all_metrics = []
all_confusion_matrices = []
all_per_class_metrics = []
all_predictions = {}  # Store predictions and probabilities for visualization

# Label encoder for numeric encoding
label_encoder = LabelEncoder()
y_train_encoded = label_encoder.fit_transform(y_train)
y_val_encoded = label_encoder.transform(y_val)
y_test_encoded = label_encoder.transform(y_test)

for model_name, model in models.items():
    log_message(f"  Evaluating {model_name}...")
    
    # Select appropriate data (scaled for LR, raw for tree-based)
    if model_name == 'lr':
        X_train_eval = X_train_scaled
        X_val_eval = X_val_scaled
        X_test_eval = X_test_scaled
    else:
        X_train_eval = X_train
        X_val_eval = X_val
        X_test_eval = X_test
    
    # Handle XGBoost label encoding
    if model_name == 'xgb' and hasattr(model, 'label_encoder_'):
        y_train_enc = model.label_encoder_.transform(y_train)
        y_val_enc = model.label_encoder_.transform(y_val)
        y_test_enc = model.label_encoder_.transform(y_test)
        y_train_pred_enc = model.predict(X_train_eval)
        y_val_pred_enc = model.predict(X_val_eval)
        y_test_pred_enc = model.predict(X_test_eval)
        # Convert back to string labels
        y_train_pred = model.label_encoder_.inverse_transform(y_train_pred_enc)
        y_val_pred = model.label_encoder_.inverse_transform(y_val_pred_enc)
        y_test_pred = model.label_encoder_.inverse_transform(y_test_pred_enc)
    else:
        # Predictions
        y_train_pred = model.predict(X_train_eval)
        y_val_pred = model.predict(X_val_eval)
        y_test_pred = model.predict(X_test_eval)
    
    # Get prediction probabilities for ROC/AUC (if available)
    try:
        y_test_proba = model.predict_proba(X_test_eval)
        all_predictions[model_name] = {
            'y_test': y_test,
            'y_test_pred': y_test_pred,
            'y_test_proba': y_test_proba,
            'y_test_encoded': y_test_encoded
        }
    except:
        # Some models (like baseline) may not have predict_proba
        all_predictions[model_name] = {
            'y_test': y_test,
            'y_test_pred': y_test_pred,
            'y_test_proba': None,
            'y_test_encoded': y_test_encoded
        }
    
    # Metrics
    train_acc = accuracy_score(y_train, y_train_pred)
    val_acc = accuracy_score(y_val, y_val_pred)
    test_acc = accuracy_score(y_test, y_test_pred)
    
    train_bal_acc = balanced_accuracy_score(y_train, y_train_pred)
    val_bal_acc = balanced_accuracy_score(y_val, y_val_pred)
    test_bal_acc = balanced_accuracy_score(y_test, y_test_pred)
    
    # Confusion matrix
    cm_test = confusion_matrix(y_test, y_test_pred, labels=valid_env_list)
    
    # Per-class metrics
    precision, recall, f1, support = precision_recall_fscore_support(
        y_test, y_test_pred, labels=valid_env_list, zero_division=0
    )
    
    # Classification report
    report = classification_report(y_test, y_test_pred, labels=valid_env_list, output_dict=True)
    
    # Macro and weighted averages
    macro_f1 = report['macro avg']['f1-score']
    weighted_f1 = report['weighted avg']['f1-score']
    
    log_message(f"    Train accuracy: {train_acc:.4f} (balanced: {train_bal_acc:.4f})")
    log_message(f"    Val accuracy: {val_acc:.4f} (balanced: {val_bal_acc:.4f})")
    log_message(f"    Test accuracy: {test_acc:.4f} (balanced: {test_bal_acc:.4f})")
    log_message(f"    Macro F1: {macro_f1:.4f}, Weighted F1: {weighted_f1:.4f}")
    
    # Overfitting check
    overfit = train_acc - val_acc
    if overfit > 0.1:
        log_message(f"    ⚠ WARNING: Possible overfitting (train-val diff: {overfit:.4f})")
    
    # Store results
    results[model_name] = {
        'train_acc': train_acc,
        'val_acc': val_acc,
        'test_acc': test_acc,
        'train_bal_acc': train_bal_acc,
        'val_bal_acc': val_bal_acc,
        'test_bal_acc': test_bal_acc,
        'macro_f1': macro_f1,
        'weighted_f1': weighted_f1,
        'overfit': overfit
    }
    
    # Store metrics for output
    all_metrics.append({
        'model': model_name,
        'train_accuracy': train_acc,
        'val_accuracy': val_acc,
        'test_accuracy': test_acc,
        'train_balanced_accuracy': train_bal_acc,
        'val_balanced_accuracy': val_bal_acc,
        'test_balanced_accuracy': test_bal_acc,
        'macro_f1': macro_f1,
        'weighted_f1': weighted_f1,
        'overfit': overfit
    })
    
    # Confusion matrix
    cm_df = pd.DataFrame(cm_test, index=valid_env_list, columns=valid_env_list)
    cm_df.index.name = 'true_environment'
    cm_df.columns.name = 'predicted_environment'
    all_confusion_matrices.append((model_name, cm_df))
    
    # Per-class metrics
    per_class = []
    for i, env in enumerate(valid_env_list):
        per_class.append({
            'model': model_name,
            'environment': env,
            'precision': precision[i],
            'recall': recall[i],
            'f1_score': f1[i],
            'support': support[i]
        })
    all_per_class_metrics.extend(per_class)
    
    # Get prediction probabilities for ROC/AUC (if available)
    try:
        y_test_proba = model.predict_proba(X_test_eval)
        all_predictions[model_name] = {
            'y_test': y_test,
            'y_test_pred': y_test_pred,
            'y_test_proba': y_test_proba,
            'y_test_encoded': y_test_encoded
        }
    except:
        # Some models (like baseline) may not have predict_proba
        all_predictions[model_name] = {
            'y_test': y_test,
            'y_test_pred': y_test_pred,
            'y_test_proba': None,
            'y_test_encoded': y_test_encoded
        }

log_message("")

# ============================================================================
# Load GO Labels for Feature Annotation
# ============================================================================

log_message("Loading GO term labels for feature annotation...")
go_labels_map = {}
try:
    if GO_LABELS_FILE.exists():
        go_labels_df = pd.read_csv(GO_LABELS_FILE, sep='\t')
        go_labels_df['category'] = go_labels_df['category'].astype(str).str.zfill(7)
        # Create mapping: category -> full_label
        for _, row in go_labels_df.iterrows():
            cat = str(row['category']).zfill(7)
            # Use full_label if available, otherwise short_label
            label = row['full_label'] if pd.notna(row['full_label']) else row['short_label']
            go_labels_map[cat] = label
        log_message(f"  ✓ Loaded {len(go_labels_map)} GO term labels")
    else:
        log_message(f"  ⚠ GO labels file not found: {GO_LABELS_FILE}")
        log_message("    Will use numeric IDs for features")
except Exception as e:
    log_message(f"  ⚠ Error loading GO labels: {e}")
    log_message("    Will use numeric IDs for features")

def get_feature_label(feature_name):
    """Get human-readable label for a feature."""
    # Handle special features
    if feature_name == 'genes_total':
        return 'Total genes (genome size)'
    elif feature_name == 'go_total':
        return 'Total annotated GO domains'
    elif feature_name.startswith('genes_total_'):
        return f'{feature_name.replace("genes_total_", "")} (genome size)'
    elif feature_name.endswith('_per_gene'):
        base = feature_name.replace('_per_gene', '')
        if base in go_labels_map:
            return f'{go_labels_map[base]} (per gene)'
        return f'{base} (per gene)'
    elif feature_name.endswith('_log'):
        base = feature_name.replace('_log', '')
        if base in go_labels_map:
            return f'{go_labels_map[base]} (log)'
        return f'{base} (log)'
    else:
        # Check if it's a GO term (7-digit format)
        if len(feature_name) == 7 and feature_name.isdigit():
            if feature_name in go_labels_map:
                return go_labels_map[feature_name]
        return feature_name

log_message("")

# ============================================================================
# 6. Visualization: Performance Metrics and Figures
# ============================================================================

if args.plot and HAS_MATPLOTLIB:
    log_message("6. Generating Visualizations")
    log_message("-" * 80)
    
    # Set style
    sns.set_style("whitegrid")
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['font.size'] = 10
    
    def save_figure(fig, filename):
        """Save figure in both PNG and PDF formats."""
        filename_with_prefix = f"{prefix}{filename}" if prefix else filename
        png_path = OUTPUT_DIR / f"{filename_with_prefix}.png"
        pdf_path = OUTPUT_DIR / f"{filename_with_prefix}.pdf"
        fig.savefig(png_path, bbox_inches='tight', dpi=300)
        fig.savefig(pdf_path, bbox_inches='tight')
        log_message(f"  ✓ Saved {filename_with_prefix}.png and .pdf")
        plt.close(fig)
    
    # Figure 1: Model Performance Comparison
    log_message("  Generating Figure 1: Model Performance Comparison...")
    try:
        metrics_df = pd.DataFrame(all_metrics)
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Test accuracy comparison
        ax = axes[0, 0]
        metrics_df_sorted = metrics_df.sort_values('test_accuracy', ascending=False)
        bars = ax.barh(metrics_df_sorted['model'], metrics_df_sorted['test_accuracy'], color='steelblue')
        ax.set_xlabel('Test Accuracy', fontsize=12)
        ax.set_title('Test Accuracy by Model', fontsize=14, weight='bold')
        ax.set_xlim([0, 1])
        for i, (idx, row) in enumerate(metrics_df_sorted.iterrows()):
            ax.text(row['test_accuracy'] + 0.01, i, f"{row['test_accuracy']:.3f}", 
                   va='center', fontsize=9)
        ax.grid(True, alpha=0.3, axis='x')
        
        # Balanced accuracy comparison
        ax = axes[0, 1]
        metrics_df_sorted = metrics_df.sort_values('test_balanced_accuracy', ascending=False)
        bars = ax.barh(metrics_df_sorted['model'], metrics_df_sorted['test_balanced_accuracy'], color='darkgreen')
        ax.set_xlabel('Test Balanced Accuracy', fontsize=12)
        ax.set_title('Test Balanced Accuracy by Model', fontsize=14, weight='bold')
        ax.set_xlim([0, 1])
        for i, (idx, row) in enumerate(metrics_df_sorted.iterrows()):
            ax.text(row['test_balanced_accuracy'] + 0.01, i, f"{row['test_balanced_accuracy']:.3f}", 
                   va='center', fontsize=9)
        ax.grid(True, alpha=0.3, axis='x')
        
        # Macro F1 comparison
        ax = axes[1, 0]
        metrics_df_sorted = metrics_df.sort_values('macro_f1', ascending=False)
        bars = ax.barh(metrics_df_sorted['model'], metrics_df_sorted['macro_f1'], color='coral')
        ax.set_xlabel('Macro F1-Score', fontsize=12)
        ax.set_title('Macro F1-Score by Model', fontsize=14, weight='bold')
        ax.set_xlim([0, 1])
        for i, (idx, row) in enumerate(metrics_df_sorted.iterrows()):
            ax.text(row['macro_f1'] + 0.01, i, f"{row['macro_f1']:.3f}", 
                   va='center', fontsize=9)
        ax.grid(True, alpha=0.3, axis='x')
        
        # Train vs Test accuracy (overfitting check)
        ax = axes[1, 1]
        x_pos = np.arange(len(metrics_df))
        width = 0.35
        ax.bar(x_pos - width/2, metrics_df['train_accuracy'], width, label='Train', color='lightblue')
        ax.bar(x_pos + width/2, metrics_df['test_accuracy'], width, label='Test', color='steelblue')
        ax.set_xlabel('Model', fontsize=12)
        ax.set_ylabel('Accuracy', fontsize=12)
        ax.set_title('Train vs Test Accuracy (Overfitting Check)', fontsize=14, weight='bold')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(metrics_df['model'], rotation=45, ha='right')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        save_figure(fig, "fig01_model_performance_comparison")
    except Exception as e:
        log_message(f"    ✗ ERROR generating Figure 1: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Figure 2: Confusion Matrices
    log_message("  Generating Figure 2: Confusion Matrices...")
    try:
        n_models = len(all_confusion_matrices)
        n_cols = min(3, n_models)
        n_rows = (n_models + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
        if n_models == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes if isinstance(axes, list) else [axes]
        else:
            axes = axes.flatten()
        
        for idx, (model_name, cm_df) in enumerate(all_confusion_matrices):
            ax = axes[idx] if n_models > 1 else axes[0]
            
            # Normalize confusion matrix to percentages
            cm_normalized = cm_df.values.astype('float') / cm_df.values.sum(axis=1)[:, np.newaxis]
            cm_normalized = np.nan_to_num(cm_normalized)
            
            sns.heatmap(cm_normalized, annot=True, fmt='.2f', cmap='Blues', 
                       xticklabels=cm_df.columns, yticklabels=cm_df.index,
                       ax=ax, cbar_kws={'label': 'Proportion'})
            ax.set_xlabel('Predicted Environment', fontsize=10)
            ax.set_ylabel('True Environment', fontsize=10)
            ax.set_title(f'{model_name.upper()} Confusion Matrix', fontsize=12, weight='bold')
            ax.tick_params(axis='x', rotation=45)
            ax.tick_params(axis='y', rotation=0)
        
        # Hide unused subplots
        for idx in range(n_models, len(axes)):
            axes[idx].set_visible(False)
        
        plt.tight_layout()
        save_figure(fig, "fig02_confusion_matrices")
    except Exception as e:
        log_message(f"    ✗ ERROR generating Figure 2: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Figure 3: ROC Curves (AUC) for models with probabilities
    log_message("  Generating Figure 3: ROC Curves (AUC)...")
    try:
        models_with_proba = {k: v for k, v in all_predictions.items() if v['y_test_proba'] is not None}
        if len(models_with_proba) > 0:
            n_models = len(models_with_proba)
            n_cols = min(3, n_models)
            n_rows = (n_models + n_cols - 1) // n_cols
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
            if n_models == 1:
                axes = [axes]
            elif n_rows == 1:
                axes = axes if isinstance(axes, list) else [axes]
            else:
                axes = axes.flatten()
            
            for idx, (model_name, pred_data) in enumerate(models_with_proba.items()):
                ax = axes[idx] if n_models > 1 else axes[0]
                
                y_test_enc = pred_data['y_test_encoded']
                y_proba = pred_data['y_test_proba']
                n_classes = len(valid_env_list)
                
                # Compute ROC curve for each class (one-vs-rest)
                auc_scores = []
                for i, env in enumerate(valid_env_list):
                    y_binary = (y_test_enc == i).astype(int)
                    if y_binary.sum() > 0:  # Only if class is present
                        fpr, tpr, _ = roc_curve(y_binary, y_proba[:, i])
                        roc_auc = auc(fpr, tpr)
                        auc_scores.append(roc_auc)
                        ax.plot(fpr, tpr, lw=2, alpha=0.7, 
                               label=f'{env} (AUC={roc_auc:.3f})')
                
                ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5, label='Random (AUC=0.5)')
                ax.set_xlabel('False Positive Rate', fontsize=10)
                ax.set_ylabel('True Positive Rate', fontsize=10)
                mean_auc = np.mean(auc_scores) if auc_scores else 0
                ax.set_title(f'{model_name.upper()} ROC Curves\n(Mean AUC={mean_auc:.3f})', 
                           fontsize=12, weight='bold')
                ax.legend(loc='lower right', fontsize=8)
                ax.grid(True, alpha=0.3)
            
            # Hide unused subplots
            for idx in range(n_models, len(axes)):
                axes[idx].set_visible(False)
            
            plt.tight_layout()
            save_figure(fig, "fig03_roc_curves")
        else:
            log_message("    ⚠ No models with probability predictions for ROC curves")
    except Exception as e:
        log_message(f"    ✗ ERROR generating Figure 3: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Figure 4: Per-Class Performance Metrics
    log_message("  Generating Figure 4: Per-Class Performance Metrics...")
    try:
        per_class_df = pd.DataFrame(all_per_class_metrics)
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # Precision by environment and model
        ax = axes[0, 0]
        pivot_precision = per_class_df.pivot(index='environment', columns='model', values='precision')
        sns.heatmap(pivot_precision, annot=True, fmt='.3f', cmap='YlOrRd', ax=ax, cbar_kws={'label': 'Precision'})
        ax.set_title('Precision by Environment and Model', fontsize=12, weight='bold')
        ax.set_xlabel('Model', fontsize=10)
        ax.set_ylabel('Environment', fontsize=10)
        ax.tick_params(axis='x', rotation=45)
        ax.tick_params(axis='y', rotation=0)
        
        # Recall by environment and model
        ax = axes[0, 1]
        pivot_recall = per_class_df.pivot(index='environment', columns='model', values='recall')
        sns.heatmap(pivot_recall, annot=True, fmt='.3f', cmap='YlGn', ax=ax, cbar_kws={'label': 'Recall'})
        ax.set_title('Recall by Environment and Model', fontsize=12, weight='bold')
        ax.set_xlabel('Model', fontsize=10)
        ax.set_ylabel('Environment', fontsize=10)
        ax.tick_params(axis='x', rotation=45)
        ax.tick_params(axis='y', rotation=0)
        
        # F1-Score by environment and model
        ax = axes[1, 0]
        pivot_f1 = per_class_df.pivot(index='environment', columns='model', values='f1_score')
        sns.heatmap(pivot_f1, annot=True, fmt='.3f', cmap='Blues', ax=ax, cbar_kws={'label': 'F1-Score'})
        ax.set_title('F1-Score by Environment and Model', fontsize=12, weight='bold')
        ax.set_xlabel('Model', fontsize=10)
        ax.set_ylabel('Environment', fontsize=10)
        ax.tick_params(axis='x', rotation=45)
        ax.tick_params(axis='y', rotation=0)
        
        # Support (sample size) by environment
        ax = axes[1, 1]
        support_by_env = per_class_df.groupby('environment')['support'].first().sort_values(ascending=False)
        bars = ax.barh(range(len(support_by_env)), support_by_env.values, color='steelblue')
        ax.set_yticks(range(len(support_by_env)))
        ax.set_yticklabels(support_by_env.index)
        ax.set_xlabel('Number of Test Samples', fontsize=10)
        ax.set_title('Test Set Sample Size by Environment', fontsize=12, weight='bold')
        for i, (env, count) in enumerate(support_by_env.items()):
            ax.text(count + max(support_by_env.values)*0.01, i, f"{count}", 
                   va='center', fontsize=9)
        ax.grid(True, alpha=0.3, axis='x')
        
        plt.tight_layout()
        save_figure(fig, "fig04_per_class_metrics")
    except Exception as e:
        log_message(f"    ✗ ERROR generating Figure 4: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    log_message("")

# ============================================================================
# 7. Feature Importance Analysis
# ============================================================================

log_message("6. Feature Importance Analysis")
log_message("-" * 80)

all_feature_importances = []

for model_name, model in models.items():
    if model_name == 'baseline':
        continue
    
    log_message(f"  Extracting feature importances from {model_name}...")
    
    if hasattr(model, 'feature_importances_'):
        # Tree-based models
        importances = model.feature_importances_
        feature_imp_df = pd.DataFrame({
            'model': model_name,
            'feature': feature_cols,
            'importance': importances
        }).sort_values('importance', ascending=False)
        
    elif hasattr(model, 'coef_'):
        # Logistic Regression
        # Average absolute coefficients across all classes
        coef_abs = np.abs(model.coef_).mean(axis=0)
        feature_imp_df = pd.DataFrame({
            'model': model_name,
            'feature': feature_cols,
            'importance': coef_abs
        }).sort_values('importance', ascending=False)
    else:
        log_message(f"    ⚠ No feature importance available for {model_name}")
        continue
    
    all_feature_importances.append(feature_imp_df)
    
    # Log top 10 features
    top_10 = feature_imp_df.head(10)
    log_message(f"    Top 10 features:")
    for idx, row in top_10.iterrows():
        log_message(f"      {row['feature']}: {row['importance']:.6f}")

log_message("")

# ============================================================================
# 8. Visualization: Feature Importance Plots
# ============================================================================

if args.plot and HAS_MATPLOTLIB and len(all_feature_importances) > 0:
    log_message("8. Generating Feature Importance Visualizations")
    log_message("-" * 80)
    
    # Figure 5: Feature Importance Comparison
    log_message("  Generating Figure 5: Top Feature Importances...")
    try:
        n_models = len(all_feature_importances)
        n_cols = min(3, n_models)
        n_rows = (n_models + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(10*n_cols, 10*n_rows))
        if n_models == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes if isinstance(axes, list) else [axes]
        else:
            axes = axes.flatten()
        
        for idx, feature_imp_df in enumerate(all_feature_importances):
            ax = axes[idx] if n_models > 1 else axes[0]
            model_name = feature_imp_df['model'].iloc[0]
            
            # Get top 20 features
            top_features = feature_imp_df.head(20).copy()
            
            # Add labels
            top_features['label'] = top_features['feature'].apply(get_feature_label)
            
            y_pos = np.arange(len(top_features))
            bars = ax.barh(y_pos, top_features['importance'].values, color='steelblue')
            ax.set_yticks(y_pos)
            
            # Use labels for y-axis, truncate if too long
            labels = top_features['label'].values
            max_label_len = 70  # Truncate long labels
            labels_truncated = [label[:max_label_len] + '...' if len(label) > max_label_len else label 
                               for label in labels]
            ax.set_yticklabels(labels_truncated, fontsize=9)
            ax.set_xlabel('Feature Importance', fontsize=11)
            ax.set_title(f'{model_name.upper()} Top 20 Features', fontsize=13, weight='bold')
            ax.invert_yaxis()
            ax.grid(True, alpha=0.3, axis='x')
            
            # Add value labels on bars
            max_imp = max(top_features['importance'].values)
            for i, (pos, imp) in enumerate(zip(y_pos, top_features['importance'].values)):
                ax.text(imp + max_imp * 0.01, pos, 
                       f'{imp:.4f}', va='center', fontsize=7)
        
        # Hide unused subplots
        for idx in range(n_models, len(axes)):
            axes[idx].set_visible(False)
        
        plt.tight_layout()
        save_figure(fig, "fig05_feature_importance")
        log_message("    ✓ Figure 5 saved with annotated GO term labels")
    except Exception as e:
        log_message(f"    ✗ ERROR generating Figure 5: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    log_message("")

# ============================================================================
# 9. Save Outputs
# ============================================================================

log_message("9. Saving Outputs")
log_message("-" * 80)

# Save metrics
metrics_df = pd.DataFrame(all_metrics)
metrics_file = OUTPUT_DIR / f"{prefix}env_prediction_metrics.tsv"
metrics_df.to_csv(metrics_file, sep='\t', index=False)
log_message(f"  ✓ Saved metrics to {metrics_file.name}")

# Save confusion matrices
for model_name, cm_df in all_confusion_matrices:
    cm_file = OUTPUT_DIR / f"{prefix}env_prediction_confusion_matrix_{model_name}.tsv"
    cm_df.to_csv(cm_file, sep='\t')
    log_message(f"  ✓ Saved confusion matrix ({model_name}) to {cm_file.name}")

# Save per-class metrics
per_class_df = pd.DataFrame(all_per_class_metrics)
per_class_file = OUTPUT_DIR / f"{prefix}env_prediction_per_class_metrics.tsv"
per_class_df.to_csv(per_class_file, sep='\t', index=False)
log_message(f"  ✓ Saved per-class metrics to {per_class_file.name}")

# Save feature importances (with labels if available)
if len(all_feature_importances) > 0:
    feature_imp_df = pd.concat(all_feature_importances, ignore_index=True)
    
    # Add labels if GO labels were loaded
    if len(go_labels_map) > 0:
        feature_imp_df['feature_label'] = feature_imp_df['feature'].apply(get_feature_label)
        # Reorder columns: model, feature, feature_label, importance
        feature_imp_df = feature_imp_df[['model', 'feature', 'feature_label', 'importance']]
        log_message(f"  ✓ Added GO term labels to feature importances")
    
    feature_imp_file = OUTPUT_DIR / f"{prefix}env_prediction_feature_importances.tsv"
    feature_imp_df.to_csv(feature_imp_file, sep='\t', index=False)
    log_message(f"  ✓ Saved feature importances to {feature_imp_file.name}")

# Save models
for model_name, model in models.items():
    model_file = OUTPUT_DIR / f"{prefix}env_prediction_model_{model_name}.pkl"
    if HAS_JOBLIB:
        joblib.dump(model, model_file)
    elif HAS_PICKLE:
        with open(model_file, 'wb') as f:
            pickle.dump(model, f)
    log_message(f"  ✓ Saved model ({model_name}) to {model_file.name}")

log_message("")

# ============================================================================
# 10. Robustness Checks
# ============================================================================

log_message("8. Robustness Checks")
log_message("-" * 80)

# Taxonomic distribution
if 'organism_taxId' in df.columns:
    log_message("  Taxonomic distribution per environment:")
    tax_dist = df.groupby('environment')['organism_taxId'].nunique().sort_values(ascending=False)
    for env, n_taxa in tax_dist.items():
        log_message(f"    {env}: {n_taxa} unique taxa")
    
    # Check for potential phylogenetic confounding
    total_taxa = df['organism_taxId'].nunique()
    log_message(f"  Total unique taxa: {total_taxa}")
    if total_taxa < len(df) * 0.1:
        log_message("    ⚠ WARNING: Low taxonomic diversity (potential phylogenetic confounding)")

# GO annotation completeness
log_message("  GO annotation completeness:")
for env in valid_env_list:
    env_df = df[df['environment'] == env]
    mean_go_total = env_df['go_total'].mean()
    log_message(f"    {env}: {mean_go_total:.1f} mean total GO counts per genome")

# Feature statistics
log_message("  Feature matrix statistics:")
log_message(f"    Mean: {X.mean().mean():.2f}")
log_message(f"    Std: {X.std().mean():.2f}")
log_message(f"    Min: {X.min().min():.0f}")
log_message(f"    Max: {X.max().max():.0f}")

log_message("")

# ============================================================================
# Final Summary
# ============================================================================

log_message("=" * 80)
log_message("Summary")
log_message("=" * 80)
log_message(f"  Total genomes: {len(df)}")
log_message(f"  Environments: {len(valid_env_list)}")
log_message(f"  Features: {len(feature_cols)}")
log_message(f"  Models trained: {len(models)}")

log_message("  Best test accuracy by model:")
for model_name in models.keys():
    if model_name in results:
        test_acc = results[model_name]['test_acc']
        test_bal_acc = results[model_name]['test_bal_acc']
        log_message(f"    {model_name}: {test_acc:.4f} (balanced: {test_bal_acc:.4f})")

log_message("")

# Write QC log
with open(QC_LOG_FILE, 'w') as f:
    f.write('\n'.join(qc_log))

log_message(f"✓ QC log written to {QC_LOG_FILE.name}")
log_message("")
log_message("Script 05 completed successfully!")

