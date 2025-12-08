# Plan: KEGG GO Analyses Python Script

## Current Understanding

### Input Files Structure
- **Location**: `results/3_GO_analyses/GO_lists_with_counts/GCF_*.txt`
- **Format**: Tab-separated, one GO term per line
  - Column 1: GO term ID (7-digit format, e.g., `0004427`)
  - Column 2: Count (integer, number of times GO term appears in genome)
- **Example**:
  ```
  0006260	10
  0003677	58
  0004427	5
  ```
- **Total genomes**: 3,088 files
- **GO terms per genome**: Variable (e.g., ~900 for GCF_000006985.1)

### Mapping Structure
- **GO → Reaction**: From `data/kegg/kegg_reaction2go.txt`
  - Format: `KEGG_REACTION:R00004 > GO:description ; GO:0004427`
  - Need to extract: GO term ID (7-digit) → Reaction ID (R00004)
- **Reaction → KO**: From KEGG REST API (`/link/ko/reaction_id`)
- **Reaction → Pathway**: From KEGG REST API (`/link/pathway/reaction_id`)

### Output Requirements
**ALL outputs must be in**: `results/3.5_KEGG_n_reaction_analyses/`

**Required outputs**:
1. `ubiquitous_counts_table_kegg_reaction.txt` - Table with ubiquitous reactions (≥95% prevalence)
2. `REACTION_lists_with_counts_/R*.txt` - Individual reaction files (genome → count)
3. `KO_lists_with_counts_/K*.txt` - Individual KO files (genome → count)
4. `intermediate_mappings/` - Mapping files for reproducibility
   - `go_to_reactions.tsv` - GO term ID → Reaction IDs
   - `reaction_to_gos.tsv` - Reaction ID → GO term IDs
   - `reaction_to_kos.tsv` - Reaction ID → KO IDs
   - `reaction_to_pathways.tsv` - Reaction ID → Pathway IDs

## Implementation Plan

### Phase 1: Safety & Setup
1. **Safety Functions**:
   - `validate_output_path(file_path)`: Ensures path is within RESULTS_DIR
   - `safe_open_write(file_path, mode)`: Wrapper for file writing with validation
   - All file writes must use `safe_open_write()`

2. **Directory Structure**:
   - Create all output directories upfront
   - Verify input files exist

3. **Logging**:
   - Detailed logging to file: `kegg_go_analyses.log` (in RESULTS_DIR)
   - Log all major steps, file I/O, API calls, progress
   - Include timestamps and memory usage where relevant

### Phase 2: Data Loading
1. **Load GO → Reaction Mapping**:
   - Parse `data/kegg/kegg_reaction2go.txt`
   - Build: `go_to_reactions[go_term_id] = set([reaction_ids])`
   - Save to `intermediate_mappings/go_to_reactions.tsv` (using GO term IDs, not descriptions)

2. **Load GO Counts from Individual Genome Files**:
   - Iterate through `GO_lists_with_counts/GCF_*.txt`
   - For each genome:
     - Load GO term → count dictionary
     - Store in: `go_counts_per_genome[genome_id][go_term] = count`
   - **Note**: Use ALL GO terms, not just ubiquitous ones

### Phase 3: KEGG API Queries
1. **Collect All Reaction IDs**:
   - From GO → Reaction mapping, get all unique reaction IDs
   - Batch into groups of 100 for API calls

2. **Query Reaction → KO**:
   - Use batch API: `/link/ko/rn:R00004+rn:R00005+...`
   - Parse results: `reaction_to_kos[reaction_id] = set([ko_ids])`
   - Save to `intermediate_mappings/reaction_to_kos.tsv`
   - Cache API responses to `kegg_api_cache/` for reproducibility

3. **Query Reaction → Pathway**:
   - Use batch API: `/link/pathway/rn:R00004+rn:R00005+...`
   - Parse results: `reaction_to_pathways[reaction_id] = set([pathway_ids])`
   - Save to `intermediate_mappings/reaction_to_pathways.tsv`

### Phase 4: Build Count Matrices
1. **Reaction-Level Counts**:
   - For each genome:
     - For each GO term with count:
       - Get reactions mapped to that GO term
       - Add GO count to each reaction
   - Result: `reaction_counts_per_genome[genome_id][reaction_id] = total_count`

2. **KO-Level Counts**:
   - For each genome:
     - For each reaction with count:
       - Get KOs mapped to that reaction
       - Add reaction count to each KO
   - Result: `ko_counts_per_genome[genome_id][ko_id] = total_count`

3. **Pathway-Level Counts** (optional):
   - Similar aggregation from reactions to pathways

### Phase 5: Write Output Files
1. **Identify Ubiquitous Reactions**:
   - Calculate prevalence: reactions present in ≥95% of genomes
   - Filter to ubiquitous set

2. **Write Ubiquitous Reaction Table**:
   - Format: `Genome\tR00130\tR00253\t...`
   - Only include ubiquitous reactions

3. **Write Individual Reaction Files**:
   - For each reaction: `REACTION_lists_with_counts_/R00004.txt`
   - Format: `genome_id\tcount` (only non-zero counts)

4. **Write Individual KO Files**:
   - For each KO: `KO_lists_with_counts_/K00001.txt`
   - Format: `genome_id\tcount` (only non-zero counts)

### Phase 6: Testing Strategy
1. **Subset Test** (first 10 genomes):
   - Run full pipeline on subset
   - Verify:
     - All outputs in correct directory
     - Counts make sense (non-negative, reasonable values)
     - Mappings are correct
     - File formats are correct
   - Check logs for errors/warnings

2. **Full Run**:
   - Once subset test passes, run on all 3,088 genomes
   - Monitor progress and memory usage
   - Log completion status

## Key Differences from Notebook

1. **Input Source**: Use `GO_lists_with_counts/` files instead of `ubiquitous_counts_table.txt`
   - This gives us ALL GO terms, not just ubiquitous ones
   - More comprehensive mapping

2. **Safety**: All file writes validated to prevent corruption of upstream data

3. **Logging**: Comprehensive logging for debugging and reproducibility

4. **Error Handling**: Robust error handling for API calls, file I/O, edge cases

5. **Progress Tracking**: Progress bars/counters for long-running operations

## Expected Outputs Summary

- `ubiquitous_counts_table_kegg_reaction.txt`: ~3,088 genomes × ~100-200 ubiquitous reactions
- `REACTION_lists_with_counts_/`: ~1,000-2,000 reaction files
- `KO_lists_with_counts_/`: ~500-1,000 KO files
- `intermediate_mappings/`: 4 TSV files with mappings
- `kegg_api_cache/`: Cached API responses (optional, for reproducibility)
- `kegg_go_analyses.log`: Detailed execution log

## Validation Checklist

- [ ] All file paths validated (within RESULTS_DIR)
- [ ] Input files exist and are readable
- [ ] GO → Reaction mapping parsed correctly
- [ ] KEGG API queries successful (with retry logic)
- [ ] Count matrices built correctly (sums make sense)
- [ ] Output files written to correct locations
- [ ] File formats match expected structure
- [ ] No files written outside RESULTS_DIR
- [ ] Logging captures all major steps

