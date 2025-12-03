#!/usr/bin/env python3
"""
Safe XLSX to CSV conversion script.
Handles multiple sheets, encoding, and preserves data integrity.
"""

import pandas as pd
import sys
import os
from pathlib import Path

def convert_xlsx_to_csv(xlsx_path, output_dir=None, sheet_name=None):
    """
    Safely convert XLSX file to CSV format.
    
    Parameters:
    -----------
    xlsx_path : str
        Path to the input XLSX file
    output_dir : str, optional
        Directory to save CSV files. If None, saves in same directory as XLSX.
    sheet_name : str or int, optional
        Specific sheet to convert. If None, converts all sheets.
    """
    xlsx_path = Path(xlsx_path)
    
    if not xlsx_path.exists():
        raise FileNotFoundError(f"File not found: {xlsx_path}")
    
    if output_dir is None:
        output_dir = xlsx_path.parent
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read Excel file
    try:
        xl_file = pd.ExcelFile(xlsx_path, engine='openpyxl')
        sheet_names = xl_file.sheet_names
        print(f"Found {len(sheet_names)} sheet(s): {sheet_names}")
    except Exception as e:
        raise RuntimeError(f"Error reading Excel file: {e}")
    
    # Determine which sheets to convert
    if sheet_name is not None:
        if sheet_name not in sheet_names and sheet_name not in range(len(sheet_names)):
            raise ValueError(f"Sheet '{sheet_name}' not found. Available sheets: {sheet_names}")
        sheets_to_convert = [sheet_name]
    else:
        sheets_to_convert = sheet_names
    
    converted_files = []
    
    # Convert each sheet
    for sheet in sheets_to_convert:
        try:
            print(f"\nReading sheet: {sheet}")
            df = pd.read_excel(xlsx_path, sheet_name=sheet, engine='openpyxl')
            print(f"  Shape: {df.shape[0]} rows × {df.shape[1]} columns")
            
            # Generate output filename
            if len(sheets_to_convert) == 1:
                # Single sheet: use base filename
                csv_filename = xlsx_path.stem + '.csv'
            else:
                # Multiple sheets: append sheet name
                safe_sheet_name = str(sheet).replace('/', '_').replace('\\', '_')
                csv_filename = f"{xlsx_path.stem}_{safe_sheet_name}.csv"
            
            csv_path = output_dir / csv_filename
            
            # Write to CSV with safe settings
            df.to_csv(
                csv_path,
                index=False,  # Don't write row indices
                encoding='utf-8',  # UTF-8 encoding for special characters
                quoting=1,  # QUOTE_ALL for safety with special characters
                na_rep='',  # Empty string for missing values
                float_format='%.10g'  # Preserve precision for floats
            )
            
            print(f"  Saved: {csv_path}")
            converted_files.append(csv_path)
            
        except Exception as e:
            print(f"  ERROR converting sheet '{sheet}': {e}", file=sys.stderr)
            continue
    
    print(f"\n✓ Conversion complete! Converted {len(converted_files)} sheet(s).")
    return converted_files

if __name__ == '__main__':
    # Default: convert the gold metadata file
    default_xlsx = Path(__file__).parent.parent / 'data' / '0_20251106_gold_metadata.xlsx'
    
    if len(sys.argv) > 1:
        xlsx_path = sys.argv[1]
    else:
        xlsx_path = default_xlsx
    
    # Optional: specify sheet name as second argument
    sheet_name = sys.argv[2] if len(sys.argv) > 2 else None
    
    try:
        convert_xlsx_to_csv(xlsx_path, sheet_name=sheet_name)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

