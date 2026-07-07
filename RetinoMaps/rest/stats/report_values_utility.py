#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
report_values_utility.py
------------------------------------------------------------------------------------------
Goal:
Grab results and paste them into the results section of the manuscript
"""

import os
import pandas as pd
from pathlib import Path
from typing import Dict, Tuple


def load_correlation_csv(csv_path: str) -> pd.DataFrame:
    """Load correlation CSV and return as DataFrame."""
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV not found: {csv_path}")
    return pd.read_csv(csv_path, index_col=0)


def extract_correlation_data(
    group_folder: str,
    corr_type: str = "full_corr",
    run_variant: str = "concat_clean"
) -> Dict[str, Dict[str, pd.DataFrame]]:
    """
    Load median, p25, p75 CSVs for both hemispheres.
    
    Args:
        group_folder: Path to group folder
        corr_type: 'full_corr' or 'partial_corr' (converted to hyphenated form for filename)
        run_variant: 'concat_clean', 'concat', 'run-01', 'run-02'
    
    Returns:
        Nested dict: {hemi: {stat: dataframe}}
    """
    corr_folder = Path(group_folder) / corr_type / "by_hemi" / "task-constrained"
    
    # Convert underscore to hyphen for filename construction
    corr_type_hyphenated = corr_type.replace("_", "-")
    
    data = {}
    for hemi in ["lh", "rh"]:
        data[hemi] = {}
        for stat in ["median", "p25", "p75"]:
            filename = (
                f"seed-task_by_macror-task_{corr_type_hyphenated}_r_{stat}_{run_variant}_{hemi}_legacy.csv"
            )
            filepath = corr_folder / filename
            data[hemi][stat] = load_correlation_csv(str(filepath))
    
    return data


def _format_paragraph(
    data: Dict[str, Dict[str, pd.DataFrame]],
    corr_type: str
) -> str:
    """
    Generate manuscript paragraph from loaded correlation data.
    
    Structure:
      - For each target (mPCS, sPCS, iPCS, sIPS, iIPS)
      - For each hemisphere (left = lh, right = rh)
      - For each seed, format: "seed_name seed: r = X.XX, [X.XX; X.XX]"
    
    Args:
        data: Nested dict from extract_correlation_data()
        corr_type: 'full_corr' or 'partial_corr'
    
    Returns:
        Formatted paragraph string ready to paste into manuscript
    """
    targets = ["mPCS", "sPCS", "iPCS", "sIPS", "iIPS"]
    seeds = ["mPCS", "sPCS", "iPCS", "sIPS", "iIPS"]
    
    target_paragraphs = []
    
    for target in targets:
        hemi_texts = []
        
        # Iterate hemispheres: lh → left, rh → right
        for hemi_label, hemi_code in [("left hemisphere", "lh"), ("right hemisphere", "rh")]:
            median_df = data[hemi_code]["median"]
            p25_df = data[hemi_code]["p25"]
            p75_df = data[hemi_code]["p75"]
            
            seed_texts = []
            for seed in seeds:
                # All data in each hemi CSV uses _ipsi suffix
                # (ipsi for lh = left side; ipsi for rh = right side)
                col_name = f"{target}_ipsi"
                
                median_val = median_df.loc[seed, col_name]
                p25_val = p25_df.loc[seed, col_name]
                p75_val = p75_df.loc[seed, col_name]
                
                seed_str = f"{seed} seed: r = {median_val:.2f}, [{p25_val:.2f}; {p75_val:.2f}]"
                seed_texts.append(seed_str)
            
            hemi_str = f"{hemi_label}: " + ", ".join(seed_texts)
            hemi_texts.append(hemi_str)
        
        target_text = f"{target} target ({'; '.join(hemi_texts)})"
        target_paragraphs.append(target_text)
    
    return ", ".join(target_paragraphs)


def format_manuscript_results(
    group_folder: str,
    run_variant: str = "concat_clean",
    output_file: str = None
) -> Tuple[str, str]:
    """
    Generate both full-correlation and partial-correlation paragraphs.
    
    Args:
        group_folder: Path to group folder
        run_variant: Run variant to use (default: concat_clean)
        output_file: Optional file path to save results
    
    Returns:
        Tuple of (full_corr_text, partial_corr_text)
    """
    print("Loading full correlation data...")
    full_corr_data = extract_correlation_data(group_folder, "full_corr", run_variant)
    
    print("Loading partial correlation data...")
    partial_corr_data = extract_correlation_data(group_folder, "partial_corr", run_variant)
    
    print("Formatting paragraphs...")
    full_text = _format_paragraph(full_corr_data, "full_corr")
    partial_text = _format_paragraph(partial_corr_data, "partial_corr")
    
    if output_file:
        with open(output_file, "w") as f:
            f.write("=" * 80 + "\n")
            f.write("FULL CORRELATION RESULTS\n")
            f.write("=" * 80 + "\n\n")
            f.write(full_text + "\n\n")
            f.write("=" * 80 + "\n")
            f.write("PARTIAL CORRELATION RESULTS\n")
            f.write("=" * 80 + "\n\n")
            f.write(partial_text + "\n")
        print(f"✓ Results saved to {output_file}")
    
    return full_text, partial_text


# ============================================================
# USAGE
# ============================================================
if __name__ == "__main__":
    group_folder = "/home/marc_be/disks/meso_shared/RetinoMaps/derivatives/pp_data/group/91k/rest"
    
    full_text, partial_text = format_manuscript_results(
        group_folder,
        run_variant="concat_clean",
        output_file="manuscript_results.txt"  # Optional: saves to file
    )
    
    print("\n" + "=" * 80)
    print("FULL CORRELATION")
    print("=" * 80)
    print(full_text)
    print("\n" + "=" * 80)
    print("PARTIAL CORRELATION")
    print("=" * 80)
    print(partial_text)

