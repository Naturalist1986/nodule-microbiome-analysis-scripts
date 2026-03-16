#!/usr/bin/env python3
"""
Create horizontal stacked bar chart for bacterial vs non-bacterial mapping.
Supports both 3-category (original) and 4-category (with MAGs) modes.
Style matches read_mapping_horizontal_stacked.svg example.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

script_dir = os.path.dirname(os.path.abspath(__file__))

# Try the MAGs version first, fall back to original
mags_tsv = os.path.join(script_dir, 'bacterial_vs_nonbacterial_with_mags_summary.tsv')
orig_tsv = os.path.join(script_dir, 'bacterial_vs_nonbacterial_summary.tsv')

if os.path.exists(mags_tsv):
    print(f'Reading: {mags_tsv}')
    df = pd.read_csv(mags_tsv, sep='\t')
    has_mags = True
else:
    print(f'MAGs TSV not found, using original: {orig_tsv}')
    df = pd.read_csv(orig_tsv, sep='\t')
    has_mags = False

treatment_names = {
    'carK': 'CAR-C',
    'carR': 'CAR-G',
    'ces':  'CES',
    'hok':  'HUK',
    'mtz':  'MTZ',
    'RH':   'RH',
}
df['Display_Name'] = df['Treatment'].map(treatment_names).fillna(df['Treatment'])

if has_mags:
    # Sort alphabetically by display name (descending so A appears at top in barh)
    df['Total_Bacterial_%'] = df['MAG_%'] + df['Other_Bacterial_%']
    df = df.sort_values('Display_Name', ascending=False)

    treatments = df['Display_Name'].values
    mags = df['MAG_%'].values
    other_bact = df['Other_Bacterial_%'].values
    nonbacterial = df['Non-Bacterial_%'].values
    unmapped = df['Unmapped_%'].values

    # Colors
    color_mags = '#D95F02'           # Dark orange for MAGs
    color_other_bact = '#FDAE6B'     # Light orange for other bacterial
    color_nonbacterial = '#6AAF6A'   # Green for non-bacterial
    color_unmapped = '#BABABA'       # Gray for unmapped

    fig, ax = plt.subplots(figsize=(10, 5))
    y_pos = np.arange(len(treatments))
    bar_height = 0.6

    # Plot 4 stacked segments
    bars1 = ax.barh(y_pos, mags, bar_height,
                    color=color_mags, label='MAGs (Selected Bins)',
                    edgecolor='white', linewidth=0.5)

    bars2 = ax.barh(y_pos, other_bact, bar_height,
                    left=mags,
                    color=color_other_bact, label='Other Bacterial',
                    edgecolor='white', linewidth=0.5)

    bars3 = ax.barh(y_pos, nonbacterial, bar_height,
                    left=mags + other_bact,
                    color=color_nonbacterial, label='Non-Bacterial',
                    edgecolor='white', linewidth=0.5)

    bars4 = ax.barh(y_pos, unmapped, bar_height,
                    left=mags + other_bact + nonbacterial,
                    color=color_unmapped, label='Unmapped',
                    edgecolor='white', linewidth=0.5)

    # Add percentage labels inside each segment
    # Use smaller font for narrow segments, skip if < 3%
    for i in range(len(treatments)):
        segments = [
            (0, mags[i], 'white'),
            (mags[i], other_bact[i], '#333333'),
            (mags[i] + other_bact[i], nonbacterial[i], 'white'),
            (mags[i] + other_bact[i] + nonbacterial[i], unmapped[i], '#333333'),
        ]
        for left, width, color in segments:
            if width < 3:
                continue
            x = left + width / 2
            fs = 7
            ax.text(x, y_pos[i], f'{width:.1f}%',
                    ha='center', va='center',
                    fontsize=fs, fontweight='bold', color=color)

    title = 'Read Mapping Distribution Across Treatments\n(MAGs | Other Bacterial | Non-Bacterial | Unmapped)'
    output_base = os.path.join(script_dir, 'bacterial_vs_nonbacterial_with_mags_stacked')

else:
    # Original 3-category mode
    df = df.sort_values('Display_Name', ascending=False)

    treatments = df['Display_Name'].values
    bacterial = df['Bacterial_%'].values
    nonbacterial = df['Non-Bacterial_%'].values
    unmapped = df['Unmapped_%'].values

    color_bacterial = '#E8853A'
    color_nonbacterial = '#6AAF6A'
    color_unmapped = '#BABABA'

    fig, ax = plt.subplots(figsize=(10, 5))
    y_pos = np.arange(len(treatments))
    bar_height = 0.6

    bars1 = ax.barh(y_pos, bacterial, bar_height,
                    color=color_bacterial, label='Bacterial (Bact)',
                    edgecolor='white', linewidth=0.5)

    bars2 = ax.barh(y_pos, nonbacterial, bar_height,
                    left=bacterial,
                    color=color_nonbacterial, label='Non-Bacterial',
                    edgecolor='white', linewidth=0.5)

    bars3 = ax.barh(y_pos, unmapped, bar_height,
                    left=bacterial + nonbacterial,
                    color=color_unmapped, label='Unmapped',
                    edgecolor='white', linewidth=0.5)

    for i in range(len(treatments)):
        if bacterial[i] > 8:
            ax.text(bacterial[i] / 2, y_pos[i],
                    f'{bacterial[i]:.1f}%',
                    ha='center', va='center',
                    fontsize=9, fontweight='bold', color='white')

        mid_nonbact = bacterial[i] + nonbacterial[i] / 2
        if nonbacterial[i] > 8:
            ax.text(mid_nonbact, y_pos[i],
                    f'{nonbacterial[i]:.1f}%',
                    ha='center', va='center',
                    fontsize=9, fontweight='bold', color='white')

        mid_unmapped = bacterial[i] + nonbacterial[i] + unmapped[i] / 2
        if unmapped[i] > 8:
            ax.text(mid_unmapped, y_pos[i],
                    f'{unmapped[i]:.1f}%',
                    ha='center', va='center',
                    fontsize=9, fontweight='bold', color='#333333')

    title = 'Read Mapping Distribution Across Treatments\n(Bacterial vs Non-Bacterial)'
    output_base = os.path.join(script_dir, 'bacterial_vs_nonbacterial_stacked')

# Common formatting
ax.set_yticks(y_pos)
ax.set_yticklabels(treatments, fontsize=11, fontweight='bold')
ax.set_xlabel('Percentage of Total Reads (%)', fontsize=12, fontweight='bold')
ax.set_ylabel('Treatment', fontsize=12, fontweight='bold')
ax.set_title(title, fontsize=13, fontweight='bold', pad=15)

ax.set_xlim(0, 100)
ax.xaxis.set_major_locator(plt.MultipleLocator(20))
ax.xaxis.grid(True, linestyle='-', alpha=0.3, color='#cccccc')
ax.set_axisbelow(True)

ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=4,
          fontsize=9, framealpha=0.9, edgecolor='#cccccc', fancybox=False)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

# Save in multiple formats
plt.savefig(f'{output_base}.png', dpi=300, bbox_inches='tight', facecolor='white')
print(f'Saved: {output_base}.png')

plt.savefig(f'{output_base}.svg', bbox_inches='tight', facecolor='white')
print(f'Saved: {output_base}.svg')

plt.savefig(f'{output_base}.pdf', bbox_inches='tight', facecolor='white')
print(f'Saved: {output_base}.pdf')

plt.close()
print('\nDone!')
