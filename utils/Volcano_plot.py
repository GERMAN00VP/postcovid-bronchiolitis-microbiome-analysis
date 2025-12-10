# Volcano plot by germán

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

def reg(x,sig_col,lfc_col,threshold):
    if not x[sig_col]:
        return 'Not Significant'
    elif x[lfc_col]>threshold:
        return 'Increased Abundance'
    else:
        return "Decreased Abundance"



def volcano_plot(df, logFC_col, q_value_col, logFC_thresh, q_thresh=0.05, ax=None, legend=True, reverse_lfc=False, passed_ss=False, taxa=False):
    """
    Plots a volcano plot for taxa, including standard error on log-fold changes,
    and annotates significant points with non-overlapping labels using adjustText.

    Parameters:
    -----------
    df : DataFrame
        DataFrame containing the taxa as the index and the columns for log-fold change,
        q-values, and standard error.
    logFC_col : str
        Name of the column in the DataFrame representing the log-fold change (logFC, base e).
    q_value_col : str
        Name of the column in the DataFrame representing the q-values (false discovery rate adjusted p-values).
    std_err_col : str
        Name of the column in the DataFrame representing the standard error of the log-fold change.
    logFC_thresh : float
        Threshold for the log-fold change.
    q_thresh : float, optional, default=0.05
        Threshold for the q-value. Points below this threshold will be considered significant.
    reverse_lfc: bool, if True reverts the effect of the LFC (LFC*-1)
    
    Returns:
    --------
    ax : matplotlib Axes
        The axes object containing the plot.
    """
    # Transform the q-values to -log10(q)
    df['-log10(q_value)'] = -np.log10(df[q_value_col])

    if reverse_lfc:
        df[logFC_col] = df[logFC_col] * -1
    
    if passed_ss:
        passed_ss = 'passed_ss_'+ logFC_col.split("lfc_")[1]
        df = df.loc[df[passed_ss]] 

    # Define significant taxa based on logFC and q-value thresholds
    df['Significant'] = (df[logFC_col].abs() >= logFC_thresh) & (df[q_value_col] <= q_thresh)
    df['Regulation'] = df.apply(reg, sig_col='Significant', lfc_col=logFC_col, threshold=logFC_thresh, axis=1)

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    # Plot points
    sns.scatterplot(x=df[logFC_col], y=df['-log10(q_value)'],
                    hue=df['Regulation'],
                    s=50, palette={'Increased Abundance': 'red', 
                                   'Decreased Abundance': 'blue', 
                                   'Not Significant': 'gray'}, ax=ax)
    
    # Add threshold lines
    ax.axvline(x=logFC_thresh, color='black', linestyle='--', label=f'LogFC threshold = {logFC_thresh}', alpha=0.5)
    ax.axvline(x=-logFC_thresh, color='black', linestyle='--', alpha=0.5)
    ax.axhline(y=-np.log10(q_thresh), color='black', linestyle='--', label=f'q-value threshold = {q_thresh}', alpha=0.5)
    
    # Add labels for significant points
    if taxa:
        texts = []
        df_sig = df.loc[df.Significant].copy()
        
        # Sort by significance (most significant first)
        df_sig = df_sig.sort_values(by=['-log10(q_value)', logFC_col], 
                                   ascending=[False, False])
        
        # Get taxa names (last part after split)
        taxa_names = [str(name).split(";")[-1] for name in df_sig[taxa]]
        
        # Create text objects
        for idx, row in df_sig.iterrows():
            x = row[logFC_col]
            y = row['-log10(q_value)']
            # Use abbreviated names if too long
            label = taxa_names.pop(0)
            if len(label) > 25:
                label = label[:22] + "..."
            texts.append(ax.text(x, y, label, fontsize=9, fontweight=600,
                                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                                         edgecolor='gray', alpha=0.8)))

        # Adjust text positions with improved parameters
        if texts:  # Only run if there are texts to adjust
            try:
                adjust_text(
                    texts,
                    arrowprops=dict(arrowstyle='-', color='black', alpha=0.5, lw=0.5),
                    autoalign='xy',          # Align in both directions
                    expand_text=(1.2, 1.2),  # Expand text spacing
                    expand_points=(1.5, 1.5), # Expand from points
                    expand_axes=False,       # Don't expand axes
                    force_text=(0.5, 1.5),   # Force between text
                    force_points=(0.8, 2.0), # Force from points
                    lim=500,                 # Maximum iterations
                    precision=0.001,         # Precision for convergence
                    only_move={'points': 'xy', 'text': 'xy'},  # Allow movement in both directions
                    avoid_self=True,         # Avoid overlapping with itself
                    avoid_points=True,       # Avoid overlapping with points
                    avoid_text=True,         # Avoid overlapping with other text
                    text_from_text=True,     # Consider text-to-text distances
                    text_from_points=True,   # Consider text-to-point distances
                    va='center',             # Vertical alignment
                    ha='center',             # Horizontal alignment
                    ax=ax                    # Specify axis
                )
            except Exception as e:
                print(f"Warning: adjust_text encountered an issue: {e}")
                # Fallback: adjust manually with minimal spacing
                for text in texts:
                    text.set_rotation(45)
                    text.set_alpha(0.8)

    # Set labels and title
    ax.set_xlabel('Log-Fold Change (logFC, ln)', fontsize=12)
    ax.set_ylabel('-log10(q-value)', fontsize=12)
    ax.set_title('Volcano Plot with Upregulated and Downregulated Taxa', fontsize=14, pad=20)
    
    # Handle legend
    if legend:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    else:
        ax.legend().set_visible(False)

    # Adjust layout
    plt.tight_layout()
    
    # Return appropriate columns
    if taxa:    
        return ax, df[[logFC_col, '-log10(q_value)', 'Regulation', "Significant", taxa]].sort_values(by=logFC_col)
    
    return ax, df[[logFC_col, '-log10(q_value)', 'Regulation', "Significant"]].sort_values(by=logFC_col)