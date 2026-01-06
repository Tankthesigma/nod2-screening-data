#!/usr/bin/env python3
"""
PHASE 12B: Monte Carlo Visualization Generator

Generates figures for Febuxostat Monte Carlo simulation results.

Figures:
1. Response Rate Bar Chart
2. P-Value Distribution
3. Effect Size Distribution
4. Power Curve
5. Forest Plot
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

# =============================================================================
# SETUP
# =============================================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")
OUTPUT_DIR = RESULTS_DIR  # Save figures alongside data

# Load data
TRIALS_FILE = os.path.join(RESULTS_DIR, "all_trials_data.csv")
POWER_FILE = os.path.join(RESULTS_DIR, "power_curve_data.csv")

# Style settings
plt.style.use('seaborn-v0_8-whitegrid')
COLORS = {
    'placebo': '#808080',      # Gray
    'febux_40mg': '#3498db',   # Blue
    'febux_80mg': '#e74c3c',   # Red
    'significant': '#27ae60',   # Green
    'nonsignificant': '#bdc3c7' # Light gray
}


# =============================================================================
# FIGURE 1: RESPONSE RATE BAR CHART
# =============================================================================

def create_response_rate_chart(df):
    """Bar chart of mean response rates with 95% CI error bars."""

    fig, ax = plt.subplots(figsize=(10, 6))

    # Calculate statistics
    arms = ['placebo', 'febux_40mg', 'febux_80mg']
    labels = ['Placebo', 'Febuxostat\n40mg', 'Febuxostat\n80mg']
    colors = [COLORS['placebo'], COLORS['febux_40mg'], COLORS['febux_80mg']]

    means = []
    ci_lows = []
    ci_highs = []

    for arm in arms:
        rates = df[f'{arm}_rate'] * 100
        mean = rates.mean()
        ci_low = np.percentile(rates, 2.5)
        ci_high = np.percentile(rates, 97.5)
        means.append(mean)
        ci_lows.append(mean - ci_low)
        ci_highs.append(ci_high - mean)

    x = np.arange(len(arms))
    bars = ax.bar(x, means, color=colors, edgecolor='black', linewidth=1.5, width=0.6)

    # Error bars
    ax.errorbar(x, means, yerr=[ci_lows, ci_highs], fmt='none', color='black',
                capsize=8, capthick=2, elinewidth=2)

    # Add value labels
    for i, (bar, mean) in enumerate(zip(bars, means)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + ci_highs[i] + 2,
                f'{mean:.1f}%', ha='center', va='bottom', fontsize=12, fontweight='bold')

    # Significance annotation
    empirical_power = (df['primary_pvalue'] < 0.05).mean() * 100
    if empirical_power > 50:
        ax.annotate('', xy=(0, means[0] + 5), xytext=(2, means[2] + 5),
                    arrowprops=dict(arrowstyle='-', color='black', lw=1.5))
        ax.text(1, max(means) + 15, f'p < 0.05 in {empirical_power:.0f}% of trials',
                ha='center', fontsize=11, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=12)
    ax.set_ylabel('CDAI-100 Response Rate (%)', fontsize=12)
    ax.set_ylim(0, max(means) + max(ci_highs) + 25)
    ax.set_title('Simulated Response Rates Across 1,000 Trials\n(Febuxostat Phase II)',
                 fontsize=14, fontweight='bold')

    # Add note
    ax.text(0.02, 0.98, 'Error bars: 95% CI across trials', transform=ax.transAxes,
            fontsize=9, va='top', style='italic')

    plt.tight_layout()
    output_file = os.path.join(OUTPUT_DIR, "response_rates_barchart.png")
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"  Created: {output_file}")
    return output_file


# =============================================================================
# FIGURE 2: P-VALUE DISTRIBUTION
# =============================================================================

def create_pvalue_distribution(df):
    """Histogram of p-values with power annotation."""

    fig, ax = plt.subplots(figsize=(10, 6))

    pvalues = df['primary_pvalue']

    # Create histogram
    n, bins, patches = ax.hist(pvalues, bins=50, range=(0, 1),
                                edgecolor='black', linewidth=0.5)

    # Color bars by significance
    for i, patch in enumerate(patches):
        if bins[i] < 0.05:
            patch.set_facecolor(COLORS['significant'])
        else:
            patch.set_facecolor(COLORS['nonsignificant'])

    # Add vertical line at p = 0.05
    ax.axvline(x=0.05, color='red', linestyle='--', linewidth=2, label='α = 0.05')

    # Calculate power
    empirical_power = (pvalues < 0.05).mean() * 100
    strong_evidence = (pvalues < 0.01).mean() * 100

    # Add power annotation box
    textstr = f'Empirical Power: {empirical_power:.1f}%\n(p < 0.05 in {int(empirical_power*10)} of 1000 trials)'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.95, 0.95, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', horizontalalignment='right', bbox=props)

    ax.set_xlabel('P-value (Primary Logistic Regression)', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title('Distribution of P-Values Across 1,000 Simulated Trials\n(Placebo vs Febuxostat 80mg)',
                 fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')

    # Add legend for colors
    sig_patch = mpatches.Patch(color=COLORS['significant'], label=f'p < 0.05 ({empirical_power:.0f}%)')
    nonsig_patch = mpatches.Patch(color=COLORS['nonsignificant'], label=f'p ≥ 0.05 ({100-empirical_power:.0f}%)')
    ax.legend(handles=[sig_patch, nonsig_patch], loc='upper right')

    plt.tight_layout()
    output_file = os.path.join(OUTPUT_DIR, "pvalue_distribution.png")
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"  Created: {output_file}")
    return output_file


# =============================================================================
# FIGURE 3: EFFECT SIZE DISTRIBUTION
# =============================================================================

def create_effect_size_distribution(df):
    """Histogram of risk differences."""

    fig, ax = plt.subplots(figsize=(10, 6))

    effect_sizes = df['risk_difference'] * 100  # Convert to percentage points

    # Create histogram
    n, bins, patches = ax.hist(effect_sizes, bins=40, color=COLORS['febux_80mg'],
                                edgecolor='black', linewidth=0.5, alpha=0.7)

    # Mean and CI
    mean_effect = effect_sizes.mean()
    ci_low = np.percentile(effect_sizes, 2.5)
    ci_high = np.percentile(effect_sizes, 97.5)

    # Add vertical lines
    ax.axvline(x=mean_effect, color='darkred', linestyle='-', linewidth=2.5,
               label=f'Mean: {mean_effect:.1f} pp')
    ax.axvline(x=ci_low, color='darkred', linestyle='--', linewidth=1.5, alpha=0.7)
    ax.axvline(x=ci_high, color='darkred', linestyle='--', linewidth=1.5, alpha=0.7)

    # Shade CI region
    ax.axvspan(ci_low, ci_high, alpha=0.2, color='red', label=f'95% CI: [{ci_low:.1f}, {ci_high:.1f}]')

    # Add zero line
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1, alpha=0.5)

    ax.set_xlabel('Effect Size (Risk Difference, percentage points)', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title('Distribution of Effect Sizes Across 1,000 Simulated Trials\n(Febuxostat 80mg - Placebo)',
                 fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')

    plt.tight_layout()
    output_file = os.path.join(OUTPUT_DIR, "effect_size_distribution.png")
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"  Created: {output_file}")
    return output_file


# =============================================================================
# FIGURE 4: POWER CURVE
# =============================================================================

def create_power_curve(power_df):
    """Line plot of empirical power vs effect size."""

    fig, ax = plt.subplots(figsize=(10, 6))

    effect_sizes = power_df['effect_size']
    powers = power_df['empirical_power']

    # Plot power curve
    ax.plot(effect_sizes, powers, 'o-', color=COLORS['febux_80mg'], linewidth=2.5,
            markersize=10, markerfacecolor='white', markeredgewidth=2)

    # Add horizontal line at 80% power
    ax.axhline(y=80, color='gray', linestyle='--', linewidth=1.5, label='80% Power Target')

    # Add vertical line at assumed effect (25%)
    ax.axvline(x=25, color='green', linestyle='--', linewidth=1.5, label='Assumed Effect (25%)')

    # Find and mark intersection
    # Interpolate to find where power = 80%
    for i in range(len(effect_sizes)-1):
        if powers.iloc[i] < 80 <= powers.iloc[i+1]:
            # Linear interpolation
            x_intercept = effect_sizes.iloc[i] + (80 - powers.iloc[i]) * \
                         (effect_sizes.iloc[i+1] - effect_sizes.iloc[i]) / \
                         (powers.iloc[i+1] - powers.iloc[i])
            ax.plot(x_intercept, 80, 'ko', markersize=12)
            ax.annotate(f'MDE ≈ {x_intercept:.0f}%', xy=(x_intercept, 80),
                       xytext=(x_intercept + 3, 85), fontsize=10,
                       arrowprops=dict(arrowstyle='->', color='black'))

    # Mark assumed effect point
    assumed_power = powers[effect_sizes == 25].values[0]
    ax.plot(25, assumed_power, 's', color='green', markersize=12)
    ax.annotate(f'{assumed_power:.0f}%', xy=(25, assumed_power),
               xytext=(25 + 2, assumed_power - 5), fontsize=11, fontweight='bold')

    ax.set_xlabel('True Effect Size (percentage points)', fontsize=12)
    ax.set_ylabel('Empirical Power (%)', fontsize=12)
    ax.set_title('Power Curve: Empirical Power vs Effect Size\n(500 trials per point)',
                 fontsize=14, fontweight='bold')
    ax.set_xlim(5, 40)
    ax.set_ylim(0, 105)
    ax.legend(loc='lower right')

    # Add grid
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    output_file = os.path.join(OUTPUT_DIR, "power_curve.png")
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"  Created: {output_file}")
    return output_file


# =============================================================================
# FIGURE 5: FOREST PLOT
# =============================================================================

def create_forest_plot(df):
    """Forest plot of effect sizes for both comparisons."""

    fig, ax = plt.subplots(figsize=(10, 5))

    # Calculate effect sizes and CIs
    comparisons = [
        ('Placebo vs 40mg\n(Exploratory)', 'febux_40mg', 'placebo'),
        ('Placebo vs 80mg\n(Primary)', 'febux_80mg', 'placebo')
    ]

    y_positions = [1, 0]

    for i, (label, drug_arm, control_arm) in enumerate(comparisons):
        # Risk difference
        rd = (df[f'{drug_arm}_rate'] - df[f'{control_arm}_rate']) * 100
        mean_rd = rd.mean()
        ci_low = np.percentile(rd, 2.5)
        ci_high = np.percentile(rd, 97.5)

        y = y_positions[i]

        # Plot CI line
        ax.plot([ci_low, ci_high], [y, y], color='black', linewidth=2)

        # Plot point estimate
        color = COLORS['febux_80mg'] if '80mg' in label else COLORS['febux_40mg']
        ax.plot(mean_rd, y, 's', color=color, markersize=15, markeredgecolor='black', markeredgewidth=1)

        # Add text annotation
        ax.text(ci_high + 1, y, f'{mean_rd:.1f} [{ci_low:.1f}, {ci_high:.1f}]',
                va='center', fontsize=10)

    # Add vertical line at 0
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)

    # Add vertical line at no effect
    ax.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)

    # Labels
    ax.set_yticks(y_positions)
    ax.set_yticklabels([c[0] for c in comparisons], fontsize=11)
    ax.set_xlabel('Risk Difference (percentage points)', fontsize=12)
    ax.set_title('Effect Sizes by Comparison\n(Mean and 95% CI across 1,000 trials)',
                 fontsize=14, fontweight='bold')

    # Set limits
    ax.set_xlim(-10, 55)
    ax.set_ylim(-0.5, 1.5)

    # Add note
    ax.text(0.02, 0.02, 'Positive values favor treatment over placebo',
            transform=ax.transAxes, fontsize=9, style='italic')

    plt.tight_layout()
    output_file = os.path.join(OUTPUT_DIR, "forest_plot.png")
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"  Created: {output_file}")
    return output_file


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 60)
    print("PHASE 12B: Monte Carlo Visualization Generator")
    print("=" * 60)

    # Load data
    print("\nLoading data...")
    trials_df = pd.read_csv(TRIALS_FILE)
    power_df = pd.read_csv(POWER_FILE)
    print(f"  Loaded {len(trials_df)} trial results")
    print(f"  Loaded {len(power_df)} power curve points")

    # Generate figures
    print("\nGenerating figures...")

    print("\n1. Response Rate Bar Chart")
    create_response_rate_chart(trials_df)

    print("\n2. P-Value Distribution")
    create_pvalue_distribution(trials_df)

    print("\n3. Effect Size Distribution")
    create_effect_size_distribution(trials_df)

    print("\n4. Power Curve")
    create_power_curve(power_df)

    print("\n5. Forest Plot")
    create_forest_plot(trials_df)

    print("\n" + "=" * 60)
    print("COMPLETE")
    print("=" * 60)
    print(f"\nAll figures saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
