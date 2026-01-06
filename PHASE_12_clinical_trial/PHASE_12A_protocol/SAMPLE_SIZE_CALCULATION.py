#!/usr/bin/env python3
"""
PHASE 12: Sample Size Calculation for Febuxostat Phase II Trial

Method: Two-proportion z-test
Primary comparison: Placebo vs Febuxostat 80mg

Assumptions:
- Alpha = 0.05 (two-sided)
- Power = 80%
- Placebo CDAI-100 response: 25%
- Febuxostat 80mg CDAI-100 response: 50%
- Dropout: 15%
"""

import math
import os

# Output directory
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

# Parameters
ALPHA = 0.05  # Two-sided
POWER = 0.80
P1 = 0.25     # Placebo response rate
P2 = 0.50     # Febuxostat 80mg response rate
DROPOUT = 0.15


def calculate_sample_size_per_arm(p1, p2, alpha, power):
    """
    Calculate sample size per arm using two-proportion z-test formula.

    N per arm = 2 * [(Z_alpha/2 + Z_beta)^2 * p_bar * (1 - p_bar)] / (p2 - p1)^2

    Where:
    - Z_alpha/2 = 1.96 for alpha=0.05 two-sided
    - Z_beta = 0.84 for power=0.80
    - p_bar = (p1 + p2) / 2
    """
    from scipy import stats

    # Z-scores
    z_alpha = stats.norm.ppf(1 - alpha/2)  # 1.96 for alpha=0.05
    z_beta = stats.norm.ppf(power)          # 0.84 for power=0.80

    # Pooled proportion
    p_bar = (p1 + p2) / 2

    # Effect size
    effect = p2 - p1

    # Sample size formula
    numerator = 2 * ((z_alpha + z_beta) ** 2) * p_bar * (1 - p_bar)
    denominator = effect ** 2

    n_per_arm = numerator / denominator

    return math.ceil(n_per_arm), z_alpha, z_beta


def main():
    print("=" * 60)
    print("PHASE 12: Sample Size Calculation")
    print("Febuxostat Phase II Trial in Crohn's Disease")
    print("=" * 60)

    print("\n## Assumptions")
    print(f"  Alpha (two-sided):     {ALPHA}")
    print(f"  Power:                 {POWER * 100:.0f}%")
    print(f"  Placebo response (p1): {P1 * 100:.0f}%")
    print(f"  Drug response (p2):    {P2 * 100:.0f}%")
    print(f"  Expected dropout:      {DROPOUT * 100:.0f}%")

    # Calculate
    n_raw, z_alpha, z_beta = calculate_sample_size_per_arm(P1, P2, ALPHA, POWER)

    print("\n## Z-scores")
    print(f"  Z_alpha/2: {z_alpha:.3f}")
    print(f"  Z_beta:    {z_beta:.3f}")

    print("\n## Results (2-arm comparison: Placebo vs 80mg)")
    print(f"  Raw N per arm:              {n_raw}")

    # Adjust for dropout
    n_adjusted = math.ceil(n_raw / (1 - DROPOUT))
    print(f"  Dropout-adjusted N per arm: {n_adjusted}")

    # Total for 2-arm
    total_2arm = n_adjusted * 2
    print(f"  Total N (2-arm):            {total_2arm}")

    # Total for 3-arm (add 40mg exploratory arm, same size)
    total_3arm = n_adjusted * 3
    print("\n## Results (3-arm trial: Placebo + 40mg + 80mg)")
    print(f"  N per arm:                  {n_adjusted}")
    print(f"  Total N (3-arm):            {total_3arm}")
    print(f"  Note: 40mg arm is exploratory (same size for balance)")

    # Effect size
    effect_size = P2 - P1
    nnt = 1 / effect_size
    print("\n## Effect Size")
    print(f"  Absolute risk reduction:    {effect_size * 100:.0f}%")
    print(f"  Number needed to treat:     {nnt:.1f}")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"\n  >>> ENROLL {total_3arm} PATIENTS TOTAL <<<")
    print(f"      ({n_adjusted} per arm x 3 arms)")
    print("\n  Arm 1: Placebo      (n={})".format(n_adjusted))
    print("  Arm 2: Febuxostat 40mg (n={}) [exploratory]".format(n_adjusted))
    print("  Arm 3: Febuxostat 80mg (n={}) [primary]".format(n_adjusted))

    # Save to file
    output_file = os.path.join(OUTPUT_DIR, "sample_size_output.txt")
    with open(output_file, 'w') as f:
        f.write("PHASE 12: Sample Size Calculation Results\n")
        f.write("=" * 50 + "\n\n")
        f.write("Assumptions:\n")
        f.write(f"  Alpha (two-sided): {ALPHA}\n")
        f.write(f"  Power: {POWER * 100:.0f}%\n")
        f.write(f"  Placebo response: {P1 * 100:.0f}%\n")
        f.write(f"  Febuxostat 80mg response: {P2 * 100:.0f}%\n")
        f.write(f"  Dropout rate: {DROPOUT * 100:.0f}%\n\n")
        f.write("Results:\n")
        f.write(f"  Raw N per arm: {n_raw}\n")
        f.write(f"  Dropout-adjusted N per arm: {n_adjusted}\n")
        f.write(f"  Total N (3-arm trial): {total_3arm}\n\n")
        f.write("Trial Arms:\n")
        f.write(f"  Arm 1: Placebo (n={n_adjusted})\n")
        f.write(f"  Arm 2: Febuxostat 40mg (n={n_adjusted}) [exploratory]\n")
        f.write(f"  Arm 3: Febuxostat 80mg (n={n_adjusted}) [primary]\n")

    print(f"\n  Saved to: {output_file}")


if __name__ == "__main__":
    main()
