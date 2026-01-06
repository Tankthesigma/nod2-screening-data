#!/usr/bin/env python3
"""
PHASE 12B: Monte Carlo Clinical Trial Simulation

Simulates 1,000 virtual Phase II trials to predict probability of success
for Febuxostat in NOD2-positive Crohn's disease.

DISCLAIMER: This is a computational simulation for research/educational purposes.
Simulated patients are not real. Results predict statistical outcomes, not actual
clinical efficacy.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# SETUP & REPRODUCIBILITY
# =============================================================================

SEED = 42
rng = np.random.default_rng(SEED)

# Output directory (subfolder for results)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "results")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Trial parameters (from Phase 12A)
N_PER_ARM = 70
N_TOTAL = N_PER_ARM * 3
N_TRIALS = 1000
N_BOOTSTRAP = 2000

# Response rates (assumed true rates)
RESPONSE_RATES = {
    'placebo': 0.25,
    'febux_40mg': 0.40,
    'febux_80mg': 0.50
}

# Stratification proportions
PROP_WT = 0.60          # 60% wild-type, 40% missense
PROP_PRIOR_BIO = 0.30   # 30% prior biologic

# Response modifiers (additive)
MODIFIERS = {
    'genotype_wt': 0.02,       # WT: +2%
    'genotype_missense': -0.02, # Missense: -2%
    'prior_bio_no': 0.03,      # No prior: +3%
    'prior_bio_yes': -0.03,    # Prior: -3%
    'cdai_moderate': 0.02,     # 220-300: +2%
    'cdai_severe': -0.02       # 301-450: -2%
}

# =============================================================================
# PART A: PATIENT GENERATOR
# =============================================================================

def generate_patient_cohort(rng_instance):
    """
    Generate virtual patient cohort with stratification factors.

    Stratification factors (match Phase 12A):
    1. Genotype: 60% WT, 40% missense
    2. Baseline CDAI bin: 220-300 vs 301-450
    3. Prior biologic: 30% yes, 70% no
    """
    patients = []

    for i in range(N_TOTAL):
        # Generate baseline characteristics
        genotype = 'WT' if rng_instance.random() < PROP_WT else 'missense'
        prior_biologic = 'yes' if rng_instance.random() < PROP_PRIOR_BIO else 'no'
        baseline_cdai = rng_instance.uniform(220, 450)
        cdai_bin = 'moderate' if baseline_cdai <= 300 else 'severe'

        patients.append({
            'patient_id': i + 1,
            'genotype': genotype,
            'prior_biologic': prior_biologic,
            'baseline_cdai': baseline_cdai,
            'cdai_bin': cdai_bin
        })

    df = pd.DataFrame(patients)

    # Stratified block randomization
    df = stratified_block_randomization(df, rng_instance)

    return df


def stratified_block_randomization(df, rng_instance):
    """
    Assign treatment arms using stratified block randomization.

    IMPORTANT: Guarantees exactly N_PER_ARM (70) patients per arm.
    First assigns within strata using blocks, then balances any overflow.
    """
    # Create stratum identifier
    df['stratum'] = df['genotype'] + '_' + df['cdai_bin'] + '_' + df['prior_biologic']

    arms = ['placebo', 'febux_40mg', 'febux_80mg']
    df['arm'] = ''

    # First pass: assign within strata
    for stratum in df['stratum'].unique():
        stratum_mask = df['stratum'] == stratum
        stratum_indices = df.index[stratum_mask].tolist()

        # Shuffle patients within stratum
        rng_instance.shuffle(stratum_indices)

        # Assign arms in blocks of 3
        for i, idx in enumerate(stratum_indices):
            df.loc[idx, 'arm'] = arms[i % 3]

    # Second pass: ensure exactly N_PER_ARM per arm
    # Count current assignments
    arm_counts = df['arm'].value_counts().to_dict()

    for arm in arms:
        current_count = arm_counts.get(arm, 0)

        if current_count > N_PER_ARM:
            # Too many in this arm - reassign excess to underrepresented arms
            excess = current_count - N_PER_ARM
            arm_indices = df[df['arm'] == arm].index.tolist()
            rng_instance.shuffle(arm_indices)

            for i in range(excess):
                # Find arm with fewest patients
                arm_counts = df['arm'].value_counts().to_dict()
                min_arm = min(arms, key=lambda a: arm_counts.get(a, 0))
                if arm_counts.get(min_arm, 0) < N_PER_ARM:
                    df.loc[arm_indices[i], 'arm'] = min_arm

    return df


# =============================================================================
# PART B: RESPONSE SIMULATOR
# =============================================================================

def simulate_responses(df, rng_instance):
    """
    Simulate treatment response for each patient.

    Base response probability by arm, modified by patient characteristics.
    Probabilities clipped to [0.01, 0.99].
    """
    responses = []

    for _, patient in df.iterrows():
        # Base response probability
        base_prob = RESPONSE_RATES[patient['arm']]

        # Apply modifiers
        if patient['genotype'] == 'WT':
            base_prob += MODIFIERS['genotype_wt']
        else:
            base_prob += MODIFIERS['genotype_missense']

        if patient['prior_biologic'] == 'no':
            base_prob += MODIFIERS['prior_bio_no']
        else:
            base_prob += MODIFIERS['prior_bio_yes']

        if patient['cdai_bin'] == 'moderate':
            base_prob += MODIFIERS['cdai_moderate']
        else:
            base_prob += MODIFIERS['cdai_severe']

        # Clip probability to [0.01, 0.99]
        prob = np.clip(base_prob, 0.01, 0.99)

        # Simulate binary response
        response = 1 if rng_instance.random() < prob else 0
        responses.append(response)

    df['response'] = responses
    df['response_prob'] = df.apply(lambda row: calculate_patient_prob(row), axis=1)

    return df


def calculate_patient_prob(row):
    """Calculate true response probability for a patient (for reference)."""
    base_prob = RESPONSE_RATES[row['arm']]

    if row['genotype'] == 'WT':
        base_prob += MODIFIERS['genotype_wt']
    else:
        base_prob += MODIFIERS['genotype_missense']

    if row['prior_biologic'] == 'no':
        base_prob += MODIFIERS['prior_bio_no']
    else:
        base_prob += MODIFIERS['prior_bio_yes']

    if row['cdai_bin'] == 'moderate':
        base_prob += MODIFIERS['cdai_moderate']
    else:
        base_prob += MODIFIERS['cdai_severe']

    return np.clip(base_prob, 0.01, 0.99)


# =============================================================================
# PART C: SINGLE TRIAL ANALYZER
# =============================================================================

def analyze_single_trial(df, rng_instance):
    """
    Analyze outcomes from one simulated trial.

    Primary analysis: Logistic regression (Placebo vs 80mg)
    - Adjusted for baseline CDAI, genotype, prior biologic
    - P-value is for the treatment indicator coefficient

    Secondary: Chi-square tests for comparison
    """
    results = {}

    # Response rates per arm
    for arm in ['placebo', 'febux_40mg', 'febux_80mg']:
        arm_data = df[df['arm'] == arm]
        results[f'{arm}_n'] = len(arm_data)
        results[f'{arm}_responders'] = arm_data['response'].sum()
        results[f'{arm}_rate'] = arm_data['response'].mean()

    # Get response arrays for chi-square tests
    placebo_resp = df[df['arm'] == 'placebo']['response']
    febux80_resp = df[df['arm'] == 'febux_80mg']['response']
    febux40_resp = df[df['arm'] == 'febux_40mg']['response']

    # Secondary: Chi-square test (Placebo vs 80mg) - computed FIRST as fallback
    contingency_primary = np.array([
        [placebo_resp.sum(), len(placebo_resp) - placebo_resp.sum()],
        [febux80_resp.sum(), len(febux80_resp) - febux80_resp.sum()]
    ])

    try:
        chi2, p_chi2, _, _ = stats.chi2_contingency(contingency_primary, correction=False)
        results['secondary_chi2_pvalue'] = p_chi2
    except ValueError:
        results['secondary_chi2_pvalue'] = 1.0

    # Exploratory: Chi-square test (Placebo vs 40mg)
    contingency_exploratory = np.array([
        [placebo_resp.sum(), len(placebo_resp) - placebo_resp.sum()],
        [febux40_resp.sum(), len(febux40_resp) - febux40_resp.sum()]
    ])

    try:
        chi2_exp, p_chi2_exp, _, _ = stats.chi2_contingency(contingency_exploratory, correction=False)
        results['exploratory_chi2_pvalue'] = p_chi2_exp
    except ValueError:
        results['exploratory_chi2_pvalue'] = 1.0

    # Primary analysis: Logistic regression (Placebo vs 80mg)
    primary_df = df[df['arm'].isin(['placebo', 'febux_80mg'])].copy()
    primary_df['treatment'] = (primary_df['arm'] == 'febux_80mg').astype(int)
    primary_df['genotype_wt'] = (primary_df['genotype'] == 'WT').astype(int)
    primary_df['prior_bio_yes'] = (primary_df['prior_biologic'] == 'yes').astype(int)
    primary_df['cdai_severe'] = (primary_df['cdai_bin'] == 'severe').astype(int)

    # Fit logistic regression
    X = primary_df[['treatment', 'genotype_wt', 'prior_bio_yes', 'cdai_severe']]
    X = sm.add_constant(X)
    y = primary_df['response']

    try:
        model = sm.Logit(y, X).fit(disp=0, method='bfgs', maxiter=100)

        # Extract treatment coefficient results
        results['primary_pvalue'] = model.pvalues['treatment']
        results['primary_odds_ratio'] = np.exp(model.params['treatment'])
        results['primary_or_ci_low'] = np.exp(model.conf_int().loc['treatment', 0])
        results['primary_or_ci_high'] = np.exp(model.conf_int().loc['treatment', 1])
        results['primary_converged'] = True
    except (np.linalg.LinAlgError, ValueError, KeyError):
        # Fallback if model fails to converge (separation, singular matrix, etc.)
        # Use chi-square p-value as backup for this trial
        results['primary_pvalue'] = results['secondary_chi2_pvalue']
        results['primary_odds_ratio'] = np.nan
        results['primary_or_ci_low'] = np.nan
        results['primary_or_ci_high'] = np.nan
        results['primary_converged'] = False

    # Effect size: Risk difference (80mg - placebo)
    results['risk_difference'] = results['febux_80mg_rate'] - results['placebo_rate']

    # Bootstrap CI for risk difference
    rd_bootstrap = bootstrap_risk_difference(primary_df, rng_instance)
    results['risk_diff_ci_low'] = rd_bootstrap['ci_low']
    results['risk_diff_ci_high'] = rd_bootstrap['ci_high']

    # Effect size for 40mg (exploratory)
    results['risk_difference_40mg'] = results['febux_40mg_rate'] - results['placebo_rate']

    return results


def bootstrap_risk_difference(df, rng_instance, n_bootstrap=N_BOOTSTRAP):
    """
    Bootstrap 95% CI for risk difference using 2000 resamples.
    """
    placebo_df = df[df['arm'] == 'placebo']['response'].values
    treatment_df = df[df['arm'] == 'febux_80mg']['response'].values

    risk_diffs = []

    for _ in range(n_bootstrap):
        # Resample with replacement
        placebo_sample = rng_instance.choice(placebo_df, size=len(placebo_df), replace=True)
        treatment_sample = rng_instance.choice(treatment_df, size=len(treatment_df), replace=True)

        rd = treatment_sample.mean() - placebo_sample.mean()
        risk_diffs.append(rd)

    return {
        'ci_low': np.percentile(risk_diffs, 2.5),
        'ci_high': np.percentile(risk_diffs, 97.5)
    }


# =============================================================================
# PART D: MONTE CARLO ENGINE
# =============================================================================

def run_monte_carlo(n_trials=N_TRIALS):
    """
    Run Monte Carlo simulation of n_trials Phase II trials.

    Uses single RNG for reproducibility.
    """
    print(f"\nRunning {n_trials} simulated trials...")
    print(f"Random seed: {SEED}")
    print("-" * 50)

    all_results = []
    example_patients = None

    for trial_num in range(n_trials):
        if (trial_num + 1) % 100 == 0:
            print(f"  Trial {trial_num + 1}/{n_trials}")

        # Generate patients
        patients = generate_patient_cohort(rng)

        # Simulate responses
        patients = simulate_responses(patients, rng)

        # Save first trial's patient data as example
        if trial_num == 0:
            example_patients = patients.copy()

        # Analyze trial
        trial_results = analyze_single_trial(patients, rng)
        trial_results['trial_num'] = trial_num + 1
        all_results.append(trial_results)

    results_df = pd.DataFrame(all_results)

    return results_df, example_patients


def calculate_summary_stats(results_df):
    """
    Calculate summary statistics from Monte Carlo results.
    """
    summary = {}

    # Empirical power (% with p < 0.05)
    summary['empirical_power'] = (results_df['primary_pvalue'] < 0.05).mean() * 100
    summary['strong_evidence'] = (results_df['primary_pvalue'] < 0.01).mean() * 100

    # Convergence rate
    summary['convergence_rate'] = results_df['primary_converged'].mean() * 100

    # Response rates with 95% CI
    for arm in ['placebo', 'febux_40mg', 'febux_80mg']:
        rates = results_df[f'{arm}_rate']
        summary[f'{arm}_mean'] = rates.mean() * 100
        summary[f'{arm}_ci_low'] = np.percentile(rates, 2.5) * 100
        summary[f'{arm}_ci_high'] = np.percentile(rates, 97.5) * 100

    # Effect size (risk difference)
    rd = results_df['risk_difference']
    summary['risk_diff_mean'] = rd.mean() * 100
    summary['risk_diff_ci_low'] = np.percentile(rd, 2.5) * 100
    summary['risk_diff_ci_high'] = np.percentile(rd, 97.5) * 100

    # Odds ratio
    or_vals = results_df['primary_odds_ratio'].dropna()
    summary['odds_ratio_mean'] = or_vals.mean()
    summary['odds_ratio_ci_low'] = np.percentile(or_vals, 2.5)
    summary['odds_ratio_ci_high'] = np.percentile(or_vals, 97.5)

    # P-value distribution stats
    summary['pvalue_median'] = results_df['primary_pvalue'].median()
    summary['pvalue_mean'] = results_df['primary_pvalue'].mean()

    return summary


# =============================================================================
# PART E: POWER CURVE ANALYSIS
# =============================================================================

def generate_power_curve(effect_sizes=[0.10, 0.15, 0.20, 0.25, 0.30, 0.35], n_trials_per=500):
    """
    Generate power curve: empirical power vs assumed effect size.

    For each effect size, runs 500 trials and calculates % with p < 0.05.
    """
    print(f"\nGenerating power curve ({n_trials_per} trials per effect size)...")

    power_data = []

    for effect_size in effect_sizes:
        print(f"  Effect size: {effect_size*100:.0f}%...", end=" ")

        # Temporarily modify response rates
        original_80mg = RESPONSE_RATES['febux_80mg']
        RESPONSE_RATES['febux_80mg'] = RESPONSE_RATES['placebo'] + effect_size

        # Create new RNG for this power analysis
        power_rng = np.random.default_rng(SEED + int(effect_size * 1000))

        successes = 0
        for _ in range(n_trials_per):
            patients = generate_patient_cohort(power_rng)
            patients = simulate_responses(patients, power_rng)
            trial_results = analyze_single_trial(patients, power_rng)

            if trial_results['primary_pvalue'] < 0.05:
                successes += 1

        empirical_power = successes / n_trials_per * 100
        print(f"Power: {empirical_power:.1f}%")

        power_data.append({
            'effect_size': effect_size * 100,
            'empirical_power': empirical_power
        })

        # Restore original rate
        RESPONSE_RATES['febux_80mg'] = original_80mg

    return pd.DataFrame(power_data)


# =============================================================================
# PART F: OUTPUT
# =============================================================================

def save_results(results_df, example_patients, summary, power_curve_df):
    """
    Save all results to output files.
    """
    print("\nSaving results...")

    # 1. Simulation summary
    summary_file = os.path.join(OUTPUT_DIR, "simulation_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("PHASE 12B: Monte Carlo Clinical Trial Simulation Results\n")
        f.write("=" * 60 + "\n\n")

        f.write("REPRODUCIBILITY\n")
        f.write("-" * 40 + "\n")
        f.write(f"Random seed: {SEED}\n")
        f.write(f"Number of trials: {N_TRIALS}\n")
        f.write(f"Patients per trial: {N_TOTAL}\n")
        f.write(f"Bootstrap resamples: {N_BOOTSTRAP}\n\n")

        f.write("ASSUMED RESPONSE RATES\n")
        f.write("-" * 40 + "\n")
        f.write(f"Placebo:        {RESPONSE_RATES['placebo']*100:.0f}%\n")
        f.write(f"Febuxostat 40mg: {RESPONSE_RATES['febux_40mg']*100:.0f}%\n")
        f.write(f"Febuxostat 80mg: {RESPONSE_RATES['febux_80mg']*100:.0f}%\n\n")

        f.write("KEY RESULTS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Empirical Power (p < 0.05): {summary['empirical_power']:.1f}%\n")
        f.write(f"Strong Evidence (p < 0.01): {summary['strong_evidence']:.1f}%\n")
        f.write(f"Model Convergence Rate:     {summary['convergence_rate']:.1f}%\n\n")

        f.write("MEAN RESPONSE RATES (95% CI)\n")
        f.write("-" * 40 + "\n")
        f.write(f"Placebo:        {summary['placebo_mean']:.1f}% ({summary['placebo_ci_low']:.1f}-{summary['placebo_ci_high']:.1f}%)\n")
        f.write(f"Febuxostat 40mg: {summary['febux_40mg_mean']:.1f}% ({summary['febux_40mg_ci_low']:.1f}-{summary['febux_40mg_ci_high']:.1f}%)\n")
        f.write(f"Febuxostat 80mg: {summary['febux_80mg_mean']:.1f}% ({summary['febux_80mg_ci_low']:.1f}-{summary['febux_80mg_ci_high']:.1f}%)\n\n")

        f.write("EFFECT SIZE (RISK DIFFERENCE: 80mg - Placebo)\n")
        f.write("-" * 40 + "\n")
        f.write(f"Mean: {summary['risk_diff_mean']:.1f} percentage points\n")
        f.write(f"95% CI: {summary['risk_diff_ci_low']:.1f} to {summary['risk_diff_ci_high']:.1f} percentage points\n\n")

        f.write("ODDS RATIO (80mg vs Placebo)\n")
        f.write("-" * 40 + "\n")
        f.write(f"Mean: {summary['odds_ratio_mean']:.2f}\n")
        f.write(f"95% CI: {summary['odds_ratio_ci_low']:.2f} to {summary['odds_ratio_ci_high']:.2f}\n\n")

        f.write("INTERPRETATION\n")
        f.write("-" * 40 + "\n")
        f.write(f"Based on {N_TRIALS} simulated trials, a Phase II trial of {N_TOTAL} patients\n")
        f.write(f"has approximately {summary['empirical_power']:.0f}% probability of detecting a\n")
        f.write(f"statistically significant benefit for Febuxostat 80mg over placebo.\n")

    print(f"  Saved: {summary_file}")

    # 2. All trials data
    trials_file = os.path.join(OUTPUT_DIR, "all_trials_data.csv")
    results_df.to_csv(trials_file, index=False)
    print(f"  Saved: {trials_file}")

    # 3. Example patient data
    patients_file = os.path.join(OUTPUT_DIR, "patient_example.csv")
    example_patients.to_csv(patients_file, index=False)
    print(f"  Saved: {patients_file}")

    # 4. Power curve data
    power_file = os.path.join(OUTPUT_DIR, "power_curve_data.csv")
    power_curve_df.to_csv(power_file, index=False)
    print(f"  Saved: {power_file}")


def print_summary(summary):
    """
    Print key results to terminal.
    """
    print("\n" + "=" * 60)
    print("MONTE CARLO SIMULATION RESULTS")
    print("=" * 60)

    print(f"\nRandom Seed: {SEED}")
    print(f"Trials Simulated: {N_TRIALS}")

    print(f"\n>>> EMPIRICAL POWER: {summary['empirical_power']:.1f}% <<<")
    print(f"    (% of trials with p < 0.05 for primary comparison)")

    print(f"\nStrong Evidence (p < 0.01): {summary['strong_evidence']:.1f}%")

    print(f"\nMean Response Rates (95% CI):")
    print(f"  Placebo:        {summary['placebo_mean']:.1f}% ({summary['placebo_ci_low']:.1f}-{summary['placebo_ci_high']:.1f}%)")
    print(f"  Febuxostat 40mg: {summary['febux_40mg_mean']:.1f}% ({summary['febux_40mg_ci_low']:.1f}-{summary['febux_40mg_ci_high']:.1f}%)")
    print(f"  Febuxostat 80mg: {summary['febux_80mg_mean']:.1f}% ({summary['febux_80mg_ci_low']:.1f}-{summary['febux_80mg_ci_high']:.1f}%)")

    print(f"\nMean Effect Size (Risk Difference):")
    print(f"  {summary['risk_diff_mean']:.1f} percentage points (95% CI: {summary['risk_diff_ci_low']:.1f} to {summary['risk_diff_ci_high']:.1f})")

    print(f"\nMean Odds Ratio (80mg vs Placebo):")
    print(f"  {summary['odds_ratio_mean']:.2f} (95% CI: {summary['odds_ratio_ci_low']:.2f} to {summary['odds_ratio_ci_high']:.2f})")

    print("\n" + "=" * 60)


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 60)
    print("PHASE 12B: Monte Carlo Clinical Trial Simulation")
    print("=" * 60)
    print("\nDISCLAIMER: This is a computational simulation for research/")
    print("educational purposes. Simulated patients are not real.")

    # Run main Monte Carlo simulation
    results_df, example_patients = run_monte_carlo(N_TRIALS)

    # Calculate summary statistics
    summary = calculate_summary_stats(results_df)

    # Generate power curve
    power_curve_df = generate_power_curve()

    # Save all results
    save_results(results_df, example_patients, summary, power_curve_df)

    # Print summary to terminal
    print_summary(summary)

    print("\nSimulation complete!")
    print(f"Results saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
