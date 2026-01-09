# G908R_febuxostat_rep1 Low Occupancy Investigation

## Summary

| Metric | Rep1 | Rep2 | Rep3 |
|--------|------|------|------|
| Pocket Occupancy | 45.8% | 100.0% | 100.0% |
| Mean Distance to Pocket | ~8 Å (after unbinding) | ~2.5 Å | ~2.5 Å |
| Starting Position | IN pocket (2.72 Å) | IN pocket (2.24 Å) | IN pocket (2.58 Å) |
| Temperature | 311.13 ± 0.64 K | 311.16 ± 0.66 K | 311.15 ± 0.66 K |
| Potential Energy | -3,519,951 kJ/mol | -3,519,787 kJ/mol | -3,519,959 kJ/mol |

## Starting Structure Comparison

All three replicates started with febuxostat properly positioned IN the binding pocket:

- **Rep1**: Minimum distance to pocket residues = 2.72 Å
- **Rep2**: Minimum distance to pocket residues = 2.24 Å
- **Rep3**: Minimum distance to pocket residues = 2.58 Å

**Conclusion**: Starting pose is NOT the cause of unbinding. All replicates began with essentially identical, well-docked configurations.

## Time-Resolved Occupancy Analysis

### Unbinding Timeline (Rep1)
- **0-4 ns**: Ligand remains stably bound (distance ~2.5-4 Å)
- **~4.2 ns (frame 418)**: Unbinding event begins
- **4.2-10 ns**: Ligand drifts progressively further from pocket
- **Final state**: Ligand ~16 Å from starting position

### Unbinding Characteristics
- **Type**: Gradual diffusion, NOT sudden ejection
- **Pattern**: Smooth increase in distance over ~2 ns
- **No reversal**: Once unbound, ligand did not rebind during remaining simulation

### Rep2 & Rep3 Behavior
- Both maintained tight binding throughout entire 10 ns
- Maximum distance never exceeded 2.81 Å
- No unbinding attempts observed

## Distance Plot

![Distance comparison plot](plots/rep1_investigation_distance.png)

The plot clearly shows Rep1 (blue) diverging from Rep2 (orange) and Rep3 (green) around frame 400-450 (~4-4.5 ns), while Rep2 and Rep3 remain tightly bound throughout.

## PBC Artifact Check

Periodic Boundary Condition (PBC) analysis confirmed this is a **genuine unbinding event**:

| Check | Result | Interpretation |
|-------|--------|----------------|
| Max frame-to-frame jump | 3.14 Å | Normal thermal motion |
| Ligand in box | YES | No image jumping |
| Total displacement | 16.1 Å | Real diffusion path |
| Smooth trajectory | YES | No discontinuities |

**Conclusion**: No PBC artifacts detected. The unbinding is physically real.

## Simulation Stability Check

All three replicates showed excellent thermodynamic stability:

| Parameter | Rep1 | Rep2 | Rep3 | Expected |
|-----------|------|------|------|----------|
| Mean Temp (K) | 311.13 | 311.16 | 311.15 | 310 ± 2 |
| Temp StdDev | 0.64 | 0.66 | 0.66 | < 1 |
| Stability | STABLE | STABLE | STABLE | - |

No temperature spikes, energy drift, or other simulation artifacts that could explain the unbinding.

## Seed/Random Number Analysis

Simulation log files do not contain explicit random seed values. However, different random seeds are standard practice for MD replicates and are the expected source of trajectory divergence.

## Root Cause Determination

### Verdict: STOCHASTIC UNBINDING EVENT

This is a **genuine, biologically meaningful unbinding event** caused by natural thermal fluctuations in the molecular dynamics simulation.

### Evidence Supporting This Conclusion:

1. **Starting structures identical**: All reps began with ligand properly docked
2. **No simulation artifacts**: Temperature stable, no PBC issues
3. **Gradual unbinding**: Smooth diffusion, not sudden ejection
4. **Consistent with MD theory**: Different random seeds explore different conformational space
5. **Rep2/Rep3 stable**: Demonstrates binding CAN be maintained; Rep1 sampled an escape pathway

### Biological Interpretation

The G908R mutation may create a **metastable binding state** where:
- The ligand can bind (as shown by Rep2, Rep3)
- But escape pathways exist that can be accessed through thermal motion
- Rep1 happened to sample one such escape pathway

This is valuable information: it suggests febuxostat binding to G908R-NOD2 may be kinetically less stable than to wild-type, even if the thermodynamic binding affinity is similar.

## Conclusion

The low occupancy in Rep1 (45.8%) is **NOT an error or artifact**. It represents:

1. **Valid sampling** of an alternative trajectory
2. **Evidence of binding instability** for this mutation-ligand pair
3. **Important biological signal** about G908R-febuxostat interaction dynamics

The fact that 2/3 replicates maintained 100% occupancy while 1/3 showed unbinding suggests:
- **Binding is possible** (reps 2, 3)
- **Binding is not fully stable** (rep 1 escaped)
- **More replicates needed** for reliable residence time estimates

## Recommendations

1. **Keep all three replicates** in the dataset - Rep1 provides valuable information about binding stability

2. **Report with nuance**: For G908R_febuxostat, report:
   - Mean occupancy: 81.9% ± 31.3%
   - Or: "2/3 replicates showed stable binding; 1/3 showed unbinding at ~4 ns"

3. **Consider extended simulations**: If residence time is important, run longer simulations (50-100 ns) or use enhanced sampling methods

4. **Compare to wild-type**: Does WT_febuxostat show any unbinding in replicates? If not, this supports G908R causing binding instability

5. **Binding kinetics study**: Consider umbrella sampling or metadynamics to quantify the binding free energy landscape

---

*Investigation completed: 2026-01-09*
*Analysis pipeline: MDAnalysis with PBC-corrected distance calculations*
