#!/usr/bin/env python3
"""
Ensemble simulation for discrete baryogenesis horizon model.
Runs 1000 simulations with varying seeds and E_0 factors to build statistics.
"""

import numpy as np
import matplotlib.pyplot as plt
from discrete_baryogenesis_horizon import run_simulation, print_summary

def run_ensemble(n_runs=1000, n_steps=6000, transition_m=500, lock_in_width=8.0):
    """Run ensemble of simulations and collect results."""
    results = []
    
    for i in range(n_runs):
        E_0_factor = 1
        # Vary seed
        seed = np.random.randint(0, 1000000)
        
        try:
            data = run_simulation(
                n_steps=n_steps,
                transition_m=transition_m,
                lock_in_sharpness=lock_in_width,
                E_0_factor=E_0_factor,
                seed=seed
            )
            results.append({
                'eta': data['eta'],
                'total_modes': data['total_modes'],
                'net_baryons': data['net_baryons'],
                'E_0_factor': E_0_factor,
                'seed': seed
            })
        except Exception as e:
            print(f"Simulation {i} failed: {e}")
            continue
    
    return results

def analyze_results(results):
    """Analyze ensemble results and compute statistics."""
    etas = np.array([r['eta'] for r in results])
    abs_etas = np.abs(etas)  # Track |η| for meaningful analysis
    total_modes = np.array([r['total_modes'] for r in results])
    net_baryons = np.array([r['net_baryons'] for r in results])
    E_0_factors = np.array([r['E_0_factor'] for r in results])
    
    # Statistics for |η|
    abs_eta_mean = np.mean(abs_etas)
    abs_eta_std = np.std(abs_etas)
    abs_eta_median = np.median(abs_etas)
    
    # Statistics for signed η
    eta_mean = np.mean(etas)
    eta_std = np.std(etas)
    eta_median = np.median(etas)
    
    total_modes_mean = np.mean(total_modes)
    total_modes_std = np.std(total_modes)
    
    net_baryons_mean = np.mean(net_baryons)
    net_baryons_std = np.std(net_baryons)
    
    print("Ensemble Statistics:")
    print(f"Number of successful runs: {len(results)}")
    print(f"|η| mean: {abs_eta_mean:.4e} ± {abs_eta_std:.4e}")
    print(f"|η| median: {abs_eta_median:.4e}")
    print(f"η mean: {eta_mean:.4e} ± {eta_std:.4e}")
    print(f"η median: {eta_median:.4e}")
    print(f"Total modes mean: {total_modes_mean:.4e} ± {total_modes_std:.4e}")
    print(f"Net baryons mean: {net_baryons_mean:.4e} ± {net_baryons_std:.4e}")
    
    return {
        'etas': etas,
        'abs_etas': abs_etas,
        'total_modes': total_modes,
        'net_baryons': net_baryons,
        'E_0_factors': E_0_factors,
        'eta_mean': eta_mean,
        'eta_std': eta_std,
        'eta_median': eta_median,
        'abs_eta_mean': abs_eta_mean,
        'abs_eta_std': abs_eta_std,
        'abs_eta_median': abs_eta_median
    }

def plot_results(results, output_file='ensemble_eta_dist.png'):
    """Create visualization of ensemble results."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. Histogram of |η| values (absolute values for meaningful distribution)
    axes[0, 0].hist(results['abs_etas'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    axes[0, 0].axvline(results['abs_eta_mean'], color='red', linestyle='--', 
                       label=f'|η| Mean: {results["abs_eta_mean"]:.2e} ± {results["abs_eta_std"]:.2e}')
    axes[0, 0].axvline(results['abs_eta_median'], color='green', linestyle='--', 
                       label=f'|η| Median: {results["abs_eta_median"]:.2e}')
    axes[0, 0].axvline(6.1e-10, color='orange', linestyle='-', 
                       label='Observed η: 6.1e-10')
    axes[0, 0].set_xlabel('|η| (absolute baryon-to-photon ratio)')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Distribution of |η| from Ensemble')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. |η| vs E_0 factor
    scatter = axes[0, 1].scatter(results['E_0_factors'], results['abs_etas'], 
                                alpha=0.6, c=results['abs_etas'], cmap='viridis')
    axes[0, 1].axhline(6.1e-10, color='orange', linestyle='-', 
                       label='Observed η: 6.1e-10')
    axes[0, 1].set_xlabel('E_0 factor (relative to T_Pl)')
    axes[0, 1].set_ylabel('|η|')
    axes[0, 1].set_title('|η| vs E_0 Factor')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=axes[0, 1], label='|η| value')
    
    # 3. Net baryons vs Total modes
    axes[1, 0].scatter(results['total_modes'], results['net_baryons'], 
                     alpha=0.6, color='purple')
    axes[1, 0].set_xlabel('Total modes')
    axes[1, 0].set_ylabel('Net baryons')
    axes[1, 0].set_title('Net Baryons vs Total Modes')
    axes[1, 0].grid(True, alpha=0.3)
    
    # 4. Convergence plot (running mean of |η|)
    running_mean_abs = np.cumsum(results['abs_etas']) / np.arange(1, len(results['abs_etas']) + 1)
    axes[1, 1].plot(running_mean_abs, color='blue', linewidth=2)
    axes[1, 1].axhline(6.1e-10, color='orange', linestyle='-', 
                       label='Observed η: 6.1e-10')
    axes[1, 1].set_xlabel('Simulation number')
    axes[1, 1].set_ylabel('Running mean of |η|')
    axes[1, 1].set_title('Convergence of |η| Mean')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plot saved to {output_file}")

def main():
    print("Running ensemble simulation with 1000 runs...")
    results = run_ensemble(n_runs=1000)
    
    print("\nAnalyzing results...")
    stats = analyze_results(results)
    
    print("\nGenerating visualization...")
    plot_results(stats)
    
    print("\nEnsemble simulation complete!")

if __name__ == "__main__":
    main()