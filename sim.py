import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import time

# ────────────────────────────────────────────────
# Global physics parameters — final tuned values
# ────────────────────────────────────────────────
EPSILON          = 0.05          # confirmed value from last successful test
EPSILON_BOOST    = 0.0145
KAPPA            = 0.0119
A_SHELL          = 0.048         # strengthened
Y_SHELL          = np.array([-3.58, 3.58])
W_SHELL          = 0.44          # narrowed
HURST            = 0.718
RHO              = 0.305
SIGMA_THERMAL    = 0.205
BOUNDARY_A       = 1.6
TUNNEL_SCALE     = 0.012
JUMP_PROB        = 0.022
PRESCISSION_FRAC = 0.075
COLLECTIVE_KE    = 7.3

ISOTOPES = {
    'U-233':  {'mass_sym': 117.5, 'k_mass': 4.92,  'compound_A': 234, 'Z': 92},
    'U-235':  {'mass_sym': 118.4, 'k_mass': 4.50,  'compound_A': 236, 'Z': 92},
    'Pu-239': {'mass_sym': 118.8, 'k_mass': 5.30,  'compound_A': 240, 'Z': 94},
    'Cf-252': {'mass_sym': 126.0, 'k_mass': 5.30,  'compound_A': 252, 'Z': 98},
}

def generate_fgn(n_steps, hurst, n_paths=1, L=None):
    if L is None:
        indices = np.arange(n_steps)
        ii, jj = np.meshgrid(indices, indices)
        d = ii - jj
        cov = 0.5 * (np.abs(d + 1)**(2*hurst) + np.abs(d - 1)**(2*hurst) - 2*np.abs(d)**(2*hurst))
        L = np.linalg.cholesky(cov + 1e-8 * np.eye(n_steps))
    z = np.random.randn(n_paths, n_steps)
    increments = z @ L.T
    return increments, L

def simulate_one_fission(epsilon, n_steps, Z, A, dw_x, dw_y, dw_z, dw_w, dw_v):
    x = np.zeros(n_steps)
    y = np.zeros(n_steps)
    z = np.zeros(n_steps)
    w = np.zeros(n_steps)
    v = np.zeros(n_steps)

    scission_step = None

    for t in range(1, n_steps):
        x[t] = x[t-1] + SIGMA_THERMAL * dw_x[t-1]

        barrier_force = -KAPPA * y[t-1]
        shell_force = 0.0
        for ys in Y_SHELL:
            dist = y[t-1] - ys
            shell_force -= A_SHELL * dist * np.exp(-dist**2 / (2 * W_SHELL**2))

        if abs(y[t-1]) > 4.2:
            shell_force -= 0.013 * (abs(y[t-1]) - 4.2)

        epsilon_local = epsilon
        if np.any(np.abs(y[t-1] - Y_SHELL) < W_SHELL):
            epsilon_local += EPSILON_BOOST

        jump = 0.0
        if np.any(np.abs(y[t-1] - Y_SHELL) < 0.5):
            if np.random.rand() < JUMP_PROB:
                jump = np.sign(np.random.randn()) * 0.3

        y[t] = y[t-1] + epsilon_local + SIGMA_THERMAL * dw_y[t-1] + barrier_force + shell_force + jump

        z[t] = z[t-1] + SIGMA_THERMAL * dw_z[t-1]
        w[t] = w[t-1] + SIGMA_THERMAL * dw_w[t-1]
        v[t] = v[t-1] + SIGMA_THERMAL * dw_v[t-1]

        if x[t] > 0.9 * BOUNDARY_A:
            if np.random.rand() < TUNNEL_SCALE * (x[t] / BOUNDARY_A):
                scission_step = t
                break

        if x[t] >= BOUNDARY_A:
            scission_step = t
            break

    if scission_step is None:
        scission_step = n_steps - 1

    y_final = y[scission_step]
    is_asymmetric = abs(y_final) > 0.5

    tke_coulomb = 0.1189 * Z**2 / A**(1/3) * (1 + 0.1 * abs(y_final))
    pre_scission = PRESCISSION_FRAC * tke_coulomb + 5 * (scission_step / n_steps)
    tke = tke_coulomb + pre_scission + COLLECTIVE_KE - 3 * z[scission_step]**2

    return {
        'is_asymmetric': is_asymmetric,
        'y_final': y_final,
        'tke': tke,
        'scission_step': scission_step
    }

def run_simulation(isotope, n_events=8000, n_steps=1000, epsilon=EPSILON, plot=True):
    params = ISOTOPES[isotope]
    mass_sym = params['mass_sym']
    k_mass   = params['k_mass']
    A_total  = params['compound_A']
    Z        = params['Z']

    dw_x_all = np.random.randn(n_events, n_steps)
    dw_y_all, L_fgn = generate_fgn(n_steps, HURST, n_events)
    dw_z_all = np.random.randn(n_events, n_steps)
    dw_w_all = np.random.randn(n_events, n_steps)
    dw_v_all = np.random.randn(n_events, n_steps) + RHO * dw_y_all

    results = []
    for i in range(n_events):
        res = simulate_one_fission(epsilon, n_steps, Z, A_total,
                                   dw_x_all[i], dw_y_all[i],
                                   dw_z_all[i], dw_w_all[i], dw_v_all[i])
        results.append(res)

    asym_count = sum(r['is_asymmetric'] for r in results)
    asym_frac  = asym_count / n_events * 100

    masses_light = []
    for r in results:
        yf = abs(r['y_final'])
        mass_approx = mass_sym - k_mass * yf
        # Proton parity staggering (approximate Z_light)
        z_light_approx = round(Z * (mass_approx / A_total))
        parity_mod = +0.55 if z_light_approx % 2 == 0 else -0.70
        mass = mass_approx + parity_mod + np.random.normal(0, 0.7)
        masses_light.append(mass)

    # Histogram
    hist, bin_edges = np.histogram(masses_light, bins=80, range=(80, 160), density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_width   = bin_edges[1] - bin_edges[0]

    peaks, props = find_peaks(hist, height=0.02)
    light_peak_A = np.nan
    light_peak_yield = np.nan
    if len(peaks) >= 1:
        idx = peaks[np.argmin(bin_centers[peaks])]
        light_peak_A = bin_centers[idx]
        light_peak_yield = props['peak_heights'][np.where(peaks == idx)[0][0]] * bin_width * 100

    valley_idx   = np.argmin(hist[20:60]) + 20
    valley_yield = hist[valley_idx] * bin_width * 100
    pv_ratio     = light_peak_yield / valley_yield if valley_yield > 0 else np.inf

    asym_tke = np.mean([r['tke'] for r in results if r['is_asymmetric']]) if asym_count > 0 else np.nan

    # ─── Plot ───────────────────────────────────────────────────
    if plot:
        plt.figure(figsize=(10, 6))
        plt.bar(bin_centers, hist, width=bin_width, alpha=0.6, color='skyblue', edgecolor='navy')
        plt.plot(bin_centers, hist, color='darkblue', lw=1.2)

        if not np.isnan(light_peak_A):
            plt.axvline(light_peak_A, color='red', ls='--', lw=2,
                        label=f'Light peak ≈ {light_peak_A:.2f} u')

        plt.title(f'{isotope} — Light Fragment Mass Yield\n'
                  f'Asym: {asym_frac:.2f}%, P/V: {pv_ratio:.1f}, TKE: {asym_tke:.2f} MeV')
        plt.xlabel('Light fragment mass (u)')
        plt.ylabel('Normalized yield')
        plt.xlim(80, 160)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'mass_yield_{isotope}.png', dpi=150)
        # plt.show()   # uncomment only if you want interactive display

    return {
        'asymmetric_fraction': asym_frac,
        'light_peak_A': light_peak_A,
        'light_peak_yield_pct': light_peak_yield,
        'valley_yield_pct': valley_yield,
        'peak_to_valley_ratio': pv_ratio,
        'avg_TKE_asym_MeV': asym_tke
    }

# ────────────────────────────────────────────────
# Main block — run ALL isotopes
# ────────────────────────────────────────────────
if __name__ == "__main__":
    np.random.seed(42)

    print("Starting full test run of all isotopes")
    print(f"ε = {EPSILON}, n_events = 8000, n_steps = 1000")
    print("This may take 10–60 minutes depending on your machine\n")

    start_time = time.time()

    for iso in ISOTOPES:
        print(f"\n{'═'*60}")
        print(f"Running {iso}")
        print(f"{'═'*60}")
        res = run_simulation(iso, n_events=8000, n_steps=1000)
        for k, v in res.items():
            print(f"{k:22}: {v:.3f}" if isinstance(v, float) else f"{k:22}: {v}")
        print(f"Plot saved: mass_yield_{iso}.png\n")

    elapsed = time.time() - start_time
    print(f"All tests completed in {elapsed/60:.1f} minutes.")
