import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
import os

# --- folder na figury ---
os.makedirs("figures", exist_ok=True)

# --- Wczytanie danych ---
with open('results.pkl', 'rb') as f:
    results = pickle.load(f)

# --- Konfiguracja ---
conditions = ['no_tail_effect', 'tail_effect']
drifts = ['slow', 'fast']

# --- Pomocnicze funkcje ---
def extract_series(runs, attr):
    series_list = [getattr(r, attr) for r in runs]
    max_len = max(len(s) for s in series_list)

    data = np.full((len(series_list), max_len), np.nan)

    for i, s in enumerate(series_list):
        data[i, :len(s)] = s

    return data

def mean_std(data):
    return np.nanmean(data, axis=0), np.nanstd(data, axis=0)

def mean_last_k(runs, attr, k=50):
    vals = []
    for r in runs:
        series = getattr(r, attr)
        vals.append(np.nanmean(series[-k:]))
    return np.array(vals)

def final_values(runs, attr):
    return np.array([getattr(r, attr)[-1] for r in runs])

def extinction_info(runs):
    extinct = [r.extinct_at is not None for r in runs]
    return sum(extinct), len(extinct)

# =========================================================
# 1. FITNESS
# =========================================================
fig, axes = plt.subplots(1, 2, figsize=(12,5))

for i, condition in enumerate(conditions):
    ax = axes[i]

    for drift in drifts:
        runs = results[condition][drift]

        data = extract_series(runs, 'mean_fitnesses')
        mean, std = mean_std(data)

        ax.plot(mean, label=drift)
        ax.fill_between(range(len(mean)), mean-std, mean+std, alpha=0.2)

    ax.set_title(condition.replace('_', ' '))
    ax.set_ylim(0, 0.7)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Fitness")
    ax.legend()

plt.suptitle("Effect of tail selection on mean fitness")
plt.savefig("figures/fitness.png", dpi=300)
plt.show()

# =========================================================
# 2. DISTANCE
# =========================================================
plt.figure(figsize=(10,6))

for condition in conditions:
    for drift in drifts:
        runs = results[condition][drift]

        data = extract_series(runs, 'distances_from_optimum')
        mean, std = mean_std(data)

        plt.plot(mean, label=f"{condition}-{drift}")
        plt.fill_between(range(len(mean)), mean-std, mean+std, alpha=0.2)

plt.xlabel("Generation")
plt.ylabel("Distance from optimum")
plt.title("Lag from optimum (tail vs no-tail conditions)")
plt.legend()
plt.savefig("figures/distance.png", dpi=300)
plt.show()

# =========================================================
# 3. VARIANCE
# =========================================================
plt.figure(figsize=(10,6))

for condition in conditions:
    for drift in drifts:
        runs = results[condition][drift]

        data = extract_series(runs, 'phenotype_variances')
        mean, std = mean_std(data)

        plt.plot(mean, label=f"{condition}-{drift}")
        plt.fill_between(range(len(mean)), mean-std, mean+std, alpha=0.2)

plt.xlabel("Generation")
plt.ylabel("Phenotypic variance")
plt.title("Phenotypic diversity under tail selection regimes")
plt.legend()
plt.savefig("figures/variance.png", dpi=300)
plt.show()

# =========================================================
# 4. MALE REPRODUCTIVE SUCCESS (CDF)
# =========================================================
plt.figure(figsize=(10,6))

for condition in conditions:
    for drift in drifts:
        runs = results[condition][drift]

        all_counts = []

        for r in runs:
            for gen_counts in r.male_offspring_series:
                if len(gen_counts) > 0:
                    all_counts.extend(gen_counts)

        if len(all_counts) == 0:
            continue

        all_counts = np.array(all_counts)
        sorted_counts = np.sort(all_counts)
        cdf = np.arange(1, len(sorted_counts)+1) / len(sorted_counts)

        plt.plot(sorted_counts, cdf,
                 label=f"{condition}-{drift}")

plt.xlabel("Number of offspring per male")
plt.ylabel("Cumulative fraction of males")
plt.title("CDF of male reproductive success (all generations)")
plt.legend()
plt.grid(alpha=0.3)
plt.savefig("figures/male_reproductive_cdf.png", dpi=300)
plt.show()

# =========================================================
# 5. EXTINCTION
# =========================================================
extinction_results = {}

for condition in conditions:
    extinction_results[condition] = {}

    for drift in drifts:
        runs = results[condition][drift]
        ext, total = extinction_info(runs)

        extinction_results[condition][drift] = (ext, total)
        print(f"{condition}-{drift}: {ext}/{total} extinct ({ext/total:.2f})")

# Fisher test
for drift in drifts:
    c1_ext, c1_tot = extinction_results[conditions[0]][drift]
    c2_ext, c2_tot = extinction_results[conditions[1]][drift]

    table = [
        [c1_ext, c1_tot - c1_ext],
        [c2_ext, c2_tot - c2_ext]
    ]

    _, p = fisher_exact(table)
    print(f"\nFisher test ({drift} drift): p = {p:.4f}")

labels = []
values = []

for condition in conditions:
    for drift in drifts:
        ext, total = extinction_results[condition][drift]
        labels.append(f"{condition.replace('_',' ')}\n{drift}")
        values.append(ext / total)

plt.figure(figsize=(6,4))
plt.bar(labels, values)
plt.ylabel("Extinction rate")
plt.title("Effect of tail selection on extinction rate")
plt.savefig("figures/extinction.png", dpi=300)
plt.show()

# =========================================================
# 6. FINAL FITNESS
# =========================================================
plt.figure(figsize=(8,5))

data_to_plot = []
labels = []

for condition in conditions:
    for drift in drifts:
        runs = results[condition][drift]
        finals = final_values(runs, 'mean_fitnesses')
        data_to_plot.append(finals)
        labels.append(f"{condition.replace('_',' ')}\n{drift}")

plt.boxplot(data_to_plot, tick_labels=labels)
plt.ylabel("Final mean fitness")
plt.title("Final fitness under tail selection regimes")
plt.savefig("figures/final_fitness.png", dpi=300)
plt.show()

# =========================================================
# 7. FINAL DISTANCE
# =========================================================
plt.figure(figsize=(8,5))

data_to_plot = []
labels = []

for condition in conditions:
    for drift in drifts:
        runs = results[condition][drift]
        finals = final_values(runs, 'distances_from_optimum')
        data_to_plot.append(finals)
        labels.append(f"{condition.replace('_',' ')}\n{drift}")

plt.boxplot(data_to_plot, tick_labels=labels)
plt.ylabel("Final lag")
plt.title("Final lag from optimum")
plt.savefig("figures/final_distance.png", dpi=300)
plt.show()

# =========================================================
# 8. CORRELATION (POPRAWIONE 🔥)
# =========================================================
fig, axes = plt.subplots(1, 2, figsize=(12,5))

for i, condition in enumerate(conditions):
    ax = axes[i]

    for drift in drifts:
        runs = results[condition][drift]

        data = extract_series(runs, 'corr_survivors_series')  # ✅ ZMIANA
        mean, std = mean_std(data)

        ax.plot(mean, label=drift)
        ax.fill_between(range(len(mean)), mean-std, mean+std, alpha=0.2)

    ax.axhline(0, linestyle='--', alpha=0.5)
    ax.set_title(condition.replace('_', ' '))
    ax.set_xlabel("Generation")
    ax.set_ylabel("Correlation (survivors)")
    ax.legend()

plt.suptitle("Honesty of tail signal (AFTER SELECTION)")
plt.savefig("figures/correlation_timeseries.png", dpi=300)
plt.show()

# =========================================================
# 9. FINAL CORRELATION (POPRAWIONE 🔥)
# =========================================================
plt.figure(figsize=(8,5))

data_to_plot = []
labels = []

for condition in conditions:
    for drift in drifts:
        runs = results[condition][drift]
        finals = final_values(runs, 'corr_survivors_series')  # ✅ ZMIANA
        data_to_plot.append(finals)
        labels.append(f"{condition.replace('_',' ')}\n{drift}")

plt.boxplot(data_to_plot, tick_labels=labels)
plt.ylabel("Final correlation (survivors)")
plt.title("Final honesty of tail signal")
plt.savefig("figures/final_correlation.png", dpi=300)
plt.show()

# =========================================================
# 10. TAIL EVOLUTION
# =========================================================
fig, axes = plt.subplots(1, 2, figsize=(12,5))

for i, condition in enumerate(conditions):
    ax = axes[i]

    for drift in drifts:
        runs = results[condition][drift]

        data = extract_series(runs, 'male_mean_tails')
        mean, std = mean_std(data)

        ax.plot(mean, label=drift)
        ax.fill_between(range(len(mean)), mean-std, mean+std, alpha=0.2)

    ax.set_title(condition.replace('_', ' '))
    ax.set_xlabel("Generation")
    ax.set_ylabel("Mean male tail")
    ax.legend()

plt.suptitle("Evolution of tail under selection regimes")
plt.savefig("figures/tail_evolution.png", dpi=300)
plt.show()

# =========================================================
# 11. MEAN CORRELATION (TAIL vs NO-TAIL)
# =========================================================

def mean_corr_per_run(runs, attr, last_k=50):
    vals = []
    for r in runs:
        series = getattr(r, attr)
        if len(series) == 0:
            continue
        vals.append(np.nanmean(series[-last_k:]))
    return np.array(vals)

print("\n=== MEAN CORRELATION (last generations) ===")

summary = {}

for condition in conditions:
    all_vals = []

    for drift in drifts:
        runs = results[condition][drift]
        vals = mean_corr_per_run(runs, 'tail_base_corr_series')

        print(f"{condition} - {drift}: mean = {np.mean(vals):.3f} ± {np.std(vals):.3f}")

        all_vals.extend(vals)

    all_vals = np.array(all_vals)
    summary[condition] = all_vals

    print(f"{condition} (ALL): mean = {np.mean(all_vals):.3f} ± {np.std(all_vals):.3f}\n")

# --- wykres ---
plt.figure(figsize=(6,5))

data = [summary[c] for c in conditions]
labels = [c.replace('_',' ') for c in conditions]

plt.boxplot(data, tick_labels=labels)
plt.axhline(0, linestyle='--', alpha=0.5)

plt.ylabel("Correlation (tail vs quality)")
plt.title("Mean honesty of tail signal")

plt.savefig("figures/mean_correlation_summary.png", dpi=300)
plt.show()

# =========================================================
# 12. NUMERIC SUMMARY + STAT TEST
# =========================================================

def collect_all_corr(runs, attr, last_k=50):
    vals = []
    for r in runs:
        series = getattr(r, attr)
        if len(series) == 0:
            continue
        vals.append(np.nanmean(series[-last_k:]))
    return np.array(vals)

tail_vals = []
no_tail_vals = []

for drift in drifts:
    tail_vals.extend(
        collect_all_corr(results['tail_effect'][drift], 'tail_base_corr_series')
    )
    no_tail_vals.extend(
        collect_all_corr(results['no_tail_effect'][drift], 'tail_base_corr_series')
    )

tail_vals = np.array(tail_vals)
no_tail_vals = np.array(no_tail_vals)

print("\n=== NUMERIC COMPARISON ===")
print(f"tail_effect:     mean = {tail_vals.mean():.4f} ± {tail_vals.std():.4f}")
print(f"no_tail_effect:  mean = {no_tail_vals.mean():.4f} ± {no_tail_vals.std():.4f}")

from scipy.stats import ttest_ind

t_stat, p_val = ttest_ind(tail_vals, no_tail_vals, equal_var=False)

print("\n=== STATISTICAL TEST ===")
print(f"t = {t_stat:.4f}")
print(f"p = {p_val:.6f}")

print("\n=== NUMERIC COMPARISON ===")

print(f"tail_effect:     mean = {tail_vals.mean():.4f} ± {tail_vals.std():.4f}")
print(f"tail_effect:     VAR  = {np.var(tail_vals):.6f}")

print(f"no_tail_effect:  mean = {no_tail_vals.mean():.4f} ± {no_tail_vals.std():.4f}")
print(f"no_tail_effect:  VAR  = {np.var(no_tail_vals):.6f}")

# =========================================================
# 13. BASE FITNESS DISTRIBUTION (VIOLIN) 
# =========================================================

generations_to_check = [50, 150]

for gen in generations_to_check:
    plt.figure(figsize=(10,6))

    data_to_plot = []
    labels = []

    for condition in conditions:
        for drift in drifts:
            runs = results[condition][drift]

            all_base_fit = []

            for r in runs:
                if hasattr(r, "individual_base_fitness_series"):
                    if len(r.individual_base_fitness_series) > gen:
                        gen_data = r.individual_base_fitness_series[gen]
                        if len(gen_data) > 0:
                            all_base_fit.extend(gen_data)

            if len(all_base_fit) > 0:
                data_to_plot.append(all_base_fit)
                labels.append(f"{condition.replace('_',' ')}\n{drift}")

    if len(data_to_plot) > 0:
        plt.violinplot(data_to_plot, showmeans=True)
        plt.xticks(range(1, len(labels)+1), labels)

        plt.ylabel("Base fitness (no tail cost)")
        plt.title(f"Base fitness distribution at generation {gen}")

        plt.savefig(f"figures/violin_base_fitness_gen_{gen}.png", dpi=300)
        plt.show()
