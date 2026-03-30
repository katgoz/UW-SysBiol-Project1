import pickle
import numpy as np
import matplotlib.pyplot as plt

# --- load ---
with open("results_grid.pkl", "rb") as f:
    results = pickle.load(f)

conditions = ['no_tail_effect', 'tail_effect']
c_values = sorted(results['tail_effect'].keys())

# -----------------------------
# helpers
# -----------------------------
def final_fitness(runs):
    return np.array([r.mean_fitnesses[-1] for r in runs])

def extinction_rate(runs):
    extinct = [r.extinct_at is not None for r in runs]
    return np.mean(extinct)

# =========================================================
# 1. FITNESS vs c
# =========================================================
plt.figure(figsize=(7,5))

for condition in conditions:
    means = []
    stds = []

    for c in c_values:
        runs = results[condition][c]   # ✅ prosto

        vals = final_fitness(runs)
        means.append(vals.mean())
        stds.append(vals.std())

    plt.errorbar(c_values, means, yerr=stds, marker='o',
                 label=condition.replace('_',' '))

plt.xlabel("Drift speed (c)")
plt.ylabel("Final mean fitness")
plt.title("Fitness vs environmental drift")
plt.legend()
plt.grid(alpha=0.3)

plt.savefig("figures/fitness_vs_c.png", dpi=300)
plt.show()


# =========================================================
# 2. EXTINCTION vs c
# =========================================================
plt.figure(figsize=(7,5))

for condition in conditions:
    rates = []

    for c in c_values:
        runs = results[condition][c]   # ✅ prosto

        rates.append(extinction_rate(runs))

    plt.plot(c_values, rates, marker='o',
             label=condition.replace('_',' '))

plt.xlabel("Drift speed (c)")
plt.ylabel("Extinction rate")
plt.title("Extinction vs environmental drift")
plt.legend()
plt.grid(alpha=0.3)

plt.savefig("figures/extinction_vs_c.png", dpi=300)
plt.show()