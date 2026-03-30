import numpy as np
import pickle
import config

from main import run_simulation
from population import Population
from environment import LinearShiftEnvironment
from selection import TwoStageSelection
from mutation import IsotropicMutation
from reproduction import AsexualReproduction, SexualReproduction

results = {}



for condition, use_tail, use_tail_cost in [
        ('no_tail_effect', False, False),
        ('tail_effect', True, True)
    ]:
    results[condition] = {}

    for drift_name, c_value in [('slow', 0.005), ('fast', 0.02)]:

        print(f"Running: {condition}, drift={drift_name}")

        runs = []

        for seed in range(20):
            np.random.seed(seed)

            pop = Population(
                config.N,
                config.n,
                config.init_scale,
                alpha_init=config.alpha0
            )

            env = LinearShiftEnvironment(
                config.alpha0.copy(),
                np.ones_like(config.alpha0) * c_value,
                config.delta
            )

            sel = TwoStageSelection(
                config.sigma,
                config.threshold,
                config.N,
                use_tail_cost=use_tail_cost
            )

            mut = IsotropicMutation(
                config.mu,
                config.mu_c,
                config.xi
            )

            rep_strategy = SexualReproduction(use_tail=use_tail)

            stats = run_simulation(
                pop, env, sel, rep_strategy, mut,
                max_generations=config.max_generations,
                frames_dir=None,
                verbose=False
            )

            runs.append(stats)

        results[condition][drift_name] = runs

with open('results.pkl', 'wb') as f:
    pickle.dump(results, f)