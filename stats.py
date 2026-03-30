# stats.py

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional
import numpy as np


@dataclass
class GenerationRecord:
    generation: int
    mean_fitness: float
    mean_phenotype: np.ndarray
    phenotype_variance: float
    distance_from_optimum: float
    population_size: int
    male_mean_tail: float = 0.0
    n_parents: int = 0
    median_offspring: float = 0.0
    max_offspring: int = 0
    extra: dict = field(default_factory=dict)


class SimulationStats:

    def __init__(self) -> None:
        self.records: List[GenerationRecord] = []
        self.extinct_at: Optional[int] = None
        self.alpha_history: List[np.ndarray] = []
        self.male_offspring_distributions: List[List[int]] = []

        # przechowuje ostatnią korelację po selekcji
        self._last_corr_survivors: float = np.nan

    # =========================================================
    #   korelacja po selekcji
    # =========================================================
    def record_survivors(self, survivors, alpha, sigma):
        males = [ind for ind in survivors if ind.get_sex() == "M"]

        if len(males) < 2:
            self._last_corr_survivors = np.nan
            return

        tails = np.array([m.get_tail() for m in males])
        phenos = np.array([m.get_phenotype() for m in males])

        diffs = phenos - alpha
        base_fit = np.exp(-np.sum(diffs * diffs, axis=1) / (2 * sigma**2))

        if np.std(tails) > 1e-8:
            corr = np.corrcoef(tails, base_fit)[0, 1]
        else:
            corr = 0.0

        self._last_corr_survivors = corr

    # =========================================================
    # Główne zbieranie statystyk
    # =========================================================
    def record(self, generation: int, population, alpha: np.ndarray,
               sigma: float, reproduction_strategy=None) -> None:

        from selection import compute_fitnesses

        individuals = population.get_individuals()
        if not individuals:
            return

        self.alpha_history.append(alpha.copy())

        phenotypes = np.array([ind.get_phenotype() for ind in individuals])
        male_inds = [ind for ind in individuals if ind.get_sex() == "M"]

        male_tails = np.array([ind.get_tail() for ind in male_inds], dtype=float) if male_inds else np.array([])
        male_mean_tail = float(male_tails.mean()) if len(male_tails) > 0 else 0.0

        fitnesses = compute_fitnesses(individuals, alpha, sigma)

        mean_phenotype = phenotypes.mean(axis=0)
        phenotype_variance = phenotypes.var(axis=0).mean()
        distance = float(np.linalg.norm(mean_phenotype - alpha))
        mean_fitness = float(fitnesses.mean())

        # --- reprodukcja stats ---
        repro = (reproduction_strategy.get_reproduction_stats()
                 if reproduction_strategy is not None else None) or {}

        if reproduction_strategy is not None and hasattr(reproduction_strategy, "get_male_offspring_counts"):
            self.male_offspring_distributions.append(
                reproduction_strategy.get_male_offspring_counts()
            )
        else:
            self.male_offspring_distributions.append([])

        # =========================================================
        # korelacja w CAŁEJ populacji (słabsza, pomocnicza)
        # =========================================================
        if len(male_inds) > 1:
            male_pheno = np.array([ind.get_phenotype() for ind in male_inds])

            diffs_m = male_pheno - alpha
            base_fit_m = np.exp(-np.sum(diffs_m * diffs_m, axis=1) / (2 * sigma**2))

            if np.std(male_tails) > 1e-8:
                corr_population = np.corrcoef(male_tails, base_fit_m)[0, 1]
            else:
                corr_population = 0.0
        else:
            corr_population = np.nan

        # =========================================================
        # zapis rekordu
        # =========================================================
        self.records.append(GenerationRecord(
            generation=generation,
            mean_fitness=mean_fitness,
            mean_phenotype=mean_phenotype,
            phenotype_variance=float(phenotype_variance),
            distance_from_optimum=distance,
            population_size=len(individuals),
            male_mean_tail=male_mean_tail,
            n_parents=repro.get('n_parents', 0),
            median_offspring=repro.get('median_offspring', 0.0),
            max_offspring=repro.get('max_offspring', 0),
            extra={
                'tail_base_corr': corr_population,
                'corr_survivors': self._last_corr_survivors
            }
        ))

    # =========================================================
    # meta
    # =========================================================
    def mark_extinct(self, generation: int) -> None:
        self.extinct_at = generation

    # =========================================================
    # properties
    # =========================================================
    @property
    def male_offspring_series(self):
        return self.male_offspring_distributions

    @property
    def tail_base_corr_series(self):
        return np.array([r.extra.get('tail_base_corr', np.nan) for r in self.records])

    @property
    def corr_survivors_series(self):
        return np.array([r.extra.get('corr_survivors', np.nan) for r in self.records])

    @property
    def generations(self) -> np.ndarray:
        return np.array([r.generation for r in self.records])

    @property
    def mean_fitnesses(self) -> np.ndarray:
        return np.array([r.mean_fitness for r in self.records])

    @property
    def distances_from_optimum(self) -> np.ndarray:
        return np.array([r.distance_from_optimum for r in self.records])

    @property
    def phenotype_variances(self) -> np.ndarray:
        return np.array([r.phenotype_variance for r in self.records])

    @property
    def population_sizes(self) -> np.ndarray:
        return np.array([r.population_size for r in self.records])

    @property
    def n_parents_series(self) -> np.ndarray:
        return np.array([r.n_parents for r in self.records])

    @property
    def median_offspring_series(self) -> np.ndarray:
        return np.array([r.median_offspring for r in self.records])

    @property
    def max_offspring_series(self) -> np.ndarray:
        return np.array([r.max_offspring for r in self.records])

    @property
    def male_mean_tails(self) -> np.ndarray:
        return np.array([r.male_mean_tail for r in self.records])

    # =========================================================
    # summary
    # =========================================================
    def survived(self) -> bool:
        return self.extinct_at is None

    def final_mean_fitness(self) -> float:
        if not self.records:
            return 0.0
        return self.records[-1].mean_fitness

    def summary(self) -> str:
        if not self.records:
            return "Brak danych."

        last = self.records[-1]
        status = (f"Wymarła w pokoleniu {self.extinct_at}"
                  if self.extinct_at is not None else "Przeżyła")

        return (
            f"Pokoleń: {last.generation + 1} | Status: {status}\n"
            f"Ostatnie śr. fitness: {last.mean_fitness:.4f} | "
            f"Odległość od optimum: {last.distance_from_optimum:.4f} | "
            f"Wariancja fenotypowa: {last.phenotype_variance:.4f} | "
            f"Śr. ogon samców: {last.male_mean_tail:.4f}"
        )
