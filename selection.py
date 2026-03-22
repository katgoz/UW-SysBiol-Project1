# selection.py
import numpy as np
from strategies import SelectionStrategy


# ---------------------------------------------------------------------------
# Funkcje pomocnicze (używane też w stats.py)
# ---------------------------------------------------------------------------

def fitness_function(phenotype: np.ndarray, alpha: np.ndarray, sigma: float, tail: float = 0.0) -> float:
    """
    Gaussowska funkcja fitness z kosztem ogona:
        phi_base(p) = exp( -||p - alpha||^2 / (2 * sigma^2) )
        phi_total   = phi_base * exp(-0.7 * tail)  tu można zmienić ten koszt posiadania ogona

    :param phenotype: fenotyp osobnika
    :param alpha: optymalny fenotyp środowiska
    :param sigma: parametr siły selekcji
    :param tail: długość ogona w [0,1]
    :return: wartość fitness
    """
    diff = phenotype - alpha
    base_fitness = float(np.exp(-np.dot(diff, diff) / (2 * sigma ** 2)))
    tail_cost = np.exp(-0.7 * tail)
    return float(base_fitness * tail_cost)


def compute_fitnesses(individuals: list, alpha: np.ndarray, sigma: float) -> np.ndarray:
    """Oblicza fitness dla całej listy osobników. Zwraca tablicę numpy (N,)."""
    return np.array([
        fitness_function(ind.get_phenotype(), alpha, sigma, ind.get_tail())
        for ind in individuals
    ])


# ---------------------------------------------------------------------------
# Strategie selekcji
# ---------------------------------------------------------------------------

class ThresholdSelection(SelectionStrategy):
    """
    Selekcja progowa: eliminuje osobniki o fitness poniżej progu.
    Zwraca ocalałych – może ich być mniej niż N.
    Reprodukcja uzupełni populację do N w następnym kroku.
    """

    def __init__(self, sigma: float, threshold: float):
        self.sigma = sigma
        self.threshold = threshold

    def select(self, individuals: list, alpha: np.ndarray) -> list:
        return [ind for ind in individuals
                if fitness_function(ind.get_phenotype(), alpha, self.sigma, ind.get_tail()) >= self.threshold]


class ProportionalSelection(SelectionStrategy):
    """
    Selekcja proporcjonalna (ruletka / Wright-Fisher):
    losuje N osobników z powtórzeniami, proporcjonalnie do fitness.
    Zwraca dokładnie N osobników.
    """

    def __init__(self, sigma: float, N: int):
        self.sigma = sigma
        self.N = N

    def select(self, individuals: list, alpha: np.ndarray) -> list:
        fitnesses = compute_fitnesses(individuals, alpha, self.sigma)
        total = fitnesses.sum()
        probs = fitnesses / total if total > 0 else np.ones(len(individuals)) / len(individuals)
        chosen = np.random.choice(len(individuals), size=self.N, replace=True, p=probs)
        return [individuals[i] for i in chosen]


class TwoStageSelection(SelectionStrategy):
    """
    Dwuetapowa selekcja (domyślna – zgodna z treścią zadania):
      Etap 1 – progowy: eliminuje osobniki z fitness < threshold
      Etap 2 – proporcjonalny: spośród ocalałych losuje N osobników
                               proporcjonalnie do fitness

    Zwraca dokładnie N osobników (lub pustą listę = wymarcie w etapie 1).
    """

    def __init__(self, sigma: float, threshold: float, N: int):
        self.sigma = sigma
        self.threshold = threshold
        self.N = N

    def select(self, individuals: list, alpha: np.ndarray) -> list:
        # Etap 1: selekcja progowa
        survivors = [ind for ind in individuals
                     if fitness_function(ind.get_phenotype(), alpha, self.sigma, ind.get_tail()) >= self.threshold]
        if not survivors:
            return []

        # Etap 2: selekcja proporcjonalna – wypełnia do N
        fitnesses = compute_fitnesses(survivors, alpha, self.sigma)
        total = fitnesses.sum()
        probs = fitnesses / total if total > 0 else np.ones(len(survivors)) / len(survivors)
        chosen = np.random.choice(len(survivors), size=self.N, replace=True, p=probs)
        return [survivors[i] for i in chosen]
