# population.py


import numpy as np
from individual import Individual


class Population:
    """
    Klasa przechowuje listę osobników (Individual)
    oraz pomaga w obsłudze różnych operacji na populacji.
    """

    def __init__(self, size, n_dim, init_scale: float = 0.1, alpha_init=None):
        """
        Inicjalizuje populację losowymi fenotypami w n-wymiarach.

        :param size: liczba osobników (N)
        :param n_dim: wymiar fenotypu (n)
        :param init_scale: odchylenie std rozkładu startowego wokół optimum
        :param alpha_init: centrum inicjalizacji
        """
        center = (
            np.array(alpha_init, dtype=float)
            if alpha_init is not None
            else np.zeros(n_dim)
        )

        self.individuals = []

        for _ in range(size):
            phenotype = np.random.normal(loc=center, scale=init_scale, size=n_dim)

            if np.random.rand() < 0.5:
                tail = 0.0
            else:
                tail = np.random.rand()

            self.individuals.append(Individual(phenotype, tail))


    def get_individuals(self):
        return self.individuals

    def set_individuals(self, new_individuals):
        self.individuals = new_individuals

    def __len__(self) -> int:
        return len(self.individuals)


if __name__ == "__main__":
    pop = Population(size=10, n_dim=3)
    for ind in pop.get_individuals():
        print(ind.get_phenotype(), ind.get_sex())
