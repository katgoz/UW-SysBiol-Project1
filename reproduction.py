# reproduction.py

import copy
import numpy as np
import copy
import numpy as np
from individual import Individual
from strategies import ReproductionStrategy


class AsexualReproduction(ReproductionStrategy):
    def __init__(self):
        self._last_counts: np.ndarray = np.array([])

    def reproduce(self, survivors: list, target_size: int) -> list:
        if not survivors:
            self._last_counts = np.array([])
            return []
        indices = np.random.randint(0, len(survivors), size=target_size)
        self._last_counts = np.bincount(indices, minlength=len(survivors))
        return [copy.deepcopy(survivors[i]) for i in indices]

    def get_reproduction_stats(self) -> dict:
        if len(self._last_counts) == 0:
            return {'n_parents': 0, 'median_offspring': 0.0, 'max_offspring': 0}
        reproducing = self._last_counts[self._last_counts > 0]
        return {
            'n_parents': int(len(reproducing)),
            'median_offspring': float(np.median(reproducing)) if len(reproducing) else 0.0,
            'max_offspring': int(self._last_counts.max()),
        }


class SexualReproduction(ReproductionStrategy):
    def __init__(self):
        self._last_counts = {}

    def acceptance_probability(self, male):
        t = male.get_tail()
        return 0.1 + 0.8 * t

    def choose_father_for_female(self, males):
        if len(males) == 0:
            return None

        candidate_count = min(5, len(males))
        candidates = list(np.random.choice(males, size=candidate_count, replace=False))

        for male in candidates:
            u = np.random.rand()
            p = self.acceptance_probability(male)

            if u < p:
                return male

        return None

    def reproduce(self, survivors: list, target_size: int) -> list:
        if not survivors:
            self._last_counts = {}
            return []

        females = [ind for ind in survivors if ind.get_sex() == "F"]
        males = [ind for ind in survivors if ind.get_sex() == "M"]

        if len(females) == 0 or len(males) == 0:
            self._last_counts = {}
            return []

        children = []
        self._last_counts = {}

        shuffled_females = list(np.random.permutation(females))

        for mother in shuffled_females:
            if len(children) >= target_size:
                break

            father = self.choose_father_for_female(males)

            if father is None:
                continue

            mother_id = id(mother)
            father_id = id(father)

            for _ in range(3):
                if len(children) >= target_size:
                    break

                mother_pheno = mother.get_phenotype()
                father_pheno = father.get_phenotype()

                mask = np.random.rand(len(mother_pheno)) < 0.5
                child_phenotype = np.where(mask, mother_pheno, father_pheno)

                child_sex = np.random.choice(["F", "M"])

                if child_sex == "M":
                    child_tail = np.clip(
                        father.get_tail() + np.random.normal(0, 0.05),
                        0.0, 1.0
                    )
                else:
                    child_tail = 0.0

                child = Individual(child_phenotype, child_sex, child_tail)
                children.append(child)

                self._last_counts[mother_id] = self._last_counts.get(mother_id, 0) + 1
                self._last_counts[father_id] = self._last_counts.get(father_id, 0) + 1

        return children

    def get_reproduction_stats(self) -> dict:
        if not self._last_counts:
            return {'n_parents': 0, 'median_offspring': 0.0, 'max_offspring': 0}

        counts = np.array(list(self._last_counts.values()))
        return {
            'n_parents': int(len(counts)),
            'median_offspring': float(np.median(counts)),
            'max_offspring': int(np.max(counts)),
        }


def asexual_reproduction(survivors: list, N: int) -> list:
    return AsexualReproduction().reproduce(survivors, N)


