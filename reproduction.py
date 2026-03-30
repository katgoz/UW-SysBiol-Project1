# reproduction.py

import copy
import numpy as np
from individual import Individual
from strategies import ReproductionStrategy
from selection import fitness_function


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
    def __init__(self, use_tail: bool = True):
        self._last_counts = {}
        self._male_counts = {}
        self._female_counts = {}
        self.use_tail = use_tail

    def acceptance_probability(self, male):
      if not self.use_tail:
          return 0.5

      t = male.get_tail()
      return 0.05 + 0.7 * t

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
        self._male_counts = {}
        self._female_counts = {}

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

                """
                w = np.random.normal(0.5, 0.25, size=len(mother_pheno))
                w = np.clip(w, 0.0, 1.0)
                child_phenotype = w * mother_pheno + (1 - w) * father_pheno
                """
                mask = np.random.rand(len(mother_pheno)) < 0.5
                child_phenotype = np.where(mask, mother_pheno, father_pheno)

                if np.random.rand() < 0.5:
                    child_tail = 0.0
                else:
                    child_tail = np.clip(
                        father.get_tail() + np.random.normal(0, 0.1),
                        1e-8, 1.0
                    )

                child = Individual(child_phenotype, child_tail)
                children.append(child)

                self._female_counts[mother_id] = self._female_counts.get(mother_id, 0) + 1
                self._male_counts[father_id] = self._male_counts.get(father_id, 0) + 1

        return children

    def get_male_offspring_counts(self):
      return list(self._male_counts.values()) if self._male_counts else [0]

    def get_reproduction_stats(self) -> dict:
      if not self._male_counts and not self._female_counts:
          return {'n_parents': 0, 'median_offspring': 0.0, 'max_offspring': 0}

      counts = list(self._male_counts.values()) + list(self._female_counts.values())
      counts = np.array(counts)

      return {
          'n_parents': int(len(counts)),
          'median_offspring': float(np.median(counts)),
          'max_offspring': int(np.max(counts)),
      }


def asexual_reproduction(survivors: list, N: int) -> list:
    return AsexualReproduction().reproduce(survivors, N)

