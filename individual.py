# individual.py

"""
Klasa opisująca pojedynczego osobnika.
Przechowuje wektor fenotypu w n-wymiarowej przestrzeni.
"""
import numpy as np

class Individual:
    def __init__(self, phenotype, sex, tail=0.0):
        self.phenotype = np.array(phenotype, dtype=float)
        self.sex = sex
        self.tail = float(np.clip(tail, 0.0, 1.0))

    def get_phenotype(self):
        return self.phenotype

    def set_phenotype(self, new_phenotype):
        self.phenotype = np.array(new_phenotype, dtype=float)

    def get_sex(self):
        return self.sex

    def set_sex(self, new_sex):
        self.sex = new_sex

    def get_tail(self):
        return self.tail

    def set_tail(self, new_tail):
        self.tail = float(np.clip(new_tail, 0.0, 1.0))


