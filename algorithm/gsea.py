from dataclasses import dataclass
from typing import Tuple

import numpy as np
import pandas as pd

from algorithm.base import NonTBAlgorithm
from data.model import GeneSet, GeneExpressionSeries


@dataclass
class GSEAAlgorithm(NonTBAlgorithm):
    query_gene_set: GeneSet

    permutations: int = 10

    def compute_distribution(self, gene_expression_series: GeneExpressionSeries) -> Tuple[pd.Series, np.ndarray]:
        gene_expression_data = gene_expression_series.matrix()
        ndarray = gene_expression_data.to_numpy().T.copy()
        (samples, genes) = ndarray.shape

        a_samples = gene_expression_series.phenotype_distribution()["lung cancer"]
        b_samples = samples - a_samples

        scores = np.zeros((genes, self.permutations))

        for i in range(self.permutations):
            a = ndarray[a_samples:, :]
            b = ndarray[:a_samples, :]

            a_mean = np.mean(a, axis=0)
            b_mean = np.mean(b, axis=0)

            a_std = np.std(a, axis=0)
            b_std = np.std(b, axis=0)

            scores[:, i] = np.abs((a_mean - b_mean) / np.sqrt((a_std ** 2) / a_samples) + ((b_std ** 2) / b_samples))

            if i != (self.permutations - 1):
                np.random.shuffle(ndarray)
                pass
        return pd.Series(scores[:, 1], index=gene_expression_data.index), scores

    def find_differently_expressed(self, gene_expression_series: GeneExpressionSeries):
        rankings, permutations = self.compute_distribution(gene_expression_series)

        for _ in range(1):
            query = list(map(str, self.query_gene_set))
            gene_indices = np.array(rankings.index.isin(query), dtype=bool)
            scores = np.where(np.tile(gene_indices, (self.permutations, 1)).T, permutations, 1 / rankings.size)
