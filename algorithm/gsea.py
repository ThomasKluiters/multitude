from dataclasses import dataclass

import numpy as np
import pandas as pd

from algorithm.base import NonTBAlgorithm
from data.model import GeneSet, GeneExpressionSeries


@dataclass
class GSEAAlgorithm(NonTBAlgorithm):
    query_gene_set: GeneSet

    permutations: int = 1000

    def compute_distribution(self, gene_expression_series: GeneExpressionSeries):
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

            scores[:, i] = (a_mean - b_mean) / np.sqrt((a_std ** 2) / a_samples) + ((b_std ** 2) / b_samples)

            if i != (self.permutations - 1):
                np.random.shuffle(ndarray)
                pass
        return pd.Series(scores[:, 1], index=gene_expression_data.index), pd.DataFrame(scores[:, 1:], index=gene_expression_data.index)

    def find_differently_expressed(self, gene_expression_series: GeneExpressionSeries):
        rankings, permutations = self.compute_distribution(gene_expression_series)
        query = list(map(str, self.query_gene_set))
        # rankings = t_test(gene_expression_series)
        #
        # query_genes = [gene.name for gene in self.query_gene_set]
        # series = rankings["statistic"].loc[query_genes].rename("Hit")
        # annotated_rankings = rankings.join(series).fillna(-1 / len(self.query_gene_set)).sort_values('statistic',
        #                                                                                              ascending=False)
        #
        # return annotated_rankings["Hit"].cumsum()
