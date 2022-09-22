from abc import ABC

from data.model import GeneExpressionSeries


class Algorithm:
    def find_differently_expressed(self, gene_expression_series: GeneExpressionSeries):
        raise NotImplementedError


class NonTBAlgorithm(Algorithm, ABC):
    pass
