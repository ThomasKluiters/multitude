import pandas as pd

from data.model import GeneExpressionSeries
from scipy import stats


def t_test(gene_expression_series: GeneExpressionSeries, **kwargs):
    phenotypes = list(set(sample.phenotype for sample in gene_expression_series.samples.values()))

    if len(phenotypes) != 2:
        raise ValueError("Data must only have two phenotypes!")

    wt = next(
        phenotype for phenotype
        in phenotypes
        if ("wt" in phenotype.lower()) or ("ctrl" in phenotype.lower())
    )
    ko = next(phenotype for phenotype in phenotypes if phenotype != wt)

    a = gene_expression_series.group_phenotype(wt)
    b = gene_expression_series.group_phenotype(ko)

    (statistic, p_values) = stats.ttest_ind(a, b, axis=1, **kwargs)

    df = pd.DataFrame({"statistic": statistic, "p": p_values}, index=a.index).sort_values("statistic")
    return (df - df.min()) / (df.max() - df.min())
