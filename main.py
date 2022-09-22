import random
import matplotlib

from algorithm.gsea import GSEAAlgorithm
from algorithm.significance import t_test
from data.loader import GeneExpressionDataSetLoader, LocalCachingGSEDownloader
from data.model import Gene
import matplotlib.pyplot as plt

if __name__ == '__main__':
    for dataset in [
        # "GSE22873",
        # "GSE6030",
        # ("GSE29048", "Gene Symbol", "genotype"),
        # "GSE70302",
        # "GSE70302",
        # "GSE58120",
        # "GSE46211",
        # "GSE49166",
        # "GSE50933",
        # "GSE62999",
        # ("GSE57917", "Gene Symbol", "onecut genotype"),
        ("GSE68571", "Gene Symbol", "disease_state")
    ]:
        loader = GeneExpressionDataSetLoader(LocalCachingGSEDownloader())
        dataset = loader.load_dataset(*dataset)
        query = {
            Gene(gene) for gene in
            (dataset.samples.get(next(iter(dataset.samples.keys()))).expression_levels.index.values) if
            ("hsp" in gene.lower())
        }

        de = GSEAAlgorithm(query_gene_set=query).find_differently_expressed(dataset)
