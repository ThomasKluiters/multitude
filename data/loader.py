import os
from dataclasses import dataclass
from typing import Dict

import requests
from lxml import etree

from data.model import SeriesMetaData, SampleMetaData, PlatformData, SampleData, GeneExpressionSeries, \
    GeneExpressionData

CACHE_DIR = f"{os.path.expanduser('~')}/.multitude"
GSE_BASE_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"


class GSEDownloader:
    def do_get(self, parameters: Dict[str, str]) -> bytes:
        response = requests.get(GSE_BASE_URL, {**parameters, "form": "xml"})
        response.raise_for_status()
        return response.content

    def download_sample(self, gsm: str) -> SampleMetaData:
        data = self.do_get({"acc": gsm})
        xml = etree.fromstring(data)
        return SampleMetaData.from_miniml_xml(xml)

    def download_series(self, gse: str) -> SeriesMetaData:
        data = self.do_get({"acc": gse})
        xml = etree.fromstring(data)

        samples = tuple(map(self.download_sample, [el.get("iid") for el in xml.xpath("*[local-name() = 'Sample']")]))

        return SeriesMetaData(name=gse, samples=samples)

    def download_platform_data(self, gpl: str) -> PlatformData:
        data = self.do_get({"acc": gpl, "view": "data"})
        xml = etree.fromstring(data, parser=etree.XMLParser(huge_tree=True))
        return PlatformData.from_miniml_xml(xml.xpath("*[local-name() = 'Platform']")[0])

    def download_sample_data(self, sample_metadata: SampleMetaData) -> SampleData:
        data = self.do_get({"acc": sample_metadata.id, "view": "data"})
        xml = etree.fromstring(data, parser=etree.XMLParser(huge_tree=True))
        return SampleData.from_miniml_xml(xml.xpath("*[local-name() = 'Sample']")[0], sample_metadata)


class LocalCachingGSEDownloader(GSEDownloader):
    def do_get(self, parameters: Dict[str, str]) -> bytes:
        file_name = f"{parameters.get('acc') + parameters.get('view', 'brief')}"

        if not os.path.isdir(CACHE_DIR):
            os.mkdir(CACHE_DIR)

        cache_file_name = f"{CACHE_DIR}/{file_name}.xml"

        if not os.path.isfile(cache_file_name):
            data = super(LocalCachingGSEDownloader, self).do_get(parameters)
            with open(cache_file_name, "wb+") as fp:
                fp.write(data)
        else:
            with open(cache_file_name, "rb") as fp:
                data = fp.read()

        return data


@dataclass
class GeneExpressionDataSetLoader:
    gse_downloader: GSEDownloader

    def load_dataset(
            self,
            gse_reference: str,
            gene_symbol_column: str,
            geno_type_characteristic: str
    ) -> GeneExpressionSeries:
        series = self.gse_downloader.download_series(gse_reference)
        sample_data = list(map(self.gse_downloader.download_sample_data, series.samples))
        platform_ids = set(sample.platform.reference for sample in series.samples)
        platforms = {
            platform_id: self.gse_downloader.download_platform_data(platform_id)
            .internal_data
            .to_dataframe(["ID", gene_symbol_column])
            .rename(columns={gene_symbol_column: "Gene"})
            .set_index("ID")
            for platform_id
            in platform_ids
        }

        gene_expression_samples = {
            sample.metadata.id: GeneExpressionData(sample.internal_data.to_dataframe()
                                                   .join(platforms[sample.metadata.platform.reference], on="ID_REF")
                                                   .dropna()
                                                   .set_index("Gene")["VALUE"]
                                                   .dropna()
                                                   .rename(sample.metadata.id),
                                                   sample.metadata.characteristics.get(geno_type_characteristic))
            for sample
            in sample_data
        }

        return GeneExpressionSeries(gene_expression_samples)
