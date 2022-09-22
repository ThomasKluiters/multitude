import os
from typing import Dict

import requests
from lxml import etree

from data.model import SeriesMetaData, SampleMetaData, PlatformData

CACHE_DIR = f"{os.path.expanduser('~')}/.multitude"
GSE_BASE_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"


class GSEDownloader:
    def do_get(self, parameters: Dict[str, str]) -> bytes:
        response = requests.get(GSE_BASE_URL, {**parameters, "form": "xml"})
        response.raise_for_status()
        return response.content

    def download_sample(self, gsm: str) -> SampleMetaData:
        data = self.do_get({"acc": gsm})
        parser = etree.XMLParser()
        xml = etree.fromstring(data, parser=parser)
        return SampleMetaData.from_miniml_xml(xml)

    def download_series(self, gse: str) -> SeriesMetaData:
        data = self.do_get({"acc": gse})
        parser = etree.XMLParser()
        xml = etree.fromstring(data, parser=parser)

        samples = tuple(map(self.download_sample, [el.get("iid") for el in xml.xpath("*[local-name() = 'Sample']")]))

        return SeriesMetaData(name=gse, samples=samples)

    def download_platform_data(self, gpl: str) -> PlatformData:
        data = self.do_get({"acc": gpl, "view": "data"})
        return data


class LocalCachingGSEDownloader(GSEDownloader):
    def do_get(self, parameters: Dict[str, str]) -> bytes:
        file_name = parameters.get("acc")

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


if __name__ == '__main__':
    for dataset in [
        "GSE22873",
        "GSE6030",
        "GSE29048",
        "GSE70302",
        "GSE70302",
        "GSE58120",
        "GSE46211",
        "GSE49166",
        "GSE50933",
        "GSE62999",
        "GSE57917",
    ]:
        loader = LocalCachingGSEDownloader()
        series = loader.download_series(dataset)

