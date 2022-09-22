from dataclasses import dataclass
from typing import Set, Tuple

import pandas as pd


@dataclass
class Gene:
    name: str


@dataclass
class Location:
    database: str
    reference: str

    @classmethod
    def from_miniml_xml(cls, element):
        return cls(element.get("database"), element.text)


@dataclass
class PlatformMetaData:
    id: str


@dataclass
class SampleMetaData:
    id: str
    location: Location
    platform: Location

    @classmethod
    def from_miniml_xml(cls, root_element):
        sample_element = root_element.xpath("*[local-name() = 'Sample']")[0]
        platform_element = root_element.xpath("*[local-name() = 'Platform']")[0]
        return cls(
            sample_element.get("iid"),
            Location.from_miniml_xml(sample_element.xpath("*[local-name() = 'Accession']")[0]),
            Location.from_miniml_xml(platform_element.xpath("*[local-name() = 'Accession']")[0])
        )


@dataclass
class SampleData:
    metadata: SampleMetaData
    platform: Location


class PlatformData(pd.DataFrame):
    pass


@dataclass
class SampleData:
    metadata: SampleMetaData


@dataclass
class SeriesMetaData:
    name: str
    samples: Tuple[SampleMetaData, ...]


GeneSet = Set[Gene]
