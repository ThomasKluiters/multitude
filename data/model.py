from collections import Counter
from dataclasses import dataclass
from io import StringIO
from typing import Set, Tuple, List, Optional, Dict

import pandas as pd


@dataclass(unsafe_hash=True)
class Gene(str):
    name: str


@dataclass(unsafe_hash=True)
class Location:
    database: str
    reference: str

    @classmethod
    def from_miniml_xml(cls, element):
        return cls(element.get("database"), element.text)


@dataclass
class InternalDataColumn:
    name: str
    description: Optional[str]
    position: int

    @classmethod
    def from_miniml_xml(cls, root_element) -> 'InternalDataColumn':
        position = int(root_element.get("position"))
        name_element = root_element.xpath("*[local-name() = 'Name']")
        description_element = root_element.xpath("*[local-name() = 'Description']")
        return cls(name_element[0].text if name_element else None,
                   description_element[0].text if description_element else "Unknown", position)


@dataclass
class InternalData:
    columns: Tuple[InternalDataColumn, ...]
    data: str

    def ordered_columns(self) -> List[InternalDataColumn]:
        return list(sorted(self.columns, key=lambda column: column.position))

    def to_dataframe(self, columns: Optional[List[str]] = None):
        if columns is not None:
            columns = [column for column in self.ordered_columns() if column.name in columns]
        else:
            columns = self.ordered_columns()
        column_indices = [col.position - 1 for col in columns]
        return pd.read_csv(StringIO(self.data), usecols=column_indices, names=[col.name for col in columns],
                           header=None, delimiter='\t', low_memory=True)

    @classmethod
    def from_miniml_xml(cls, root_element) -> 'InternalData':
        columns = tuple(map(InternalDataColumn.from_miniml_xml, root_element.xpath("*[local-name() = 'Column']")))
        data_element = root_element.xpath("*[local-name() = 'Internal-Data']")[0]
        return cls(columns, data_element.text)


@dataclass
class PlatformMetaData:
    id: str


@dataclass
class SampleMetaData:
    id: str
    location: Location
    platform: Location
    characteristics: Dict[str, str]

    @classmethod
    def from_miniml_xml(cls, root_element):
        sample_element = root_element.xpath("*[local-name() = 'Sample']")[0]
        platform_element = root_element.xpath("*[local-name() = 'Platform']")[0]
        characteristics = dict(
            (characteristic.get("tag", "").strip(), characteristic.text.strip())
            for characteristic
            in root_element.xpath("//*[local-name() = 'Characteristics']")
        )
        return cls(
            sample_element.get("iid"),
            Location.from_miniml_xml(sample_element.xpath("*[local-name() = 'Accession']")[0]),
            Location.from_miniml_xml(platform_element.xpath("*[local-name() = 'Accession']")[0]),
            characteristics
        )


@dataclass
class SampleData:
    metadata: SampleMetaData
    platform: Location


@dataclass
class PlatformData:
    internal_data: InternalData

    @classmethod
    def from_miniml_xml(cls, root_element):
        return cls(InternalData.from_miniml_xml(root_element.xpath("*[local-name() = 'Data-Table']")[0]))


@dataclass
class SampleData:
    internal_data: Optional[InternalData]
    metadata: SampleMetaData

    @classmethod
    def from_miniml_xml(cls, root_element, metadata: SampleMetaData):
        data_table_element = root_element.xpath("*[local-name() = 'Data-Table']")
        internal_data = InternalData.from_miniml_xml(data_table_element[0]) if data_table_element else None
        return cls(internal_data, metadata)


@dataclass
class SeriesMetaData:
    name: str
    samples: Tuple[SampleMetaData, ...]


@dataclass
class GeneExpressionData:
    expression_levels: pd.Series
    phenotype: str


@dataclass
class GeneExpressionSeries:
    samples: Dict[str, GeneExpressionData]

    def matrix(self) -> pd.DataFrame:
        sorted_samples = reversed(
            sorted([sample for sample in self.samples.values()], key=lambda sample: sample.phenotype))
        return pd.concat([sample.expression_levels for sample in sorted_samples], axis=1)

    def phenotype_distribution(self):
        return Counter(sample.phenotype for sample in self.samples.values())

    def phenotypes(self) -> List[str]:
        return list(set(sample.phenotype for sample in self.samples.values()))

    def group_phenotype(self, phenotype: str) -> pd.DataFrame:
        return pd.concat([data.expression_levels for data in self.select_phenotype(phenotype).values()], axis=1)

    def select_phenotype(self, phenotype: str) -> Dict[str, GeneExpressionData]:
        return {key: value for key, value in self.samples.items() if value.phenotype == phenotype}


GeneSet = Set[Gene]
