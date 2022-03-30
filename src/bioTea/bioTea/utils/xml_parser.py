"""xml parsing to get useful objects"""
# To annotate the class with itself
from __future__ import annotations

from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import BinaryIO
import xmltodict

from bioTea.utils.errors import UnsupportedChip


@dataclass
class GeoPlatform:
    accession: str
    manufacturer: str
    title: str

    def __str__(self) -> str:
        return (
            f"GEO - Platform ({self.accession}) - {self.title} from {self.manufacturer}"
        )


@dataclass
class GeoSample:
    accession: str
    source: str
    extracted_molecule: str
    organism: str
    conditions: list[str]
    suppl_data_ftp: str

    def __str__(self) -> str:
        return f"GEO - Sample ({self.accession}): {self.__dict__}"

    @classmethod
    def sample_from_xmldict(cls: GeoSample, raw_dict: dict) -> GeoSample:
        """Generate a new instance of the class from a dict parsed from an xml"""
        # Fail if the sample is more than one channel
        if (ch_count := int(raw_dict["Channel-Count"])) != 1:
            raise UnsupportedChip(f"Too many channels to parse ({ch_count})")

        if type((char := raw_dict["Channel"]["Characteristics"])) is str:
            conditions = {"type": char}
        else:
            # This is a list of small dicts. We fuse them together.
            def dictify(dict):
                return {dict["@tag"]: dict["#text"]}
            char = [dictify(x) for x in char]
            conditions = {}
            [conditions.update(x) for x in char]

        print(conditions)

        return cls(
            accession=raw_dict["@iid"],
            organism=raw_dict["Channel"]["Organism"]["#text"],
            source=raw_dict["Channel"]["Source"],
            extracted_molecule=raw_dict["Channel"]["Molecule"],
            conditions=conditions,
            suppl_data_ftp=raw_dict["Supplementary-Data"]["#text"],
        )


@dataclass
class GeoSeries:
    accession: str
    platform: GeoPlatform
    authors: list[str]
    samples: list[GeoSample]

    def __str__(self) -> str:
        authors = ", ".join(self.authors)
        samples = ", ".join([str(x) for x in self.samples])
        return f"GEO - Series ({self.accession}) from {authors}. Platform {self.platform} with {len(self.samples)} samples: {samples}"


class XmlMinimlExtractor:
    def __init__(self, xml_bytes: BinaryIO) -> None:
        self.xdict = xmltodict.parse(xml_bytes)["MINiML"]

    def extract_authors(self) -> list[str]:
        raw_authors = self.xdict["Contributor"]
        print(raw_authors)
        authors = []
        for author in raw_authors:
            if "Person" in author.keys():
                authors.append("{} {}".format(author["Person"]["First"], author["Person"]["Last"]))

        return authors

    def extract_samples(self) -> list[GeoSample]:
        raw_samples = self.xdict["Sample"]
        sample_objects = []
        for sample in raw_samples:
            sample_obj = GeoSample.sample_from_xmldict(sample)
            sample_objects.append(sample_obj)

        return sample_objects

    def extract_platform(self) -> GeoPlatform:
        return GeoPlatform(
            accession=self.xdict["Platform"]["@iid"],
            manufacturer=self.xdict["Platform"]["Manufacturer"],
            title=self.xdict["Platform"]["Title"],
        )

    def extract_id(self) -> str:
        return self.xdict["Series"]["@iid"]


# Have xmltodict do the heavy lifting
def miniml_xml_to_series(xml_path: PathLike) -> GeoSeries:
    # Follow the paths to make the objects
    with Path(xml_path).open("rb") as bytes:
        extractor = XmlMinimlExtractor(bytes)

    project = GeoSeries(
        accession=extractor.extract_id(),
        authors=extractor.extract_authors(),
        platform=extractor.extract_platform(),
        samples=extractor.extract_samples(),
    )

    return project
