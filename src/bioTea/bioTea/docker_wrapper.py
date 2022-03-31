from __future__ import annotations

from dataclasses import dataclass
import logging
from pathlib import Path
from bioTea.utils.errors import ImageNotFoundError
from bioTea.utils.path_checker import is_path_exists_or_creatable_portable

from packaging.version import parse, LegacyVersion, Version

import docker
import requests
import json

log = logging.getLogger(__name__)

COMPATIBLE_VERSIONS = ["bleeding"]
"""Compatible GATTACA versions that can be ran by bioTEA"""

@dataclass
class GattacaVersion:
    raw_version: str

    @property
    def realversion(self):
        return parse(self.raw_version)
    
    def __lt__(self, other: GattacaVersion) -> bool:
        return self.realversion < other.realversion
    
    def __eq__(self, other: object) -> bool:
        return self.realversion == other.realversion
    
    def __str__(self) -> str:
        return self.raw_version


def get_installed_versions(client: docker.DockerClient = docker.from_env()) -> list[GattacaVersion]:
    local_images = client.images.list("cmalabscience/gattaca")
    local_images = [local_image.tags[0][22:] for local_image in local_images]
    return [GattacaVersion(version) for version in local_images]


def pull_gattaca_version(version: GattacaVersion, client: docker.DockerClient = docker.from_env()):
    log.debug(f"Looking to see if {version} is available...")
    if not version in get_all_versions():
        log.error(f"Version {version} not found remotely.")
        raise ImageNotFoundError(str(version))
    
    log.info("Pulling remote image {version}...")
    client.images.get(f"cmalabscience/gattaca:{version}")
    log.info("Pulled image.")


def delete_gattaca_version(version: GattacaVersion, client: docker.DockerClient = docker.from_env()):
    if not version in get_installed_versions():
        log.error(f"Version {version} not found locally.")
        raise ImageNotFoundError(str(version))

    log.info(f"Removing image {version}...")
    client.images.remove(image = f"cmalabscience/gattaca:{version}")
    log.debug(f"Done removing image.")

    return True


def get_all_versions() -> list[GattacaVersion]:
    endpoint = "https://registry.hub.docker.com/v1/repositories/cmalabscience/gattaca/tags"
    res = requests.get(endpoint, timeout=20)
    
    images = json.loads(res.text)

    versions = []
    for image in images:
        version = GattacaVersion(image["name"])
        if type(version.realversion) is LegacyVersion:
            continue
        versions.append(version)

    return versions


def get_latest_version() -> GattacaVersion:
    return sorted(get_all_versions())[-1]


def run_gattaca(
    command: str,
    arguments: str,
    input_anchor: Path,
    output_anchor: Path,
    log_anchor: Path,
    version: str = "latest",
    re_log: bool = True,
    console_level: str = "info",
    logfile_level: str = "debug"
):
    assert is_path_exists_or_creatable_portable(input_anchor)
    assert is_path_exists_or_creatable_portable(output_anchor)
    assert is_path_exists_or_creatable_portable(log_anchor)

    client = docker.from_env()

    pass

