from __future__ import annotations

from dataclasses import dataclass
import logging
from pathlib import Path
import os
from datetime import datetime
from bioTea.utils.tools import ConsoleWindow

from packaging.version import parse, LegacyVersion
import docker
import docker.errors
from docker.types import Mount
import requests
import json

from bioTea.utils.errors import ImageNotFoundError
from bioTea.utils.path_checker import is_path_exists_or_creatable_portable

log = logging.getLogger(__name__)

COMPATIBLE_VERSIONS = ["bleeding"]
"""Compatible GATTACA versions that can be ran by bioTEA"""

REPO = "cmalabscience/gattaca"

POSSIBLE_LOG_LEVELS = ("info", "debug", "error", "warning", "disable")


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


def get_installed_versions(
    client: docker.DockerClient = docker.from_env(),
) -> list[GattacaVersion]:
    local_images = client.images.list(REPO)
    local_images = [local_image.tags[0][22:] for local_image in local_images]
    return [GattacaVersion(version) for version in local_images]


def pull_gattaca_version(
    version: GattacaVersion, client: docker.DockerClient = docker.from_env()
):
    log.debug(f"Looking to see if {version} is available...")
    if not version in get_all_versions():
        log.error(f"Version {version} not found remotely.")
        raise ImageNotFoundError(str(version))

    log.info("Pulling remote image {version}...")
    client.images.get(f"{REPO}:{version}")
    log.info("Pulled image.")


def delete_gattaca_version(
    version: GattacaVersion, client: docker.DockerClient = docker.from_env()
):
    if not version in get_installed_versions():
        log.error(f"Version {version} not found locally.")
        raise ImageNotFoundError(str(version))

    log.info(f"Removing image {version}...")
    client.images.remove(image=f"{REPO}:{version}")
    log.debug(f"Done removing image.")

    return True


def get_all_versions() -> list[GattacaVersion]:
    endpoint = (
        "https://registry.hub.docker.com/v1/repositories/cmalabscience/gattaca/tags"
    )
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
    arguments: list[str],
    input_anchor: Path,
    output_anchor: Path,
    log_anchor: Path,
    version: str = "latest",
    console_level: str = "info",
    logfile_level: str = "debug",
    log_name: str = "auto",
) -> int:
    assert logfile_level in POSSIBLE_LOG_LEVELS
    assert console_level in POSSIBLE_LOG_LEVELS

    client = docker.from_env()

    if version == "latest":
        version = get_latest_version()
        log.info(f"Latest version: {version}")

    if version not in get_installed_versions(client=client):
        pull_gattaca_version(version, client=client)

    for path in [input_anchor, output_anchor, log_anchor]:
        assert is_path_exists_or_creatable_portable(
            path
        ), f"Path {path} is inaccessible."
        assert not path.is_file(), f"Path {path} points to a file, not a folder."

        if not path.exists():
            log.debug("Making anchor point: {path}")
            os.makedirs(path)

    try:
        image = client.images.get(f"{REPO}:{version}")
    except docker.errors.ImageNotFound:
        log.error(f"Cannot find local image: {version}")
        return 0

    if log_name == "auto":
        now = datetime.today().strftime("%Y-%m-%d-%H:%M:%S")
        log_name = f"GATTACA_{now}"

    try:
        # The composed command is like such:
        # UUID GUID command logname loglevel_console loglevel_file (args)
        composed_command = (
            os.getuid()
            + os.getgid()
            + command
            + log_name
            + console_level
            + logfile_level
            + " ".join(arguments)
        )
        container = client.containers.run(
            image,
            command=composed_command,
            remove=True,
            stderr=True,
            detach=True,
            mounts=[
                Mount("/GATTACA/target", output_anchor, type="bind"),
                Mount("/GATTACA/input", input_anchor, type="bind", read_only=True),
                Mount("/GATTACA/logs", log_anchor, type="bind"),
            ],
        )
    except Exception as e:
        log.error(f"Launching container failed: {e}")
        raise e

    with ConsoleWindow(10, name="GATTACA container", line_prefix="> ") as window:
        for line in container.logs(stream=True):
            window.print(line)

    log.debug("Container exited.")
    return 1
