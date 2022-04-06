from __future__ import annotations

from dataclasses import dataclass
import logging
from pathlib import Path
import os
from datetime import datetime
from abc import ABC
from telnetlib import GA
from typing import Any, Callable

from packaging.version import parse, LegacyVersion
import docker
import docker.errors
from docker.types import Mount
import requests
import json

from bioTea.utils.errors import ImageNotFoundError
from bioTea.utils.path_checker import is_path_exists_or_creatable_portable
from bioTea.utils.tools import ConsoleWindow

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


# Some checks
def na_or(check) -> Callable:
    def _wrapped_check(argument):
        if type(argument) == str and argument.upper() in ("NA"):
            return True
        else:
            return check(argument)

    return _wrapped_check


def is_(argtype) -> Callable:
    def _wrapped_check(argument):
        return issubclass(type(argument), argtype)

    return _wrapped_check


def is_in(arglist) -> Callable:
    def _wrapped_check(argument):
        return argument in arglist

    return _wrapped_check


def is_valid_design_string(argument):
    # TODO: implement this check
    return is_(str)(argument)


def is_valid_color(argument):
    # TODO: implement this
    return is_(str)(argument)


def is_list_of(check) -> Callable:
    def _wrapped_check(argument):
        if not type(argument) == list:
            return False

        return all([check(x) for x in argument])

    return _wrapped_check


class GattacaArgument:
    """Class used to mark an argument in a GattacaInterface dict."""

    def __init__(self, check: Callable, default: Any) -> None:
        self.check = check
        self.default = default

        assert self.check(
            default
        ), "Invalid default GATTACA value. Someone coded it wrong."

    def __call__(self, argument=None):
        if argument == None:
            argument == self.default

        if not self.check(argument):
            raise ValueError(f"Argument check failed. Invalid argument {argument}")


class RequiredGattacaArgument(GattacaArgument):
    def __init__(self, check: Callable) -> None:
        self.check = check
        # The default does not matter. It HAS to be overridden.
        self.default = None

    def __call__(self, argument):
        if not self.check(argument):
            raise ValueError(f"Argument check failed. Invalid argument {argument}")


class GattacaInterface(ABC):
    """Abstract class that models an interface with GATTACA.

    Defines the arguments that can be passed to a GATTACA command, and handles
    parsing them to an object that can be given to `run_gattaca` to run the
    concrete command.
    """

    possible_args: dict = None
    """The possible arguments to the interface.

    This is a dictionary of "value_name" = GattacaArgument.
    If a Callable, it is used to test the arg before passing it to GATTACA.
    If a RequiredArgument,
    """

    # I am not 100% sure this is the correct way to use this, but it fails
    # if "possible_args" is not defined, so I'm happy.
    def __init_subclass__(cls, /, possible_args, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.possible_args = possible_args

    @classmethod
    def parse_arguments(self, **kwargs) -> str:
        """Check the input args and parse them to a GATTACA-compliant string."""
        required_args = [
            key
            for key, val in self.possible_args.items()
            if type(val) is RequiredGattacaArgument
        ]

        assert all(
            [x in required_args for x in kwargs.keys()]
        ), "Missing required args: {}".format(
            ", ".join([x for x in kwargs.keys() if x not in required_args])
        )
        assert all(
            [x in self.possible_args.keys() for x in kwargs.keys()]
        ), "Unrecognized argument(s) {}".format(
            ", ".join([x for x in kwargs.keys() if x not in self.possible_args.keys()])
        )

        for key, value in kwargs.items():
            # These will raise an error if the check fails. So, if this passes,
            # all passed args are OK.
            self.possible_args[key](value)

        # Here, we are sure of three things:
        # 1. All required arguments are overridden by `kwargs`
        # 2. All arguments in `kwargs` are in the possible_arguments.
        # 3. All arguments in `kwargs` are valid (they pass the checks)
        # Therefore, we can update the default values with the kwargs safely.
        defaults = {
            key: self.possible_args[key].default for key in self.possible_args.keys()
        }
        defaults.update(kwargs)

        # Update the keywords for python to R.
        def transmogrify_bool(value):
            if type(value) == bool:
                if value:
                    return "TRUE"
                return "FALSE"
            return value

        def transmogrify_none(value):
            if value is None:
                return "NA"
            return value

        def transmogrify(value):
            value = transmogrify_bool(value)
            value = transmogrify_none(value)
            return value

        defaults = {k: transmogrify(v) for k, v in defaults.items()}

        # GATTACA wants a space-separated, quoted list of string arguments.
        return " ".join([f'"{key}={val}"' for key, val in defaults.items()])


## >> Interface definitions
# To define a new interface, or check that one is compilant, look in the
# `entrypoint.R` of the GATTACA module, at the `defaults` list.
# Convert this list to a dictionary, replacing every `NULL` with RequiredArgument


class PrepAffyInterface(GattacaInterface):
    possible_args: dict = {
        "output.file": RequiredGattacaArgument(is_path_exists_or_creatable_portable),
        "remove.controls": GattacaArgument(is_(bool), True),
        "n_plots": GattacaArgument(is_(int), 1e10),
        # Plot options
        "use_pdf": GattacaArgument(is_(bool), True),
        "plot_width": GattacaArgument(is_(int), 16),
        "plot_height": GattacaArgument(is_(int), 9),
        "png_ppi": GattacaArgument(is_(int), 250),
        "enumerate_plots": GattacaArgument(is_(bool), True),
    }


class PrepAgilInterface(GattacaInterface):
    possible_args: dict = {
        "output.file": RequiredGattacaArgument(is_path_exists_or_creatable_portable),
        "remove.controls": GattacaArgument(is_(bool), True),
        "n_plots": GattacaArgument(is_(int), 1e10),
        "grep_pattern": GattacaArgument(is_(str), "*.(txt|TXT)"),
        # Plot options
        "use_pdf": GattacaArgument(is_(bool), True),
        "plot_width": GattacaArgument(is_(int), 16),
        "plot_height": GattacaArgument(is_(int), 9),
        "png_ppi": GattacaArgument(is_(int), 250),
        "enumerate_plots": GattacaArgument(is_(bool), True),
    }


class AnalyzeInterface(GattacaInterface):
    possible_args: dict = {
        "input.file": RequiredGattacaArgument(is_path_exists_or_creatable_portable),
        "output.dir": RequiredGattacaArgument(is_path_exists_or_creatable_portable),
        "experimental_design": RequiredGattacaArgument(is_valid_design_string),
        "contrasts": RequiredGattacaArgument(is_valid_design_string),
        "min_log2_expression": GattacaArgument(is_(float), 4.0),
        "fc_threshold": GattacaArgument(is_(float), 0.5),
        "min_groupwise_presence": GattacaArgument(is_(float), 0.8),
        "slowmode": GattacaArgument(is_(bool), False),
        "show_data_snippets": GattacaArgument(is_(bool), True),
        "annotation_database": GattacaArgument(is_(bool), True),
        "dryrun": GattacaArgument(is_(bool), False),
        "renormalize": GattacaArgument(is_(bool), False),
        "run_limma_analysis": GattacaArgument(is_(bool), True),
        "run_rankprod_analysis": GattacaArgument(is_(bool), True),
        "batches": GattacaArgument(na_or(is_valid_design_string), "NA"),
        "extra_limma_vars": GattacaArgument(
            na_or(is_list_of(is_valid_design_string)), "NA"
        ),
        "group_colors": GattacaArgument(
            is_list_of(is_valid_color),
            [
                "cornflowerblue",
                "firebrick3",
                "olivedrab3",
                "darkgoldenrod1",
                "purple",
                "magenta3",
            ],
        ),
        # Plot options
        "use_pdf": GattacaArgument(is_(bool), True),
        "plot_width": GattacaArgument(is_(int), 16),
        "plot_height": GattacaArgument(is_(int), 9),
        "png_ppi": GattacaArgument(is_(int), 250),
        "enumerate_plots": GattacaArgument(is_(bool), True),
    }


class AnnotateInterface(GattacaInterface):
    possible_args: dict = {
        "expression_data_path": RequiredGattacaArgument(
            is_path_exists_or_creatable_portable
        ),
        "output_path": RequiredGattacaArgument(is_path_exists_or_creatable_portable),
        "database_name": GattacaArgument(is_(str)),
    }


## <<


def run_gattaca(
    command: str,
    arguments: dict[str],
    interface: GattacaInterface,
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

    log.debug(f"Parsing arguments with interface '{type(interface)}'")
    try:
        parsed_args = interface.parse_arguments(**arguments)
    except ValueError as e:
        log.exception(
            "Invalid arguments passed to interface. Please open an issue with the bioTEA logs."
        )
        return 0

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
            + parsed_args
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
