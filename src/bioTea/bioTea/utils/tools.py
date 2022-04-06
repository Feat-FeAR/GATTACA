from __future__ import annotations
from pdb import line_prefix

import sys
import logging
from typing import Callable, Optional, TypeAlias, Union
from pathlib import Path
import ftplib
from bioTea.utils.path_checker import is_path_exists_or_creatable_portable
import requests
import zipfile
import os
import gzip
import tarfile
from collections import deque

from tqdm import tqdm
from tqdm.utils import CallbackIOWrapper

from bioTea.utils.errors import InvalidGeoId, InvalidPathError

log = logging.getLogger(__name__)

PathLike: TypeAlias = Union[Path, str]

CompressedFile = Union[tarfile.TarFile, zipfile.ZipFile, gzip.GzipFile]


class TempWorkdir:
    """Change the workdir temporarily.

    Probably super unsafe for threading.
    """

    def __init__(self, tempwd: PathLike) -> None:
        tempwd = Path(tempwd)

        self.last_wd = os.getcwd()
        self.new_wd = tempwd

    def __enter__(self):
        os.chdir(self.new_wd)

    def __exit__(self):
        os.chdir(self.last_wd)


def ask_user(
    prompt: str = "Continue?",
    options: dict = {"yes": lambda: None, "no": lambda: sys.exit()},
    allow_partial: bool = True,
    ignore_case: bool = True,
) -> None:
    """Ask the user a question with a number of options.

    Execute a fuction based on the answer given. Returns the value of the function ran.

    Args:
        prompt (str): The prompt that the user sees.
        options (dict, optional): A dict with keys as the possible answers and as values callables that will be called when the function is ran. Defaults to {"yes": lambda: None, "no": lambda: sys.exit()}.
        allow_partial (bool, optional): Allow partial answers (detected by startswith - careful with keys that start the same)?. Defaults to True.
        ignore_case (bool, optional): Ignore case of options and input?. Defaults to True.
    """
    if ignore_case:
        options = dict((k.lower(), v) for k, v in options.iteritems())

    prompt = "{} [{}]:".format(prompt, "/".join(options.keys()))

    while True:
        response = input(prompt)
        if ignore_case:
            response = response.lower()

        if allow_partial:
            for k, v in options.iteritems():
                if k.startswith(response):
                    return v()
        else:
            for k, v in options.iteritems():
                if k == response:
                    return v()

        print("Please reply with one of '{}'".format(" ,".join(options.keys())))


def user_input(
    prompt: str,
    test: Callable,
    retry: bool = True,
    retry_prompt: str = f"Invalid input. Please try again.",
) -> str:
    """Ask for input from the user.

    Tests if the answer is valid with a function, and can be made to retry until the function passes.

    Args:
        prompt (str): The prompt shown to the user at every iteration.
        test (Callable): The test applied to the input. Needs to return a bool.
        retry (bool, optional): Should the function retry? If false, raises a ValueError if the test did not pass. Defaults to True.
        retry_prompt (str, optional): Text to show before giving the prompt again to the user. Defaults to f"Invalid input. Please try again.".

    Raises:
        ValueError: If the test fails and `retry` is not true.

    Returns:
        str: The (tested) user input.
    """
    while True:
        raw_input = input(prompt)

        if test(input):
            return input

        if not retry:
            raise ValueError(
                f"Invalid check {test} for input {raw_input}, and cannot retry."
            )

        print(retry_prompt)


def download_ftp(
    ftp_url: str, destination: PathLike, filter: Callable = lambda x: True
) -> dict:
    """Downloads a file or a set of files from a ftp url

    Args:
        ftp_url: A string. May or may not start with `ftp://`
        destination: PathLike. The destination folder.
        filter: A function that takes as string. Applied to all filenames to
            download. Will dowload only those that make the function return true.

    Returns:
        A dict with keys the names of the downloaded paths and values the
        paths to the files on disk
    """
    destination = Path(destination)
    result = {}

    if not destination.exists():
        log.warn(f"Making destination: {destination}")
        os.makedirs(destination)

    if ftp_url.startswith("ftp://"):
        ftp_url = ftp_url[6:]
    host = ftp_url.split("/")[0]
    path = "/".join(ftp_url.split("/")[1:])

    log.info("Connecting to remote server....")
    log.debug(f"Connecting to '{host}'...")
    ftp = ftplib.FTP(host, timeout=240)
    log.debug(f"Logging in to '{host}'...")
    ftp.login()

    log.debug(f"Moving to '{path}'")
    ftp.cwd(path)

    files = []
    ftp.retrlines("LIST", callback=files.append)
    files = [x.split()[-1] for x in files]

    if not files:
        log.info("No files to download.")
        return {}

    all_files = len(files)
    files = [file for file in files if filter(file)]

    log.info(
        "Found {} files: {} ({} filtered out)".format(
            len(files), ", ".join(files), all_files - len(files)
        )
    )

    if not files:
        log.warning("Filter function removed all files. Nothing to do.")
        return {}

    destinations = [destination / file for file in files]

    log.debug(f"Starting download of {destinations}")
    for file, destination in tqdm(zip(files, destinations)):
        with destination.open("wb+") as outfile, tqdm(
            unit="B", unit_scale=True, unit_divisor=1024
        ) as pbar:
            outfile = CallbackIOWrapper(pbar.update, outfile, "write")

            def cb(data):
                outfile.write(data)

            ftp.retrbinary(f"RETR {file}", cb)

        result.update({file: destination})

    log.info("Finished downloading files.")
    return result


def ping(url: str) -> bool:
    """Ping some url and return true if it exist"""
    if not url.startswith("http://") or url.startswith("https://"):
        url = f"http://{url}"
    response = requests.get(url)
    return response.status_code == 200


# {obsid} is the obscured id, with the last three chars substituted by "nnn":
# e.g.: GPL1234 |> GPL1nnn
# {id} is the actual id
GEO_ENDPOINTS = {
    "GPL": {
        "miniml": "ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/{obsid}/{id}/miniml/",
        "soft": "ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/{obsid}/{id}/soft/",
        "suppl": "ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/{obsid}/{id}/suppl/",
    },
    "GDS": {
        "miniml": "ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/{obsid}/{id}/miniml/",
        "soft": "ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/{obsid}/{id}/soft",
        "suppl": "ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/{obsid}/{id}/suppl",
    },
    "GSM": {
        "miniml": "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/{obsid}/{id}/miniml/",
        "soft": "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/{obsid}/{id}/soft",
        "suppl": "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/{obsid}/{id}/suppl",
    },
    "GSE": {
        "miniml": "ftp://ftp.ncbi.nlm.nih.gov/geo/series/{obsid}/{id}/miniml/",
        "soft": "ftp://ftp.ncbi.nlm.nih.gov/geo/series/{obsid}/{id}/soft",
        "suppl": "ftp://ftp.ncbi.nlm.nih.gov/geo/series/{obsid}/{id}/suppl",
    },
}


def make_geo_ftp(geo_id: str, type: str = "miniml") -> str:
    """Get a geo FTP url from an id.

    Type can be 'miniml', 'soft' or 'suppl'

    Args:
        geo_id (str): The geo id to format with.
        type (str, optional): The type of file needed. Can be "miniml", "soft" or "suppl" for minimal files, soft files or supplementary files, respectively. Defaults to "miniml".

    Raises:
        InvalidGeoId: If the GEO id cannot be used to generate an ftp url.

    Returns:
        str: A string with the GEO ftp url pointing to the right destination.
    """
    assert type in ["miniml", "soft", "suppl"]
    geo_id = geo_id.upper()
    obsid = f"{geo_id[:-3]}nnn"

    try:
        endpoints = GEO_ENDPOINTS[geo_id[:3]]
    except KeyError:
        raise InvalidGeoId(f"Invalid ID root '{geo_id[:3]}'.")

    return endpoints[type].format(obsid=obsid, id=geo_id)


def contains_all(x: list, y: list) -> bool:
    if len(x) != len(y):
        return False

    return all(val in y for val in x)


class ConsoleWindow:
    def __init__(
        self,
        height: int,
        name: Optional[str] = None,
        line_prefix: str = "",
        clear: bool = False,
    ) -> None:
        self.height = height
        self.name = name
        self.line_prefix = line_prefix
        self.clear = clear

        self.buffer = deque([""] * self.height, maxlen=self.height)

    def print(self, line) -> None:
        self.buffer.append(line)
        self.update()

    @staticmethod
    def clearlines(n):
        # Use ANSI codes to clear `n` lines
        for _ in range(n):
            # The first goes up one line, the second clears it.
            print("\033[1A", end="\x1b[2K")

    def update(self) -> None:
        self.clearlines(self.height)
        for line in self.buffer:
            print(self.line_prefix + line)

    def __enter__(self) -> ConsoleWindow:
        if self.name:
            print(self.line_prefix + f"-- {self.name} -- " + "\n" * self.height)
        else:
            print(self.line_prefix + "\n" * self.height)

        return self

    def __exit__(self, type, value, traceback) -> None:
        self.update()

        if self.clear:
            if self.name:
                self.clearlines(self.height + 1)
            else:
                self.clearlines(self.height)


def make_path_valid(path: Path, dir: bool = False):
    """Make a path valid to write in

    If the path points to a real directory or file, makes sure that we can read
    and write there. Otherwise, makes the directory/parent directory so we
    can use the path directly.

    Args:
        path (Path): The path to test
        dir (bool, optional): Is the path theoretically pointing to  a directory?. Defaults to False.

    Raises:
        InvalidPathError: If the path is invalid.

    Returns:
        Path: The same path in input.
    """
    if is_path_exists_or_creatable_portable(path) and not path.exists():
        if dir:
            os.makedirs(path, exist_ok=True)
        else:
            os.makedirs(path.parent, exist_ok=True)

    if not is_path_exists_or_creatable_portable(path):
        raise InvalidPathError(f"Invalid path: {path}")

    return path
