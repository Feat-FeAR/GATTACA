import imp
import logging
from pathlib import Path
import os
import tempfile
from typing import Tuple

from bioTea.utils.option_parser import BioTeaOptions
from bioTea.utils.tools import PathLike, download_ftp, get_minimal_from_geo, make_geo_ftp

from tqdm import tqdm

log = logging.getLogger(__name__)


TEA = """
             ;,'
     _o_    ;:;' __    _     _______________
 ,-.'---`.__ ;  / /_  (_)___/_  __/ ____/   |
((j`=====',-'  / __ \/ / __ \/ / / __/ / /| |
 `-\     /    / /_/ / / /_/ / / / /___/ ___ |
    `-=-'    /_.___/_/\____/_/ /_____/_/  |_|
"""


def stage_new_analysis(staging_path: PathLike) -> Tuple[Path, Path]:
    """Stage a new analysis folder, with an output and temp folder.

    Args:
        staging_path (PathLike): The path to the output folder

    Returns:
        Tuple[Path, Path]: A tuple with first element the (resolved) output path, and second element the path to the temporary folder.
    """
    log.info("Staging new analysis project...")
    staging_path = Path(staging_path).resolve()
    if not staging_path.exists():
        os.makedirs(staging_path)
    log.debug(f"Staged {staging_path}.")
    log.info("Making a new temporary directory...")
    temp_dir = Path(tempfile.mkdtemp()).resolve()
    log.debug(f"Staged tempdir {temp_dir}")

    return (staging_path, temp_dir)


def retrieve_geo_data(output_folder: PathLike, geo_id: str):
    """Launch a complete bioTEA analysis given some options.

    Args:
        options (BioTeaOptions): _description_
    """

    staging_path, temp_dir = stage_new_analysis(output_folder)
    geo_series = get_minimal_from_geo(geo_id, temp_dir)

    log.info("Retrieving sample raw data...")
    downloaded = download_ftp(make_geo_ftp(geo_id, "suppl"), staging_path / "raw_data")

    log.info("Sanity Check: Testing congruency with MINiML file...")
    if not downloaded.keys() in [sample.suppl_data_ftp.split("\\")[-1] for sample in geo_series.samples]:
        print("PANIC")

    print("Pour done.")
