from os import PathLike
from pathlib import Path
from typing import NewType
import yaml

BioTeaOptions = NewType("BioTeaOptions", dict)


def parse_options(file: PathLike) -> BioTeaOptions:
    """Parse an options file to a BioTeaOptions object.

    Runs basic checks on all the options.

    Args:
        file (PathLike): The path to the options file.

    Returns:
        BioTeaOptions: The generated object.
    """
    file = Path(file)

    # TODO :: Add checks for the user input.
    return BioTeaOptions(yaml.safe_load(file))
