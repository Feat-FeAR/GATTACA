from itertools import combinations
from pathlib import Path
import os
from typing import Optional, Union
from itertools import permutations

import pandas as pd
import typer
from colorama import Fore

from bioTea.utils.path_checker import is_path_exists_or_creatable_portable
from bioTea import pour

TEA = """                                                    __/\__
             ;,'                               . _  \\\\''//
     _o_    ;:;' __    _     _______________   -( )-/_||_\\
 ,-.'---`.__ ;  / /_  (_)___/_  __/ ____/   |   .'. \_()_/
((j`=====',-'  / __ \/ / __ \/ / / __/ / /| |    |   | . \\
 `-\     /    / /_/ / / /_/ / / / /___/ ___ |    ϕ---| .  \\
    `-=-'    /_.___/_/\____/_/ /_____/_/  |_|   .'. ,\_____'.
                              W I Z A R D
"""

import logging

from bioTea.utils.tools import user_input
from bioTea.utils.errors import ErrorManager, BioTeaError, RaiserHandler

log = logging.getLogger(__name__)


def interactive_metadata_to_gattaca_options(
    metadata: Union[Path, pd.DataFrame],
    primary_var: Optional[str] = None,
    use_all_contrasts: bool = False,
):
    log.debug("Generating a GATTACA_options.yaml file.")
    typer.echo(
        Fore.LIGHTYELLOW_EX
        + "WARNING:"
        + Fore.RESET
        + "Due to limitations in the yaml module, the helpful comments in the "
        + "yaml options will be lost by generating the file in this way. Please refer to the manual directly for help "
        + "if you need to edit the file manually later."
    )

    if type(metadata) is Path:
        metadata = pd.read_csv(metadata)

    if primary_var is not None:
        if primary_var in metadata.columns:
            log.warn(
                "The specified primary variable is not valid. Falling back to prompting the user."
            )
            primary_var = None

    while primary_var is None:
        primary_var = typer.prompt(
            "Select the test variable",
            type=typer.Choice(metadata.columns),
            show_choices=True,
        )

        if not typer.confirm(f"Selection: {primary_var}. Is this correct?"):
            primary_var = None

    contrasts = None
    all_contrasts = permutations(metadata.columns, 2)
    if len(metadata.columns) == 2:
        base = typer.prompt(
            "Select the 'control' variable:",
            type=typer.Choice(metadata.columns),
            show_choices=True,
        )
        test = [x for x in metadata.columns if x != base][0]
        contrasts = [f"{test}-{base}"]

    while contrasts is None:
        typer.echo(
            f"Specify the contrasts of interest, as 'test-control test-control ...'. E.g. '{all_contrasts[0]}'"
        )
        contrasts = typer.prompt(
            "Please make a selection",
            type=typer.Choice(metadata.columns),
            show_choices=True,
            confirmation_prompt=True,
        )
        contrasts = contrasts.split(" ")

        if not all([x in all_contrasts for x in contrasts]):
            invalid_contrasts = ", ".join([x not in all_contrasts for x in contrasts])
            log.error(
                f"Invalid contrasts selected: {invalid_contrasts}. Please try again."
            )
            contrasts = None


def wizard_unhandled():
    print(TEA)

    print(
        "Welcome to the bioTEA wizard. This helper will guide you through the whole bioTEA analysis."
    )
    print("You can exit at any time by pressing CTRL+C")
    print(
        "Select a folder (that will be created if non-existant) to use as a working directory:"
    )
    staging_path = user_input(
        "> ",
        is_path_exists_or_creatable_portable,
        retry_prompt="The path seems invalid, or is non-writable. Try again.",
    )

    staging_path, temp_path = pour.stage_new_analysis(staging_path)

    print("Please specify the GEO id of the data. It should start with 'GSE':")
    geo_id = user_input(
        "> ",
    )

    return ""


def wizard(*args, **kwargs):
    wizard_manager = ErrorManager(wizard_unhandled)

    wizard_manager.add_handler(BioTeaError, RaiserHandler())

    wizard_manager.run(args, kwargs)
