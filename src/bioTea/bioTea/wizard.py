from pathlib import Path
import os

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
