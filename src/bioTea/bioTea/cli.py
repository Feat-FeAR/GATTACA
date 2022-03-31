from genericpath import samefile
import logging
from pathlib import Path
from bioTea.docker_wrapper import get_all_versions, get_installed_versions

import typer
from tqdm import tqdm
from colorama import Fore

from bioTea import __version__
from bioTea.pour import retrieve_geo_data
from bioTea.wizard import wizard

log = logging.getLogger(__name__)

# CLI structure
# biotea
#   - info: Get version information and more.
#       - biotea: Get info about the biotea tool, its version (and DOI?)
#       - containers: Get available containers, as well as those locally installed.
#   - update: Check for container and tool updates.
#   - wizard: run the wizard
#   - retrieve: Retrieve (and format) data from GEO
#   - prepare
#       - affymetrix: Prep affymetrix data for analysis
#       - agilent: Prep agilent data for analysis
#   - analyze: Analyze with GATTACA an expression file
#   - annotate: Annotate an expression matrix
#       - generate: Generate annotations for some organism.

cli_root = typer.Typer(no_args_is_help=True)
info = typer.Typer()
prepare = typer.Typer()
annotate = typer.Typer()

cli_root.add_typer(info, name = "info")
cli_root.add_typer(prepare, name = "prepare")
cli_root.add_typer(annotate, name = "annotate")

@cli_root.callback()
def context_handler():
    log.debug(f"Starting bioTEA.")

@info.callback(invoke_without_command=True)
def generic_info(ctx: typer.Context):
    """Get information on the status of the tool."""
    if ctx.invoked_subcommand :
        return
    pass

@info.command(name = "containers")
def info_containers():
    """Get information on the downloaded and available GATTACA containers."""
    local_versions = get_installed_versions()
    remote_versions = get_all_versions()

    c = lambda x: Fore.GREEN + str(x) + Fore.RESET
    col_remote_vers = [
        c(ver) if ver in local_versions else str(ver) for ver in remote_versions 
    ]

    local_versions = [str(x) for x in local_versions]
    typer.echo(Fore.LIGHTBLUE_EX + "--- Container Info ---" + Fore.RESET)
    typer.echo("Locally installed: {}".format(", ".join(sorted(local_versions))))
    typer.echo("Remotely available: {}".format(", ".join(sorted(col_remote_vers))))
    typer.echo(
        Fore.LIGHTGREEN_EX + "Note: " + Fore.RESET +
        "Remote containers installed locally are highlighted in green."
    )
    typer.echo(Fore.LIGHTBLUE_EX + "----------------------" + Fore.RESET)

@info.command(name = "biotea")
def info_biotea():
    """Get information on the version of bioTEA."""
    typer.echo(Fore.LIGHTBLUE_EX + "--- BioTEA Info ---" + Fore.RESET)
    typer.echo(Fore.LIGHTGREEN_EX + "Version: " + Fore.RESET + __version__)

@cli_root.command(name = "update")
def update_tool():
    """Check bioTEA and the container repos for updates.
    
    This command also updates the latest container, if needed.
    """
    pass

@cli_root.command(name = "wizard")
def run_wizard():
    """Run the bioTEA wizard.

    The wizard helps in setting up, running, and exploring a GATTACA analysis.
    """
    pass

@cli_root.command(name = "retrieve")
def retrieve(output_path: Path, geo_id: str):
    """Retrieve data from GEO regarding a GEO series.

    Also helps setting the options for the GATTACA analysis.
    """
    geo_series = retrieve_geo_data(output_folder = output_path, geo_id = geo_id)
    log.info("Writing metadata...")
    with (output_path / "metadata.csv").open("w+") as fileout:
        geo_series.generate_metadata().to_csv(fileout, doublequote=True, index=False)
    
    log.info(f"Done retrieving data for {geo_id}.")

@prepare.command(name = "agilent")
def prepare_agilent():
    """Prepare agilent expression data for analysis."""
    pass

@prepare.command(name = "affymetrix")
def prepare_affymetrix():
    """Prepare affymetrix expression data for analysis."""
    pass

@cli_root.command(name = "analyze")
def run_gattaca_analysis():
    """Run Differential Gene Expression with GATTACA."""
    pass

@annotate.callback(invoke_without_command=True)
def annotate_callback(ctx: typer.Context):
    """Annotate some expression data or DEA output with annotation data."""
    if ctx.invoked_subcommand :
        return
    pass

@annotate.command()
def generate_annotations():
    """Generate annotations to use with GATTACA."""
    pass

