import logging
from pathlib import Path

import typer
from tqdm import tqdm

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
#       - affymetrix: Prep affy data for analysis
#       - agilent: Prep agilent data for analysis
#   - analize: Analize with GATTACA an expression file
#   - annotate: Annotate an expression matrix
#       - generate: Generate annotations for some organism.

cli_root = typer.Typer(no_args_is_help=True)
info = typer.Typer()
prepare = typer.Typer()
annotate = typer.Typer()

cli_root.add_typer(info, name = "info")
cli_root.add_typer(prepare, name = "prepare")
cli_root.add_typer(annotate, name = "annotate")

@info.callback(invoke_without_command=True)
def generic_info(ctx: typer.Context):
    """Get information on the status of the tool."""
    if ctx.invoked_subcommand :
        return
    pass

@info.command(name = "containers")
def info_containers():
    """Get information on the downloaded and available GATTACA containers."""
    pass

@info.command(name = "biotea")
def info_biotea():
    """Get information on the version of bioTEA."""
    pass

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
    retrieve_geo_data(output_folder=output_path, geo_id = geo_id)

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

