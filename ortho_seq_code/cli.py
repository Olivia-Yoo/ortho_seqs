import click
from ortho_seq_code.orthogonal_polynomial import ortho_poly_command
from ortho_seq_code.gui.gui import gui_run
from ortho_seq_code.plotclass import rf1d_run

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group()
# @click.pass_context
def cli():
    pass


cli.add_command(ortho_poly_command, name="orthogonal-polynomial")
cli.add_command(gui_run, name="gui")
cli.add_command(rf1d_run, name="rf1d-viz")

if __name__ == "__main__":
    cli()
