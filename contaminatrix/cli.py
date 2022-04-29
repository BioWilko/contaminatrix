import click
from types import SimpleNamespace
import pysam
import os
import pathlib
from . import contamination_detection_funcs
from . import scheme_downloader


@click.command()
@click.option(
    "--min-coverage",
    help="Minimum coverage required at a position for evaluation",
    default=10,
    type=click.INT,
)
@click.option(
    "--min-difference",
    help="Minumum difference between read groups required to flag as sus",
    type=click.FloatRange(min=0, max=1),
    default=0.1,
)
@click.option(
    "--scheme-directory",
    help="Directory where primer scheme files can be found/downloaded to (Default: ~/artic/primer-schemes/)",
    type=click.Path(),
    default=str(pathlib.Path.home()) + "/artic/primer-schemes",
)
@click.option(
    "--detailed-report",
    default=False,
    help="Print detailed report including position details",
    is_flag=True,
)
@click.argument("scheme")
@click.argument("bamfile")
def main(*_, **args):
    args = SimpleNamespace(**args)
    if not os.path.exists(args.bamfile + ".bai"):
        pysam.index(args.bamfile)

    args.bedfile, args.reference = scheme_downloader.get_scheme(
        args.scheme, args.scheme_directory
    )

    sus_positions = contamination_detection_funcs.find_sus_positions(args)

    contamination_detection_funcs.print_positions(args, sus_positions)


if __name__ == "__main__":
    main()
