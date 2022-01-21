import argparse
import logging

log = logging.getLogger(__name__)

def main(args=None):

    parser = argparse.ArgumentParser(prog='chartools', description='Command line tool for manipulating charseq data')
    parser.add_argument(
        '--loglevel', default='info', help='Log level',
        choices=['debug', 'info', 'warning', 'error', 'critical'],
    )
    subparsers = parser.add_subparsers(help='Sub-commands')

    from .pairup import add_subcommand_pairup
    add_subcommand_pairup(subparsers)

    from .select import add_subcommand_select
    add_subcommand_select(subparsers)

    from .tobed import add_subcommand_tobed
    add_subcommand_tobed(subparsers)

    from .rnastats import add_subcommand_rnastats, add_subcommand_salmonstats, add_subcommand_rnastats2, add_subcommand_rnastats3
    add_subcommand_rnastats(subparsers)
    add_subcommand_rnastats2(subparsers)
    add_subcommand_rnastats3(subparsers)
    add_subcommand_salmonstats(subparsers)

    from .bam2single import add_subcommand_bam2single
    add_subcommand_bam2single(subparsers)

    from .cistrans import add_subcommand_cistrans, add_subcommand_deambig, add_subcommand_withinbp 
    add_subcommand_cistrans(subparsers)
    add_subcommand_withinbp(subparsers)
    add_subcommand_deambig(subparsers)

    from .pairup_simple import add_subcommand_pairup_simple
    add_subcommand_pairup_simple(subparsers)

    # Parse all command line arguments
    args = parser.parse_args(args)

    # This is not a good way to handle the cases
    # where help should be printed.
    # TODO: there must be a better way?
    if hasattr(args, 'func'):
        # Call the desired subcommand function
        logging.basicConfig(level=args.loglevel.upper())
        args.func(args)
        return 0
    else:
        parser.print_help()
        return 0
