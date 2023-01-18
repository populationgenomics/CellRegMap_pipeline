"""
just testing click
"""

import click


@click.command()
@click.option('--genes')
def main(genes: str):
    """
    does this work right?
    """
    print(genes)
    print(genes.split())


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
