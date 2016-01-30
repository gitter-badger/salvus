#! -*- coding: utf-8 -*-

import click

from code_generation import code_generator

@click.group()
def code_generation():
    """Test."""
    pass

@code_generation.command()
@click.option('--dimension', type=click.Choice(['2', '3']), required=True, help='Dimension of element.')
@click.option('--polynomial_order', required=True, help='Lagrange polynomial order.')
def HyperCube_generate_gll_basis(dimension, polynomial_order):
    """Autogenerate the appropriate code for a tensorized gll basis on a 2D/3D hypercube element.

    :param dimension: Dimension of element.
    :param polynomial_order: Lagrange polynomial order.
    """

    if dimension == '2':
        code_generator.tensorized_basis_2D(int(polynomial_order))
    else:
        pass


@click.group()
def solver_operation():
    """Test."""
    pass

@click.group()
def optimization_tools():
    """Test."""
    pass