import asyncio

from fractal.client.fractal import cli


def run():
    asyncio.run(cli())
