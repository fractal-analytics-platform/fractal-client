import asyncio

from fractal.client.fractal import cli


def run():
    asyncio.run(cli())


if __name__ == "__main__":
    run()
