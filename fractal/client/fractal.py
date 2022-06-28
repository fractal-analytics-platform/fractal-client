import asyncio

from .client import cli

if __name__ == "__main__":
    asyncio.run(cli())
    exit(0)
