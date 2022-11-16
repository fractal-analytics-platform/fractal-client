import argparse as ap
from sys import argv

import uvicorn

parser = ap.ArgumentParser()
parser.add_argument("--host", default="127.0.0.1")
parser.add_argument("-p", "--port", default=8000, type=int)
parser.add_argument("--reload", default=False, action="store_true")


def run():
    args = parser.parse_args(argv[1:])
    uvicorn.run(
        "fractal_server.main:app",
        host=args.host,
        port=args.port,
        reload=args.reload,
    )


if __name__ == "__main__":
    run()
