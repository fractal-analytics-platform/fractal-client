#!/usr/bin/env python
import sys

from fractal_server.__main__ import run

if __name__ == "__main__":

    if len(sys.argv) > 1:
        port = int(sys.argv[1])
    else:
        port = 8000

    sys.exit(run(port=8001))
