from sys import argv

import uvicorn


def run(port=8000):
    try:
        port = int(argv[1])
    except IndexError:
        pass
    uvicorn.run("fractal_server.main:app", port=port, reload=True)


if __name__ == "__main__":
    run()
