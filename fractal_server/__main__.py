import uvicorn


def run():
    uvicorn.run("fractal_server.main:app", reload=True)


if __name__ == "__main__":
    run()
