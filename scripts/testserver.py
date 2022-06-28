def run_testserver():
    import uvicorn

    uvicorn.run("fractal.server.main:app", reload=True)
