def run_testserver():
    import uvicorn

    uvicorn.run("fractal_server.main:app", reload=True)
