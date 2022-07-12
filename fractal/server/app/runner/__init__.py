async def submit_workflow(msg):
    from devtools import debug

    debug(f"received {msg}")

    import parsl
    from parsl.config import Config

    parsl.load(Config())

    @parsl.python_app()
    def dummy():
        print("doing stuff in parsl's python app")

    dummy().result()
