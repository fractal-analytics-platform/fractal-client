class BaseInterface:
    def __init__(self, retcode: int):
        self.retcode = retcode

    def show(self, *args, **kwargs):
        raise NotImplementedError("Implement in subclasses")

    def __repr__(self):
        return f"<Interface (retcode={self.retcode})>"


class PrintInterface(BaseInterface):
    """
    Basic interface that simply prints a string to screen
    """

    def __init__(self, retcode: int, output: str):
        super().__init__(retcode)
        self.output = output

    def show(self, *args, **kwargs):
        print(self.output)
