class Interface:
    def __init__(self, retcode: int, data=None, extra=None):
        self.retcode = retcode
        self.data = data
        self.extra = extra

    def show(self):
        print(str(self.data))
        if self.extra is not None:
            print(str(self.extra))

    def __repr__(self):
        return f"<Interface (retcode={self.retcode})>"
