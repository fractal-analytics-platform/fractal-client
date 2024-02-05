import json


def pretty_print(data) -> None:
    if isinstance(data, dict):
        print(json.dumps(data, indent=4))
    elif isinstance(data, list):
        for datum in data:
            pretty_print(datum)
    else:
        print(str(data))


class Interface:
    def __init__(self, retcode: int, data=None, extra=None):
        self.retcode = retcode
        self.data = data
        self.extra = extra

    def show(self):
        pretty_print(self.data)
        if self.extra is not None:
            pretty_print(self.extra)

    def __repr__(self):
        return f"<Interface (retcode={self.retcode})>"
