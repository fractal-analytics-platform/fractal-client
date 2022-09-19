import logging
from json.decoder import JSONDecodeError
from sys import exit


def check_response(res, expected_status_code=200, coerce=False):
    """
    Check the validity of the http response from fractal server

    If the status code of the response is not the one expected, print the
    error to stderr and terminate with exit status 1.

    On success, optionally coerce to a pydantic model
    """
    logging.debug(res.status_code)
    try:
        data = res.json()
    except JSONDecodeError:
        data = {}
    if res.status_code != expected_status_code:
        logging.error(f"Server returned {res.status_code}")
        logging.error(
            f"Original request: {res._request.method} {res._request.url}"
        )
        logging.error(
            f"Original payload: {res._request._content.decode('utf-8')}"
        )
        logging.error(f"Server error message: {data}\n")
        logging.error("Terminating.\n")
        exit(1)

    if coerce:
        return coerce(**data)
    else:
        return data