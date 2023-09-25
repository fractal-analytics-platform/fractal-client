import logging
from json.decoder import JSONDecodeError
from sys import exit
from typing import Union

from httpx import Response


def check_response(
    res: Response,
    expected_status_code: Union[int, list[int]] = 200,
):
    """
    Check the validity of the http response from fractal server

    If the status code of the response is not one of the expected values, print
    the error to stderr and terminate with exit status 1.

    Args:
        res: Response from `fractal-server`.
        expected_status_code: Expected status code(s).
    """

    try:
        data = res.json()
    except JSONDecodeError:
        data = {}

    # Also allow a list of expected status codes
    if isinstance(expected_status_code, list):
        expected_status_codes = expected_status_code
    else:
        expected_status_codes = [expected_status_code]

    logging.debug(res.status_code)
    if res.status_code not in expected_status_codes:

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

    return data
