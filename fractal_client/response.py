import json
import logging
import sys
from json.decoder import JSONDecodeError
from typing import Any

from httpx import Response

LONG_PAYLOAD_VALUE_LIMIT = 300


def check_response(
    res: Response,
    expected_status_code: int | list[int] = 200,
    redact_long_payload: bool = False,
) -> list[Any] | dict[str, Any] | str | int | float | bool:
    """
    Check the validity of the http response from fractal server

    If the status code of the response is not one of the expected values, print
    the error to stderr and terminate with exit status 1.
    Some errors (422 errors due for failed validation of request body) are
    handled in a specific way, to make their error message more readable.

    Args:
        res: Response from `fractal-server`.
        expected_status_code: Expected status code(s).
        redact_long_payload: If `True`, redact payload values of more than
            `LONG_PAYLOAD_VALUE_LIMIT` characters.

    Returns:
        The output of `res.json()`.
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

    logging.debug(f"\nResponse status code:\n    {res.status_code}")
    if res.status_code not in expected_status_codes:

        logging.error(f"Server returned {res.status_code}")
        logging.error(
            f"Original request: {res._request.method} {res._request.url}"
        )

        payload = res._request._content.decode("utf-8")
        if redact_long_payload and len(payload) > 0:
            payload_dict = json.loads(payload)
            for key, value in payload_dict.items():
                if len(str(value)) > LONG_PAYLOAD_VALUE_LIMIT:
                    payload_dict[key] = "[value too long - redacted]"
            payload = json.dumps(payload_dict)

        logging.error(f"Original payload: {payload}")

        error_msg = data

        # Detect whether the error is due to failed request-body validation,
        # and make the error message more readable
        if (
            res.status_code == 422
            and isinstance(data, dict)
            and list(data.keys()) == ["detail"]
            and isinstance(data["detail"], list)
            and len(data["detail"]) == 1
            and isinstance(data["detail"][0], dict)
            and set(data["detail"][0].keys()) == {"msg", "type", "loc"}
        ):
            msg = data["detail"][0]["msg"]
            _type = data["detail"][0]["type"]
            loc = data["detail"][0]["loc"]
            error_msg = f"\n\tmsg: {msg}\n\ttype: {_type}\n\tloc: {loc}"

        logging.error(f"Server error message: {error_msg}\n")
        logging.error("Terminating.\n")
        sys.exit(1)

    return data
