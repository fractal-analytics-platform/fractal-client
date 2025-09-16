import pytest
from httpx import Request
from httpx import Response

from fractal_client.response import check_response


def test_check_response(caplog):
    JSON = {"something": "else"}
    response = Response(status_code=200, json=JSON)
    checked_response = check_response(response)
    assert checked_response == JSON

    REQUEST_BODY = {"cache_dir": "xxxxxx"}
    RESPONSE_BODY = {
        "detail": [
            {
                "loc": [
                    "body",
                    "cache_dir",
                ],
                "msg": (
                    "String attribute 'cache_dir' must be an absolute path"
                    " (given 'xxxxxx')."
                ),
                "type": "value_error",
            },
        ],
    }
    response = Response(
        status_code=422,
        json=RESPONSE_BODY,
        request=Request("GET", "http://example.org", json=REQUEST_BODY),
    )
    caplog.clear()
    with pytest.raises(SystemExit):
        check_response(response)
    assert "Original request" in caplog.records[-4].getMessage()
    assert "Original payload" in caplog.records[-3].getMessage()
    print(caplog.records[-2].getMessage())
    assert "msg: String attribute " in caplog.records[-2].getMessage()
    assert "type: value_error" in caplog.records[-2].getMessage()
    assert "loc: ['body', 'cache_dir']" in caplog.records[-2].getMessage()
    assert "Terminating" in caplog.records[-1].getMessage()

    # Test accepted status codes
    RESPONSE_BODY = {"some": "response"}
    response = Response(
        status_code=123,
        json=RESPONSE_BODY,
        request=Request("GET", "http://example.org", json={"some": "request"}),
    )
    out = check_response(response, expected_status_code=123)
    assert out == RESPONSE_BODY
    out = check_response(response, expected_status_code=[123])
    assert out == RESPONSE_BODY
