#!/bin/bash

poetry run uvicorn fractal.server.main:app --workers 2
