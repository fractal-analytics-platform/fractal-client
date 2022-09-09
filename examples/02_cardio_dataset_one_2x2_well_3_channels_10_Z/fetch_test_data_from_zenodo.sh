#!/bin/bash

DOI="10.5281/zenodo.7057076"
CLEAN_DOI=${DOI/\//_}
zenodo_get $DOI -o ../images/${CLEAN_DOI}
