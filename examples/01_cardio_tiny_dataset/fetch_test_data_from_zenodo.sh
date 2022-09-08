#!/bin/bash

DOIs="10.5281/zenodo.7059515"

for DOI in $DOIs; do
    CLEAN_DOI=${DOI/\//_}
    echo $CLEAN_DOI
    zenodo_get $DOI -o images/${CLEAN_DOI}
done
