#!/bin/sh
# convert *.yaml to *.json

for yaml in *.yaml; do
    demes parse -j $yaml > ${yaml%%.yaml}.json
done
