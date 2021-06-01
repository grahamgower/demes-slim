#!/bin/sh
# convert *.yaml to *.json

for yaml in *.yaml; do
    python -c \
        "import sys, demes; \
        demes.dump(demes.load(sys.argv[1]), sys.stdout, format='json', simplified=False)" \
        $yaml \
        | python -m json.tool \
        | perl -pe 's/\bInfinity\b/"Infinity"/g' \
        > ${yaml%%.yaml}.json
done
