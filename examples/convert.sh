#!/bin/sh
# convert *.yaml to *.json

for yaml in *.yaml; do
    python -c \
        "import sys, demes; \
        demes.dump(demes.load(sys.argv[1]).in_generations(), \
                   sys.stdout, format='json', simplified=False)" \
        $yaml \
        | python -m json.tool \
        | sed 's/Infinity/"Infinity"/' \
        > ${yaml%%.yaml}.json
done
