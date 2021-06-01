#!/bin/sh

die() {
    echo $@
    exit 1
}

for file in *.eidos; do
    eidos $file || die "$file: failed"
done

for Q in 1 10.0; do
    for json in ../examples/*.json; do
        echo "Q=$Q; $json"
        states_file=`mktemp`
        slim -d Q=$Q -d json=\"$json\" test_demes_schedule_events.slim \
            > $states_file \
            || die "test_demes_schedule_events.slim: failed for $json"
        yaml=${json%%.json}.yaml
        ./verify_slim_against_demes.py -Q $Q $states_file $yaml \
            || die "verify_slim_against_demes.py: failed for $yaml"
        rm $states_file
    done
done
