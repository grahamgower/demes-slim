# About

This repository contains code for loading a
[Demes](https://popsim-consortium.github.io/demes-spec-docs/)
demographic model into
[SLiM](https://messerlab.org/slim/).

There are two conceptually related, but independent, bits of code here. 1) A Python
script that loads a Demes YAML file and prints out a SLiM script, and 2) SLiM/Eidos
code for loading a fully-resolved Demes JSON file and registering SLiM events.
The code here is experimental and has not been rigorously tested.

# Usage

## `demes2slim.py`

This script uses the
[`demes` Python library](https://popsim-consortium.github.io/demes-docs/)
to load a Demes YAML file and then prints the SLiM script for the model
to stdout. This can be piped directly into the `slim` executable.
As SLiM runs the script, the events which characterise the model are
printed for debugging purposes.

```
$ python demes2slim.py
usage: demes2slim.py model.yaml output.trees
```

```
$ python demes2slim.py examples/gutenkunst_ooa.yaml /tmp/out.trees | slim
// Initial random seed:
24302768653795

// RunInitializeCallbacks():
initializeTreeSeq();
initializeMutationRate(0);
initializeMutationType(1, 0.5, "f", 0);
initializeGenomicElementType(1, m1, 1);
initializeGenomicElement(g1, 0, 99);
initializeRecombinationRate(0);

// Starting run at generation <start>:
1 

1 early() {
    dbg(self);
    sim.addSubpop(0, iround(7300 / Q));
    reschedule(7300 / Q, self);
}
73001 early() {
    dbg(self);
    sim.addSubpopSplit(1, iround(12300 / Q), 0);
    p0.setSubpopulationSize(0);
}
76201 early() {
    dbg(self);
    sim.addSubpopSplit(2, iround(2100 / Q), 1);
    sim.addSubpopSplit(3, iround(12300 / Q), 1);
    p2.setMigrationRates(3, Q * 0.00025);
    p3.setMigrationRates(2, Q * 0.00025);
    p1.setSubpopulationSize(0);
}
80953 early() {
    dbg(self);
    sim.addSubpopSplit(4, iround(1000 / Q), 2);
    sim.addSubpopSplit(5, iround(510 / Q), 2);
    p2.setMigrationRates(3, 0);
    p3.setMigrationRates(2, 0);
    p4.setMigrationRates(3, Q * 3e-05);
    p3.setMigrationRates(4, Q * 3e-05);
    p5.setMigrationRates(3, Q * 1.9e-05);
    p3.setMigrationRates(5, Q * 1.9e-05);
    p5.setMigrationRates(4, Q * 9.6e-05);
    p4.setMigrationRates(5, Q * 9.6e-05);
    p2.setSubpopulationSize(0);
}
80953:81801 early() {
    dbg(self);
    // deme CEU: exponential growth
    start_size = 1000;
    end_size = 29725;
    r = log(end_size / start_size);
    gx = (sim.generation - self.start) / (self.end - self.start);
    size = iround(start_size * exp(r * gx) / Q);
    p4.setSubpopulationSize(size);
    // deme CHB: exponential growth
    start_size = 510;
    end_size = 54090;
    r = log(end_size / start_size);
    gx = (sim.generation - self.start) / (self.end - self.start);
    size = iround(start_size * exp(r * gx) / Q);
    p5.setSubpopulationSize(size);
}
81801 late() {
    dbg(self);
    sim.treeSeqOutput(trees_file);
    sim.simulationFinished();
}
```

```
$ tskit info /tmp/out.trees
sequence_length:  100.0
trees:            1
samples:          192230
individuals:      96115
nodes:            350769
edges:            350768
sites:            0
mutations:        0
migrations:       0
populations:      6
provenances:      1
```

## `demes.slim`

The `demes.slim` file contains SLiM/Eidos functions for parsing a Demes
JSON file, and programmatically registering the model's events.
See `examples/convert.sh` for how to turn a Demes YAML file into a
fully-resolved Demes JSON file.

`demes.slim` does not form a complete SLiM script (it lacks an
`initialize()` block). Instead, the code is intended to be included in a
SLiM script with `source("demes.slim")`. A complete example is provided in
`initialize.slim`. The `initialize.slim` script can also be used from
the command line as follows.

```
$ slim -d model=\"examples/browning_america.json\" initialize.slim
// Initial random seed:
24422252123378

// RunInitializeCallbacks():
initializeTreeSeq();
initializeMutationRate(0);
initializeMutationType(1, 0.5, "f", 0);
initializeGenomicElementType(1, m1, 1);
initializeGenomicElement(g1, 0, 99);
initializeRecombinationRate(0);

// Starting run at generation <start>:
1 

1 late() {
    // ancestral: conjured from the void
    sim.addSubpop(0, 731);
}
732 early() {
    // AMH: split from ancestral
    sim.addSubpopSplit(1, 1447, 0);
}
732 late() {
    // ancestral: extinction
    p0.setSubpopulationSize(0);
}
1120 early() {
    // OOA: split from AMH
    sim.addSubpopSplit(2, 186, 1);
    // AFR: split from AMH
    sim.addSubpopSplit(3, 1447, 1);
    // AFR -> OOA
    p2.setMigrationRates(3, 0.0015);
    // OOA -> AFR
    p3.setMigrationRates(2, 0.0015);
}
1120 late() {
    // AMH: extinction
    p1.setSubpopulationSize(0);
}
1232 early() {
    // EUR: split from OOA
    sim.addSubpopSplit(4, 100, 2);
    // EAS: split from OOA
    sim.addSubpopSplit(5, 51, 2);
    // AFR -| OOA
    p2.setMigrationRates(3, 0);
    // OOA -| AFR
    p3.setMigrationRates(2, 0);
    // AFR -> EUR
    p4.setMigrationRates(3, 0.00025);
    // EUR -> AFR
    p3.setMigrationRates(4, 0.00025);
    // AFR -> EAS
    p5.setMigrationRates(3, 7.8e-05);
    // EAS -> AFR
    p3.setMigrationRates(5, 7.8e-05);
    // EUR -> EAS
    p5.setMigrationRates(4, 0.000311);
    // EAS -> EUR
    p4.setMigrationRates(5, 0.000311);
}
1232:1324 early() {
    // EUR: exponential growth
    start_size = 100;
    end_size = 3404;
    r = log(end_size / start_size);
    gx = (sim.generation - self.start) / (self.end - self.start);
    size = iround(start_size * exp(r * gx));
    p4.setSubpopulationSize(size);
    // EAS: exponential growth
    start_size = 51;
    end_size = 4585;
    r = log(end_size / start_size);
    gx = (sim.generation - self.start) / (self.end - self.start);
    size = iround(start_size * exp(r * gx));
    p5.setSubpopulationSize(size);
}
1232 late() {
    // OOA: extinction
    p2.setSubpopulationSize(0);
}
1323 early() {
    // ADMIX: admixture of AFR EUR EAS
    sim.addSubpop(6, 3000);
    p6.setMigrationRates(c(3, 4, 5), c(0.167, 0.333, 0.5));
}
1323:1324 early() {
    // ADMIX: exponential growth
    start_size = 3000;
    end_size = 5466;
    r = log(end_size / start_size);
    gx = (sim.generation - self.start) / (self.end - self.start);
    size = iround(start_size * exp(r * gx));
    p6.setSubpopulationSize(size);
}
1323 late() {
    // ADMIX
    p6.setMigrationRates(c(3, 4, 5), c(0, 0, 0));
}
1324 late() {
    sim.treeSeqOutput(trees_file);
}
```

```
$ tskit info /tmp/slim_out.trees
sequence_length:  100.0
trees:            1
samples:          29804
individuals:      14902
nodes:            51809
edges:            51804
sites:            0
mutations:        0
migrations:       0
populations:      7
provenances:      1
```
