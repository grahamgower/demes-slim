# About

Demes-SLiM is SLiM/Eidos code for loading a
[Demes](https://popsim-consortium.github.io/demes-spec-docs/)
demographic model into
[SLiM](https://messerlab.org/slim/).

Note: a development version of SLiM is currently required (the minumum
version will be SLiM v3.7, once that is released).

# WARNING

The Demes spec has evolved since this code was written and the behaviour
of Demes-SLiM does not currently match the Demes spec. Please don't use
this code just yet!

# Usage

The `demes.slim` and `demes.eidos` files contain SLiM/Eidos functions for
parsing a Demes JSON file, and programmatically registering the model's events.
See `examples/convert.sh` for how to turn a Demes YAML file into a
fully-resolved Demes JSON file.

`demes.slim` does not form a complete SLiM script (it lacks an
`initialize()` block). Instead, the code is intended to be included in a
SLiM script by copying the `demes.slim` and `demes.eidos` files into the same
folder as the user's SLiM script and including `source("demes.slim")` in the
`initialize()` block (`demes.slim` already sources `demes.eidos`).
A complete example is provided in the `initialize.slim` file,
which can be used from the command line as follows.

```
$ slim -d JSON_FILE=\"examples/browning_america.json\" initialize.slim
```

## API

The Demes-SLiM functionality is split into two main parts: an Eidos function
`demes_load()` that loads the JSON file into a dictionary; and a SLiM function
`demes_schedule_events()` that schedules the SLiM events to create
subpopulations, change their sizes, and set migrations, etc.

The dictionary returned by `demes_load()` contains complete information about
the demographic model. So for example, a user wanting to draw a mutation in a
specific deme, could query the model dictionary for that deme's start time and
reschedule their script block(s) accordingly. The `get_deme_id()` function is
provided to obtain the integer SLiM ID of a deme from the deme's name.


---
### `function (object<Dictionary>$)demes_load(string$ json_file, [numeric$ scaling_factor = 1.0],	[numeric$ burn_in = 10.0])`

  Load a fully-resolved JSON-format Demes model from a file.
  If the file cannot be found, or is otherwise invalid, the function does not
  return (it calls `stop()`).
  The returned dictionary contains nested fields following
  [the Demes data model](https://popsim-consortium.github.io/demes-spec-docs/),
  with some changes:

  * Deme sizes are divided by the scaling factor and rounded to integers.
  * Times are converted to integer SLiM generations, which count forwards from
    the first generation (in Demes, time values increase towards the past).
  * The model dictionary is given an "end_time" field, defined as the generation
    in which the simulation should end (time=0 in most Demes models).
    The simulation then spans the open-closed interval (INF, model.end_time],
    where INF is approximated using a burn-in phase.
  * Each deme dictionary is given an "end_time" field, defined as the generation
    in which the deme goes extinct. Each deme then spans the open-closed
    interval (deme.start_time, deme.end_time]. Thus, early() events referencing
    the deme are valid for generations [deme.start_time + 1, deme.end_time].
  * Each epoch dictionary is given a "start_time" field. Each epoch then spans the
    open-closed interval (epoch.start_time, epoch.end_time].
    Parents in the epoch.start_time generation have properties following the
    previous (older) epoch (or their ancestor deme(s) if this is the first
    epoch), and offspring generated in the epoch.start_time generation have
    properties of the current epoch.

  `json_file` --- The filename of the fully-resolved JSON-format Demes
    file that is to be loaded. Such a file can be generated from a Demes
    YAML file using [demes-python](https://popsim-consortium.github.io/demes-docs/):
        `demes parse --json model.yaml > mode.json`.

  `scaling_factor` --- Scale the model by the given value. Model times
    and deme sizes will be divided by this value, and migration rates
    will be multiplied. It is the user's responsibility to multiply
    the mutation rate, recombination rate, and selection coefficients
    by the scaling factor (where relevant).

  `burn_in` --- The amount of time (in units of N generations) that will
    be used at the start of the simulation to generate the initial variation
    and genealogy in the root subpopulation(s).

---
### `function (Ni$)get_deme_id(object<Dictionary>$ model, string$ deme_name)`

  Return the integer deme ID for the named deme. This zero-based ID
  corresponds to the index of the deme in the model's list of demes.
  After `demes_schedule_events()`, this ID will also be used as the SLiM
  subpopulation ID for the deme (i.e. the index into the sim.subpopulations
  vector). Returns `NULL` if the deme is not found in the model.

  `model` --- The Demes model, as returned by `demes_load()`.

  `deme_name` --- The name of the deme being queried.

---
### `function (void)demes_schedule_events(object<Dictionary>$ model, [integer$ verbosity=1])`

  Construct and register the events that characterise the Demes model.
  This includes adding and splitting subpopulations, setting selfing
  and cloning rates, changing subpopulation sizes, and setting migration
  rates. This function should be called in SLiM generation 1.
  
  Note: It is the user's responsibility to handle sampling and simulation
  output (e.g. `sim.treeSeqOutput()`) in the appropriate generation(s),
  (typically a `late()` event at the model's end_time).

  `model` --- The Demes model, as returned by `demes_load()`.

  `verbosity` --- If set to 1 or higher, each scheduled event will be printed
    as the simulation executes.
