import dataclasses
import math
import itertools
import string
import sys
import textwrap

import demes


_slim_template = """\
initialize() {
    defineConstant("verbosity", 1);
    defineConstant("Q", $scaling_factor); // scaling factor
    defineConstant("burn_in", $burn_in); // multiples of N generations
    defineConstant("mutation_rate", $mutation_rate);
    defineConstant("recombination_rate", $recombination_rate);
    defineConstant("contig_length", $contig_length);
    defineConstant("trees_file", "$trees_file");
    initializeTreeSeq();
    initializeMutationRate(Q * mutation_rate);
    initializeMutationType("m1", 0.5, "f", 0.0);
    initializeGenomicElementType("g1", m1, 1.0);
    initializeGenomicElement(g1, 0, contig_length - 1);
    initializeRecombinationRate((1-(1-2*recombination_rate)^Q)/2);
}

function (void)dbg(object$$ current_sb) {
    if (verbosity > 0 & sim.generation == current_sb.start) {
        cat(current_sb.start);
        if (current_sb.start != current_sb.end) {
            cat(":" + current_sb.end);
        }
        catn(" " + current_sb.type + "() " + current_sb.source);
    }
}

function (integer)iround(numeric num) {
    return asInteger(round(num));
}

// Reschedule script blocks according to the burn in and scaling factor.
function (void)reschedule(numeric$$ N, object$$ current_sb) {
    g_start = sim.generation + max(1, iround(N * burn_in));
    for (sb in sim.scriptBlocks) {
        if (sb.type != "initialize" & sb != current_sb) {
            start = g_start + iround((sb.start - $t_offset) / Q);
            end = g_start + iround((sb.end - $t_offset) / Q);
            sim.rescheduleScriptBlock(sb, start, end);
        }
    }
}
"""


def first_event_time(graph: demes.Graph) -> float:
    """Return the earliest finite time in the graph."""
    times = set(deme.start_time for deme in graph.demes)
    times.update(deme.epochs[0].end_time for deme in graph.demes)
    times.update(migration.start_time for migration in graph.migrations)
    times.update(migration.end_time for migration in graph.migrations)
    times.update(pulse.time for pulse in graph.pulses)
    times.discard(math.inf)
    return max(times)


@dataclasses.dataclass
class ScriptBlock:
    start: float
    end: float = -1
    type: str = ""
    lines: list[str] = dataclasses.field(default_factory=list)

    @property
    def source(self):
        return "\n".join(self.lines)

    # += syntax
    def __iadd__(self, other):
        assert isinstance(other, str)
        self.lines.append(other)
        return self


def demes_to_slim(
    graph: demes.Graph,
    trees_file,
    mutation_rate=0,
    recombination_rate=0,
    burn_in=10,
    contig_length=100,
    scaling_factor=1.0,
    insert_debug=True,
) -> str:
    """
    Build a SLiM script for a demes graph.
    """
    graph = graph.in_generations()
    deme_id = {deme.name: j for j, deme in enumerate(graph.demes)}

    script_blocks = []

    # Initial populations.
    sb = ScriptBlock(-1)
    roots = [deme for deme in graph.demes if len(deme.ancestors) == 0]
    for root in roots:
        j = deme_id[root.name]
        first_epoch = root.epochs[0]
        sb += f"sim.addSubpop({j}, iround({first_epoch.start_size} / Q));"
        if first_epoch.selfing_rate > 0:
            sb += f"p{j}.setSelfingRate({first_epoch.selfing_rate});"
        if first_epoch.cloning_rate > 0:
            sb += f"p{j}.setCloningRate({first_epoch.cloning_rate});"

    # Initial migrations.
    for migration in graph.migrations:
        if math.isinf(migration.start_time):
            j = deme_id[migration.dest]
            k = deme_id[migration.source]
            sb += f"p{j}.setMigrationRates({k}, Q * {migration.rate});"

    # Reschedule.
    N_root = sum(root.epochs[0].start_size for root in roots)
    sb += f"reschedule({N_root} / Q, self);"
    script_blocks.append(sb)

    for deme in graph.demes:
        j = deme_id[deme.name]
        if not math.isinf(deme.start_time):
            size = round(deme.epochs[0].start_size)
            if len(deme.ancestors) == 1:
                k = deme_id[deme.ancestors[0]]
                sb = ScriptBlock(deme.start_time)
                sb += f"sim.addSubpopSplit({j}, iround({size} / Q), {k});"
                script_blocks.append(sb)
            else:
                assert len(deme.ancestors) > 1
                ids = tuple(deme_id[anc] for anc in deme.ancestors)
                ratios = tuple(deme.proportions)
                sb = ScriptBlock(deme.start_time, type="early")
                sb += f"sim.addSubpop({j}, iround({size} / Q));"
                sb += f"p{j}.setMigrationRates(c{ids}, c{ratios});"
                script_blocks.append(sb)
                # Turn migrations off again.
                zeros = tuple([0] * len(ratios))
                sb = ScriptBlock(deme.start_time, type="late")
                sb += f"p{j}.setMigrationRates(c{ids}, c{zeros});"
                script_blocks.append(sb)

        for i, epoch in enumerate(deme.epochs):
            if math.isinf(epoch.start_time):
                # Root epochs already accounted for.
                continue
            if epoch.size_function == "constant":
                if i > 0 and epoch.start_size != deme.epochs[i - 1].end_size:
                    sb = ScriptBlock(epoch.start_time)
                    sb += f"p{j}.setSubpopulationSize(iround({epoch.start_size} / Q));"
                    script_blocks.append(sb)
            elif epoch.size_function == "exponential":
                sb = ScriptBlock(epoch.start_time, epoch.end_time)
                sb += f"// deme {deme.name}: exponential growth"
                sb += f"start_size = {epoch.start_size};"
                sb += f"end_size = {epoch.end_size};"
                sb += "r = log(end_size / start_size);"
                sb += "gx = (sim.generation - self.start) / (self.end - self.start);"
                sb += "size = iround(start_size * exp(r * gx) / Q);"
                sb += f"p{j}.setSubpopulationSize(size);"
                script_blocks.append(sb)
            else:
                raise ValueError(f"unexpected size_function: {epoch.size_function}")

            if (i == 0 and epoch.selfing_rate > 0) or (
                i > 0 and epoch.selfing_rate != deme.epochs[i - 1].selfing_rate
            ):
                sb = ScriptBlock(epoch.start_time)
                sb += f"p{j}.setSelfingRate({epoch.selfing_rate});"
                script_blocks.append(sb)
            if (i == 0 and epoch.cloning_rate > 0) or (
                i > 0 and epoch.cloning_rate != deme.epochs[i - 1].cloning_rate
            ):
                sb = ScriptBlock(epoch.start_time)
                sb += f"p{j}.setCloningRate({epoch.cloning_rate});"
                script_blocks.append(sb)

    for migration in graph.migrations:
        j = deme_id[migration.dest]
        k = deme_id[migration.source]
        if not math.isinf(migration.start_time):
            sb = ScriptBlock(migration.start_time)
            sb += f"p{j}.setMigrationRates({k}, Q * {migration.rate});"
            script_blocks.append(sb)
        if migration.end_time != 0:
            sb = ScriptBlock(migration.end_time)
            sb += f"p{j}.setMigrationRates({k}, 0);"
            script_blocks.append(sb)

    for pulse in graph.pulses:
        j = deme_id[pulse.dest]
        k = deme_id[pulse.source]
        sb = ScriptBlock(pulse.time, type="early")
        sb += f"p{j}.setMigrationRates({k}, {pulse.proportion});"
        script_blocks.append(sb)
        sb = ScriptBlock(pulse.time, type="late")
        sb += f"p{j}.setMigrationRates({k}, 0);"
        script_blocks.append(sb)

    for deme in graph.demes:
        j = deme_id[deme.name]
        last_epoch = deme.epochs[-1]
        if last_epoch.end_time != 0:
            # population goes extinct
            sb = ScriptBlock(last_epoch.end_time)
            sb += f"p{j}.setSubpopulationSize(0);"
            script_blocks.append(sb)

    # end simulation
    sb = ScriptBlock(0, type="late")
    sb += "sim.treeSeqOutput(trees_file);"
    sb += "sim.simulationFinished();"
    script_blocks.append(sb)

    script = ""
    if graph.description is not None:
        wrapper = textwrap.TextWrapper(initial_indent="// ", subsequent_indent="// ")
        script += wrapper.fill(graph.description)
        script += "\n"
    if len(graph.doi) > 0:
        if script:
            script += "//\n"
        for doi in graph.doi:
            script += f"// {doi}\n"

    if script:
        script += "\n"

    # Time at which the simulation begins (generations ago).
    t_start = first_event_time(graph)
    # Add an arbitrary offset to indicate burn in.
    # The script blocks will be rescheduled anyway.
    t_offset = max(2, int(burn_in * N_root))

    script += string.Template(_slim_template).substitute(
        burn_in=burn_in,
        mutation_rate=mutation_rate,
        recombination_rate=recombination_rate,
        contig_length=contig_length,
        trees_file=trees_file,
        scaling_factor=scaling_factor,
        t_offset=t_offset,
    )

    def keyfunc(sb):
        if sb.start < 0:
            start = 1
        else:
            start = round(t_start + t_offset - sb.start)
        if sb.end >= 0:
            end = round(t_start + t_offset - sb.end)
        else:
            end = sb.end
        return (start, end, sb.type)

    script_blocks.sort(key=keyfunc)

    for (start, end, type_), group in itertools.groupby(script_blocks, keyfunc):
        sb_group = list(group)
        if start == 1:
            assert len(sb_group) == 1
        s = f"\n{start}"
        if end > 0:
            s += f":{end}"
        if type_:
            s += f" {type_}()"
        s += " {\n"
        if insert_debug:
            s += "    dbg(self);\n"
        s += textwrap.indent("\n".join(sb.source for sb in sb_group), 4 * " ")
        s += "\n}\n"
        script += s

    return script


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"usage: {sys.argv[0]} model.yaml output.trees")
        exit(1)

    yaml_file = sys.argv[1]
    trees_file = sys.argv[2]
    graph = demes.load(yaml_file)
    s = demes_to_slim(graph, trees_file=trees_file)
    print(s)
