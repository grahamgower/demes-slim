#!/usr/bin/env python3
import collections
import copy
import math
import sys
import typing

import demes


class State(typing.NamedTuple):
    time: int
    pid: int
    ancestor_pid: int
    size: int
    selfing_rate: float
    cloning_rate: float
    immigrant_pids: tuple[int]
    immigrant_fracs: tuple[float]


def parse_slim_states(filename, pfx="#state"):
    """
    1: late() {
        for (p in sim.subpopulations) {
            ancestor_id = p.getValue("demes_ancestor_id");
            if (isNULL(ancestor_id)) {
                ancestor_id = -1;
            }
            catn(c("#state",
                   sim.generation,
                   p.id,
                   ancestor_id,
                   length(p.genomes),
                   p.selfingRate,
                   p.cloningRate,
                   p.immigrantSubpopIDs,
                   p.immigrantSubpopFractions));
        }
    }
    """
    states = []
    with open(filename) as f:
        for line in f:
            if not line.startswith(pfx):
                continue
            fields = line.split()[1:]
            assert len(fields) >= 6
            n = len(fields[6:]) // 2
            immigrant_pids = tuple(map(int, fields[6 : 6 + n]))
            immigrant_fracs = tuple(map(float, fields[6 + n :]))
            state = State(
                time=int(fields[0]),
                pid=int(fields[1]),
                ancestor_pid=int(fields[2]),
                size=int(fields[3]),
                selfing_rate=float(fields[4]),
                cloning_rate=float(fields[5]),
                immigrant_pids=immigrant_pids,
                immigrant_fracs=immigrant_fracs,
            )
            states.append(state)
    return states


def group_states_by_pid(states):
    gstates = collections.defaultdict(list)
    for state in states:
        gstates[state.pid].append(state)
    return gstates


def get_ancestors(states):
    ancestor_pid = states[0].ancestor_pid
    assert all(state.ancestor_pid == ancestor_pid for state in states)
    if ancestor_pid == -1:
        ancestors = list(states[0].immigrant_pids)
        proportions = list(states[0].immigrant_fracs)
    else:
        ancestors = [ancestor_pid]
        proportions = [1.0]
    return ancestors, proportions


def epoch_of_deme_at_time(deme: demes.Deme, time: float) -> int:
    """
    Return the epoch of the deme at the given time.
    """
    for epoch in deme.epochs:
        if epoch.start_time >= time >= epoch.end_time:
            break
    else:
        raise ValueError(f"deme {deme.name} doesn't exist at time {time}")

    return epoch


def round_half_up(x):
    """
    Return x rounded to the nearest integer, with halves rounded up to match
    SLiM's behaviour. Python's round() function rounds halves to the nearest
    even number.
    """
    floor_x = math.floor(x)
    if x - floor_x < 0.5:
        return floor_x
    else:
        return math.ceil(x)


# Redefine round() to have round-half-up behaviour.
round = round_half_up


def size_at_time(epoch: demes.Epoch, time: float) -> int:
    """
    Return the population size of the epoch at the given time.
    """
    assert epoch.start_time >= time >= epoch.end_time

    if math.isclose(time, epoch.end_time) or epoch.start_size == epoch.end_size:
        N = epoch.end_size
    else:
        assert epoch.size_function == "exponential"
        dt = (epoch.start_time - time) / epoch.time_span
        r = math.log(epoch.end_size / epoch.start_size)
        N = epoch.start_size * math.exp(r * dt)
    return N


def graph_to_scaled_discrete(graph: demes.Graph, scaling_factor=1.0) -> demes.Graph:
    graph = graph.in_generations()
    b = demes.Builder(description=graph.description, doi=graph.doi)
    for deme in graph.demes:
        if math.isinf(deme.start_time):
            start_time = math.inf
        else:
            start_time = round(deme.start_time / scaling_factor)
        epochs = []
        for epoch in deme.epochs:
            end_time = round(epoch.end_time / scaling_factor)
            # Require an integer number of diploid individuals.
            start_size = 2 * round(epoch.start_size / 2 / scaling_factor)
            end_size = 2 * round(epoch.end_size / 2 / scaling_factor)
            epochs.append(
                dict(
                    end_time=end_time,
                    start_size=start_size,
                    end_size=end_size,
                    selfing_rate=epoch.selfing_rate,
                    cloning_rate=epoch.cloning_rate,
                )
            )

        b.add_deme(
            deme.name,
            description=deme.description,
            start_time=start_time,
            ancestors=deme.ancestors,
            proportions=deme.proportions,
            epochs=epochs,
        )

    for migration in graph.migrations:
        b.add_migration(
            source=migration.source,
            dest=migration.dest,
            start_time=round(migration.start_time / scaling_factor),
            end_time=round(migration.end_time / scaling_factor),
            rate=scaling_factor * migration.rate,
        )

    for pulse in graph.pulses:
        b.add_pulse(
            source=pulse.source,
            dest=pulse.dest,
            time=round(pulse.time / scaling_factor),
            proportion=pulse.proportion,
        )

    return b.resolve()


# XXX: When demes exports a migration matrix function, use that instead.
def migration_matrices(graph: demes.Graph):
    """
    Return a list of migration matrices, and a list of end times that
    partition them. The start time for the first matrix is inf.
    """
    uniq_times = set(migration.start_time for migration in graph.migrations)
    uniq_times.update(migration.end_time for migration in graph.migrations)
    uniq_times.discard(math.inf)
    end_times = sorted(uniq_times, reverse=True)
    n = len(graph.demes)
    mm_list = [[[0] * n for _ in range(n)] for _ in range(len(end_times))]
    deme_id = {deme.name: j for j, deme in enumerate(graph.demes)}
    for migration in graph.migrations:
        start_time = math.inf
        for k, end_time in enumerate(end_times):
            if start_time <= migration.end_time:
                break
            if end_time < migration.start_time:
                source_id = deme_id[migration.source]
                dest_id = deme_id[migration.dest]
                if mm_list[k][dest_id][source_id] > 0:
                    raise ValueError(
                        "multiple migrations defined for "
                        f"source={migration.source}, dest={migration.dest} "
                        f"between start_time={start_time}, end_time={end_time}"
                    )
                mm_list[k][dest_id][source_id] = migration.rate
            start_time = end_time
    return mm_list, end_times


def check_states_against_graph(states, graph):
    gstates = group_states_by_pid(states)
    deme_names = [deme.name for deme in graph.demes]
    assert len(gstates) == len(deme_names), f"{len(gstates)} != {len(deme_names)}"

    def time_warp(t, model_end_time=states[-1].time):
        """Translate SLiM generation to Demes time."""
        if t <= 2:
            return math.inf
        else:
            return model_end_time - t

    for pid, local_states in gstates.items():
        name = deme_names[pid]
        t_min = min(state.time for state in local_states)
        t_max = max(state.time for state in local_states)
        start_time = time_warp(t_min)
        end_time = time_warp(t_max)
        ancestor_ids, proportions = get_ancestors(local_states)
        ancestors = [deme_names[i] for i in ancestor_ids]

        if math.isinf(start_time) and len(ancestors) > 0:
            raise AssertionError(
                f"subpop {pid} introduced at generation {t_min} has no ancestors"
            )

        deme = graph[name]
        if math.isinf(deme.start_time):
            assert math.isinf(start_time)
        else:
            assert deme.start_time == start_time, f"{deme.start_time} != {start_time}"
        assert deme.end_time == end_time, f"{deme.end_time} != {end_time}"
        assert len(deme.ancestors) == len(
            ancestors
        ), f"{len(deme.ancestors)} != {len(ancestors)}"
        assert len(deme.proportions) == len(
            proportions
        ), f"{len(deme.proportions)} != {len(proportions)}"
        assert sorted(zip(deme.ancestors, deme.proportions)) == sorted(
            zip(ancestors, proportions)
        ), f"{deme.ancestors} / {deme.proportions} != {ancestors} / {proportions}"

        for state in local_states:
            time = time_warp(state.time)
            epoch = epoch_of_deme_at_time(deme, time)
            deme_size = round(size_at_time(epoch, time))
            # XXX: 5% relative tolerance. This is quite horrible.
            assert math.isclose(deme_size, state.size, rel_tol=0.05), (
                f"{pid} ({name}): {state.time} ({time}): "
                f"{deme_size} != {state.size}"
            )
            assert math.isclose(epoch.selfing_rate, state.selfing_rate), (
                f"{pid} ({name}): {state.time} ({time}): "
                f"{epoch.selfing_rate} != {state.selfing_rate}"
            )
            assert math.isclose(epoch.cloning_rate, state.cloning_rate), (
                f"{pid} ({name}): {state.time} ({time}): "
                f"{epoch.cloning_rate} != {state.cloning_rate}"
            )

    pulses_at_time = dict()
    for pulse in graph.pulses:
        if pulse.time not in pulses_at_time:
            pulses_at_time[pulse.time] = []
        pulses_at_time[pulse.time].append(pulse)

    deme_id = {deme.name: j for j, deme in enumerate(graph.demes)}

    mm_list, end_times = migration_matrices(graph)
    start_time = math.inf
    i = 0
    for migration_matrix, end_time in zip(mm_list, end_times):
        while i < len(states):
            state = states[i]
            time = time_warp(state.time)
            if time < end_time:
                break
            i += 1
            dest_id = state.pid
            if time == start_time:
                # XXX: Migrations shouldn't start until *after* the start_time.
                #      This generation is tricky to verify, so skip it.
                continue
            if time == graph.demes[dest_id].start_time:
                # XXX: The migrants at the deme's start_time might represent
                #      the deme's ancestors.
                #      This generation is tricky to verify, so skip it.
                continue
            mrow = list(migration_matrix[dest_id])
            for pulse in pulses_at_time.get(time, []):
                if deme_id[pulse.dest] == dest_id:
                    mrow[deme_id[pulse.source]] += pulse.proportion
            assert sum(m > 0 for m in mrow) == len(state.immigrant_fracs), (
                f"{dest_id} ({deme_names[dest_id]}): {state.time} ({time}): "
                f"{mrow} inconsistent with {state.immigrant_pids} / {state.immigrant_fracs} "
                f"pulses={pulses_at_time.get(time, [])}\n"
                f"{state}\n"
                f"[{start_time}, {end_time}]:\n{migration_matrix}"
            )
            assert all(
                math.isclose(mrow[i], rate)
                for i, rate in zip(state.immigrant_pids, state.immigrant_fracs)
            ), (
                f"{dest_id} ({deme_names[dest_id]}): {state.time} ({time}): "
                f"{mrow} inconsistent with {state.immigrant_pids} / {state.immigrant_fracs} "
                f"pulses={pulses_at_time.get(time, [])}\n"
                f"{state}\n"
                f"[{start_time}, {end_time}]:\n{migration_matrix}"
            )

        start_time = end_time


def parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        description="validate demes_schedule_events() against the demes graph"
    )
    parser.add_argument(
        "-Q",
        "--scaling_factor",
        type=float,
        default=1.0,
        help="Scaling factor that was applied to speed up the simulation",
    )
    parser.add_argument(
        "slim_out_filename",
        metavar="slim-states.out",
        help="Output from test_demes_schedule_events.slim",
    )
    parser.add_argument(
        "yaml_filename",
        metavar="model.yaml",
        help="The YAML-formatted demes model being simulated.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    graph = demes.load(args.yaml_filename)
    discrete_graph = graph_to_scaled_discrete(graph, scaling_factor=args.scaling_factor)
    states = parse_slim_states(args.slim_out_filename)
    check_states_against_graph(states, discrete_graph)
