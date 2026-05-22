from pathlib import Path
import ast
import json
import math
import os
import re
import pandas as pd

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib_cache").resolve()))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


GRAPH_PATH = Path("GBM39_amplicon1_graph.txt")
CYCLES_PATH = Path("GBM39_amplicon1_cycles.txt")
CHROMOSOME = "chr7"
PLOT_START = 54729879
PLOT_END = 56219661
MAX_COPY_NUMBER = 116.21741886841724


def print_warning(path, line_number, line, reason):
    print(f"WARNING: {path.name}:{line_number}: {reason}: {line}")


def parse_float(value):
    if value in {"", "None"} or pd.isna(value):
        return None
    return float(value)


def parse_int(value):
    if value in {"", "None"} or pd.isna(value):
        return None
    return int(value)


def parse_endpoint(endpoint):
    match = re.match(r"^(?P<chrom>[^:]+):(?P<pos>-?\d+)(?P<strand>[+-])$", endpoint)
    if match is None:
        return None
    return {"chrom": match.group("chrom"), "pos": int(match.group("pos")), "strand": match.group("strand")}


def make_edge_row(raw_line, line_number, graph_section, edge_type, raw_edge_string, endpoint1, endpoint2, copy_number, extra_fields):
    row = {
        "line_number": line_number,
        "graph_section": graph_section,
        "raw_line": raw_line,
        "raw_edge_string": raw_edge_string,
        "chr1": endpoint1["chrom"],
        "pos1": endpoint1["pos"],
        "strand1": endpoint1["strand"],
        "chr2": endpoint2["chrom"],
        "pos2": endpoint2["pos"],
        "strand2": endpoint2["strand"],
        "copy_number": copy_number,
        "edge_type": edge_type,
    }
    row.update(extra_fields)
    return row


def parse_sequence_edge(path, line_number, line):
    fields = line.split("\t")
    if len(fields) != 7:
        print_warning(path, line_number, line, "expected 7 fields in SequenceEdge section")
        return None
    endpoint1 = parse_endpoint(fields[1])
    endpoint2 = parse_endpoint(fields[2])
    if endpoint1 is None or endpoint2 is None:
        print_warning(path, line_number, line, "could not parse sequence endpoint")
        return None
    extra_fields = {
        "average_coverage": parse_float(fields[4]),
        "size": parse_int(fields[5]),
        "number_reads_mapped": parse_float(fields[6]),
        "number_of_read_pairs": None,
        "homology_size": None,
        "homology_or_insertion_sequence": None,
    }
    return make_edge_row(line, line_number, "SequenceEdge", fields[0], f"{fields[1]}->{fields[2]}", endpoint1, endpoint2, parse_float(fields[3]), extra_fields)


def parse_breakpoint_edge(path, line_number, line):
    fields = line.split("\t")
    if len(fields) != 6:
        print_warning(path, line_number, line, "expected 6 fields in BreakpointEdge section")
        return None
    if "->" not in fields[1]:
        print_warning(path, line_number, line, "missing breakpoint separator")
        return None
    endpoint_text1, endpoint_text2 = fields[1].split("->", 1)
    endpoint1 = parse_endpoint(endpoint_text1)
    endpoint2 = parse_endpoint(endpoint_text2)
    if endpoint1 is None or endpoint2 is None:
        print_warning(path, line_number, line, "could not parse breakpoint endpoint")
        return None
    extra_fields = {
        "average_coverage": None,
        "size": None,
        "number_reads_mapped": None,
        "number_of_read_pairs": parse_float(fields[3]),
        "homology_size": fields[4] if fields[4] != "None" else None,
        "homology_or_insertion_sequence": fields[5] if fields[5] != "None" else None,
    }
    return make_edge_row(line, line_number, "BreakpointEdge", fields[0], fields[1], endpoint1, endpoint2, parse_float(fields[2]), extra_fields)


def parse_graph_file(path):
    rows = []
    section = None
    for line_number, raw_line in enumerate(path.read_text().splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("SequenceEdge:"):
            section = "SequenceEdge"
            continue
        if line.startswith("BreakpointEdge:"):
            section = "BreakpointEdge"
            continue
        if section == "SequenceEdge":
            row = parse_sequence_edge(path, line_number, line)
        elif section == "BreakpointEdge":
            row = parse_breakpoint_edge(path, line_number, line)
        else:
            print_warning(path, line_number, line, "line before recognized graph section")
            row = None
        if row is not None:
            rows.append(row)
    return pd.DataFrame(rows)


def parse_oriented_segments(segment_text):
    parsed_segments = []
    tokens = [token.strip() for token in segment_text.split(",") if token.strip()]
    for token in tokens:
        match = re.match(r"^(?P<segment_id>\d+)(?P<orientation>[+-])$", token)
        if match is None:
            return None
        parsed_segments.append({"segment_id": int(match.group("segment_id")), "orientation": match.group("orientation")})
    return parsed_segments


def parse_interval_line(path, line_number, line):
    fields = line.split("\t")
    if len(fields) != 5:
        print_warning(path, line_number, line, "expected 5 fields in Interval line")
        return None
    return {"line_number": line_number, "raw_line": line, "interval_id": parse_int(fields[1]), "chr": fields[2], "start": parse_int(fields[3]), "end": parse_int(fields[4])}


def parse_segment_line(path, line_number, line):
    fields = line.split("\t")
    if len(fields) != 5:
        print_warning(path, line_number, line, "expected 5 fields in Segment line")
        return None
    return {"line_number": line_number, "raw_line": line, "segment_id": parse_int(fields[1]), "chr": fields[2], "start": parse_int(fields[3]), "end": parse_int(fields[4])}


def parse_cycle_line(path, line_number, line):
    match = re.match(r"^Cycle=(?P<cycle_id>\d+);Copy_count=(?P<copy_count>[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?);Segments=(?P<segments>.*)$", line)
    if match is None:
        print_warning(path, line_number, line, "could not parse Cycle line")
        return None
    parsed_segments = parse_oriented_segments(match.group("segments"))
    if parsed_segments is None:
        print_warning(path, line_number, line, "could not parse oriented segment list")
        return None
    return {
        "line_number": line_number,
        "raw_line": line,
        "cycle_id": parse_int(match.group("cycle_id")),
        "copy_count": parse_float(match.group("copy_count")),
        "segment_list": match.group("segments"),
        "segment_ids": json.dumps([item["segment_id"] for item in parsed_segments]),
        "orientations": json.dumps([item["orientation"] for item in parsed_segments]),
        "oriented_segments": json.dumps(parsed_segments),
        "segment_count": len(parsed_segments),
    }


def parse_cycles_file(path):
    interval_rows = []
    segment_rows = []
    cycle_rows = []
    for line_number, raw_line in enumerate(path.read_text().splitlines(), start=1):
        line = raw_line.strip()
        if not line or line == "List of cycle segments":
            continue
        if line.startswith("Interval\t"):
            row = parse_interval_line(path, line_number, line)
            if row is not None:
                interval_rows.append(row)
            continue
        if line.startswith("Segment\t"):
            row = parse_segment_line(path, line_number, line)
            if row is not None:
                segment_rows.append(row)
            continue
        if line.startswith("Cycle="):
            row = parse_cycle_line(path, line_number, line)
            if row is not None:
                cycle_rows.append(row)
            continue
        print_warning(path, line_number, line, "unrecognized cycles file line")
    return pd.DataFrame(interval_rows), pd.DataFrame(segment_rows), pd.DataFrame(cycle_rows)


def classify_edge(row):
    if row["graph_section"] == "SequenceEdge":
        return "sequence_edge"
    if row["graph_section"] != "BreakpointEdge":
        return "unclear"
    if pd.isna(row["chr1"]) or pd.isna(row["chr2"]) or pd.isna(row["pos1"]) or pd.isna(row["pos2"]):
        return "unclear"
    if row["pos1"] < 0 or row["pos2"] < 0:
        return "unclear"
    if row["chr1"] == row["chr2"] and row["pos1"] == row["pos2"]:
        return "self_loop"
    if row["chr1"] == row["chr2"]:
        return "breakpoint_intrachromosomal"
    return "breakpoint_interchromosomal"


def classify_graph_edges(graph_edges):
    classified = graph_edges.copy()
    classified["plotting_category"] = classified.apply(classify_edge, axis=1)
    return classified


def build_sequence_segment_rows(classified_edges):
    rows = []
    sequence_edges = classified_edges[classified_edges["plotting_category"] == "sequence_edge"].copy()
    for _, row in sequence_edges.iterrows():
        rows.append({
            "chromosome": row["chr1"],
            "start": min(row["pos1"], row["pos2"]),
            "end": max(row["pos1"], row["pos2"]),
            "copy_number": row["copy_number"],
            "segment_id": pd.NA,
            "source": "SequenceEdge",
            "raw_edge": row["raw_edge_string"],
            "draw_as_horizontal_segment": True,
        })
    return pd.DataFrame(rows)


def attach_matching_segment_ids(segment_plot, cycle_segments):
    result = segment_plot.copy()
    for _, segment in cycle_segments.iterrows():
        mask = (result["chromosome"] == segment["chr"]) & (result["start"] == segment["start"]) & (result["end"] == segment["end"])
        result.loc[mask, "segment_id"] = segment["segment_id"]
    return result


def build_cycle_segment_rows(cycle_segments, sequence_segment_rows):
    rows = []
    for _, segment in cycle_segments.iterrows():
        exact_match = sequence_segment_rows[(sequence_segment_rows["chromosome"] == segment["chr"]) & (sequence_segment_rows["start"] == segment["start"]) & (sequence_segment_rows["end"] == segment["end"])]
        copy_number = exact_match["copy_number"].iloc[0] if len(exact_match) == 1 else pd.NA
        rows.append({
            "chromosome": segment["chr"],
            "start": segment["start"],
            "end": segment["end"],
            "copy_number": copy_number,
            "segment_id": segment["segment_id"],
            "source": "cycle_segments",
            "raw_edge": segment["raw_line"],
            "draw_as_horizontal_segment": True,
        })
    return pd.DataFrame(rows)


def build_segment_plot_df(classified_edges, cycle_segments):
    sequence_segment_rows = attach_matching_segment_ids(build_sequence_segment_rows(classified_edges), cycle_segments)
    cycle_segment_rows = build_cycle_segment_rows(cycle_segments, sequence_segment_rows)
    return pd.concat([sequence_segment_rows, cycle_segment_rows], ignore_index=True)


def should_draw_breakpoint_arc(row):
    return row["plotting_category"] in {"breakpoint_intrachromosomal", "breakpoint_interchromosomal", "self_loop"} and row["edge_type"] == "discordant"


def build_breakpoint_plot_df(classified_edges):
    rows = []
    breakpoint_edges = classified_edges[classified_edges["graph_section"] == "BreakpointEdge"].copy()
    for _, row in breakpoint_edges.iterrows():
        valid_positions = row["pos1"] >= 0 and row["pos2"] >= 0
        rows.append({
            "chr1": row["chr1"],
            "pos1": row["pos1"],
            "strand1": row["strand1"],
            "chr2": row["chr2"],
            "pos2": row["pos2"],
            "strand2": row["strand2"],
            "copy_number": row["copy_number"],
            "breakpoint_type": row["plotting_category"],
            "raw_edge": row["raw_edge_string"],
            "edge_type": row["edge_type"],
            "draw_as_arc": should_draw_breakpoint_arc(row),
            "arc_start": min(row["pos1"], row["pos2"]) if valid_positions else pd.NA,
            "arc_end": max(row["pos1"], row["pos2"]) if valid_positions else pd.NA,
        })
    breakpoint_plot = pd.DataFrame(rows)
    breakpoint_plot["arc_height"] = pd.NA
    arc_rows = breakpoint_plot[breakpoint_plot["draw_as_arc"]].copy()
    arc_rows["arc_span"] = arc_rows["arc_end"] - arc_rows["arc_start"]
    arc_rows = arc_rows.sort_values(["arc_span", "copy_number"], ascending=[False, False])
    base_height = MAX_COPY_NUMBER * 1.08
    height_step = MAX_COPY_NUMBER * 0.14
    for rank, index in enumerate(arc_rows.index, start=1):
        breakpoint_plot.loc[index, "arc_height"] = base_height + rank * height_step
    return breakpoint_plot


def build_cycle_plot_df(cycles):
    rows = []
    for _, row in cycles.iterrows():
        segment_ids = json.loads(row["segment_ids"])
        orientations = json.loads(row["orientations"])
        rows.append({
            "cycle_id": row["cycle_id"],
            "copy_count": row["copy_count"],
            "ordered_segment_list": row["segment_list"],
            "segment_ids": json.dumps(segment_ids),
            "orientations": json.dumps(orientations),
            "segment_count": row["segment_count"],
        })
    return pd.DataFrame(rows)


def get_sequence_segments(segment_plot):
    return segment_plot[(segment_plot["source"] == "SequenceEdge") & (segment_plot["draw_as_horizontal_segment"])].copy()


def get_cycle_segments(segment_plot):
    return segment_plot[(segment_plot["source"] == "cycle_segments") & (segment_plot["segment_id"].notna())].copy()


def get_valid_arcs(breakpoint_plot):
    arc_rows = breakpoint_plot[
        (breakpoint_plot["draw_as_arc"])
        & (breakpoint_plot["breakpoint_type"] == "breakpoint_intrachromosomal")
        & (breakpoint_plot["chr1"] == CHROMOSOME)
        & (breakpoint_plot["chr2"] == CHROMOSOME)
        & (breakpoint_plot["pos1"] >= 0)
        & (breakpoint_plot["pos2"] >= 0)
    ].copy()
    arc_rows["arc_span"] = (arc_rows["arc_end"] - arc_rows["arc_start"]).abs()
    return arc_rows.sort_values(["arc_span", "copy_number"], ascending=[False, False])


def parse_segment_ids(value):
    return [int(item) for item in ast.literal_eval(value)]


def format_position(value):
    return f"{value / 1_000_000:.2f} Mb"


def build_cycle_membership(cycle_plot):
    membership = {}
    for _, row in cycle_plot.iterrows():
        cycle_id = int(row["cycle_id"])
        for segment_id in parse_segment_ids(row["segment_ids"]):
            membership.setdefault(segment_id, []).append(cycle_id)
    return membership


def assign_arc_display_heights(arc_rows):
    display = arc_rows.copy()
    display["display_height"] = pd.NA
    heights = [142, 156, 170, 184, 198]
    for rank, index in enumerate(display.index):
        display.loc[index, "display_height"] = heights[min(rank, len(heights) - 1)]
    return display


def draw_arc(ax, start, end, height, color, linewidth, alpha):
    midpoint = (start + end) / 2
    radius = (end - start) / 2
    baseline = MAX_COPY_NUMBER + 5
    xs = []
    ys = []
    for step in range(121):
        theta = math.pi * step / 120
        xs.append(midpoint - radius * math.cos(theta))
        ys.append(baseline + (height - baseline) * math.sin(theta))
    ax.plot(xs, ys, color=color, linewidth=linewidth, alpha=alpha, solid_capstyle="round")


def draw_sequence_segments(ax, sequence_segments):
    for _, row in sequence_segments.iterrows():
        linewidth = 7 if row["copy_number"] >= 20 else 5
        alpha = 0.95 if row["copy_number"] >= 20 else 0.75
        ax.plot([row["start"], row["end"]], [row["copy_number"], row["copy_number"]], color="#2468a8", linewidth=linewidth, alpha=alpha, solid_capstyle="butt")


def draw_cycle_highlights(ax, cycle_segments, cycle_membership):
    colors = {1: "#f2b134", 2: "#6abf69"}
    y_ranges = {1: (-13, -8), 2: (-20, -15)}
    for _, row in cycle_segments.iterrows():
        segment_id = int(row["segment_id"])
        for cycle_id in cycle_membership.get(segment_id, []):
            ymin, ymax = y_ranges[cycle_id]
            ax.fill_between([row["start"], row["end"]], ymin, ymax, color=colors[cycle_id], alpha=0.85, linewidth=0)


def draw_segment_labels(ax, cycle_segments, labeled):
    label_y = MAX_COPY_NUMBER + 9
    for _, row in cycle_segments.iterrows():
        segment_id = int(row["segment_id"])
        midpoint = (row["start"] + row["end"]) / 2
        if labeled or segment_id in {1, 2, 4}:
            ax.text(midpoint, label_y, f"S{segment_id}", ha="center", va="bottom", fontsize=9, color="#1e2a32", fontweight="bold")


def draw_breakpoint_arcs(ax, arc_rows, labeled):
    display_arcs = assign_arc_display_heights(arc_rows)
    for rank, (_, row) in enumerate(display_arcs.iterrows(), start=1):
        draw_arc(ax, row["arc_start"], row["arc_end"], float(row["display_height"]), "#b04a1f", 2.1, 0.9)
        if labeled:
            midpoint = (row["arc_start"] + row["arc_end"]) / 2
            ax.text(midpoint, float(row["display_height"]) + 4, f"BP{rank}\nCN {row['copy_number']:.1f}", ha="center", va="bottom", fontsize=8, color="#7a2d12")


def draw_cycle_annotation_box(ax, cycle_plot):
    lines = ["Cycles"]
    for _, row in cycle_plot.iterrows():
        lines.append(f"Cycle {int(row['cycle_id'])}: CN {row['copy_count']:.1f}, {row['ordered_segment_list']}")
    ax.text(0.015, 0.965, "\n".join(lines), transform=ax.transAxes, ha="left", va="top", fontsize=10, bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#c9c9c9", "alpha": 0.95})


def style_main_axis(ax):
    ax.set_xlim(PLOT_START, PLOT_END)
    ax.set_ylim(-24, 185)
    ax.set_title("GBM39 Amplicon 1 Recreated Structure", fontsize=15, pad=14)
    ax.set_xlabel("chr7 genomic coordinate")
    ax.set_ylabel("copy number")
    ticks = [PLOT_START, 54800000, 55000000, 55200000, 55400000, 55600000, 55800000, 56000000, PLOT_END]
    ax.set_xticks(ticks)
    ax.set_xticklabels([format_position(tick) for tick in ticks], rotation=35, ha="right")
    ax.set_yticks([0, 25, 50, 75, 100, 125, 150, 175])
    ax.axhline(0, color="#323232", linewidth=0.8)
    ax.grid(axis="y", color="#e7e7e7", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def add_legend(ax):
    handles = [
        Line2D([0], [0], color="#2468a8", linewidth=6, label="Sequence edge copy number"),
        Line2D([0], [0], color="#b04a1f", linewidth=2, label="Discordant breakpoint arc"),
        Patch(facecolor="#f2b134", alpha=0.85, label="Cycle 1 segments"),
        Patch(facecolor="#6abf69", alpha=0.85, label="Cycle 2 segments"),
    ]
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=9)


def add_cycle_track_labels(ax):
    ax.text(PLOT_START, -10.5, "Cycle 1", ha="left", va="center", fontsize=8, color="#7a5a00")
    ax.text(PLOT_START, -17.5, "Cycle 2", ha="left", va="center", fontsize=8, color="#2d6f31")


def build_figure(segment_plot, breakpoint_plot, cycle_plot, labeled=False):
    sequence_segments = get_sequence_segments(segment_plot)
    cycle_segments = get_cycle_segments(segment_plot)
    arc_rows = get_valid_arcs(breakpoint_plot)
    cycle_membership = build_cycle_membership(cycle_plot)
    fig, ax = plt.subplots(figsize=(12, 6))
    draw_cycle_highlights(ax, cycle_segments, cycle_membership)
    draw_sequence_segments(ax, sequence_segments)
    draw_segment_labels(ax, cycle_segments, labeled)
    draw_breakpoint_arcs(ax, arc_rows, labeled)
    draw_cycle_annotation_box(ax, cycle_plot)
    add_cycle_track_labels(ax)
    style_main_axis(ax)
    add_legend(ax)
    fig.tight_layout()
    return fig


def save_final_figures(segment_plot, breakpoint_plot, cycle_plot):
    main_fig = build_figure(segment_plot, breakpoint_plot, cycle_plot, labeled=False)
    main_fig.savefig("recreated_amplicon_plot.png", dpi=300)
    main_fig.savefig("recreated_amplicon_plot.pdf")
    plt.close(main_fig)
    labeled_fig = build_figure(segment_plot, breakpoint_plot, cycle_plot, labeled=True)
    labeled_fig.savefig("recreated_amplicon_plot_labeled.png", dpi=300)
    plt.close(labeled_fig)


def save_tables(graph_edges, cycle_segments, cycles, segment_plot, breakpoint_plot, cycle_plot):
    graph_edges.to_csv("graph_edges.csv", index=False)
    cycle_segments.to_csv("cycle_segments.csv", index=False)
    cycles.to_csv("cycles.csv", index=False)
    segment_plot.to_csv("segment_plot_data.csv", index=False)
    breakpoint_plot.to_csv("breakpoint_plot_data.csv", index=False)
    cycle_plot.to_csv("cycle_plot_data.csv", index=False)


def print_final_summary(graph_edges, breakpoint_plot, cycle_plot):
    sequence_count = int((graph_edges["graph_section"] == "SequenceEdge").sum())
    breakpoint_count = int((graph_edges["graph_section"] == "BreakpointEdge").sum())
    valid_arcs = get_valid_arcs(breakpoint_plot)
    ambiguous_edges = breakpoint_plot[(breakpoint_plot["breakpoint_type"] == "unclear") | (breakpoint_plot["pos1"] < 0) | (breakpoint_plot["pos2"] < 0)]
    print("GBM39 Project 2 final workflow complete")
    print(f"Sequence edges: {sequence_count}")
    print(f"Breakpoint edges: {breakpoint_count}")
    print(f"Valid breakpoint arcs drawn: {len(valid_arcs)}")
    print(f"Ambiguous source edges excluded from arcs: {len(ambiguous_edges)}")
    print(f"Plot focus: {CHROMOSOME}:{PLOT_START}-{PLOT_END}")
    for _, row in cycle_plot.iterrows():
        print(f"Cycle {int(row['cycle_id'])}: copy_count={row['copy_count']:.4f}, segments={row['ordered_segment_list']}")


def run_workflow():
    graph_edges = parse_graph_file(GRAPH_PATH)
    intervals, cycle_segments, cycles = parse_cycles_file(CYCLES_PATH)
    classified_edges = classify_graph_edges(graph_edges)
    segment_plot = build_segment_plot_df(classified_edges, cycle_segments)
    breakpoint_plot = build_breakpoint_plot_df(classified_edges)
    cycle_plot = build_cycle_plot_df(cycles)
    save_tables(graph_edges, cycle_segments, cycles, segment_plot, breakpoint_plot, cycle_plot)
    save_final_figures(segment_plot, breakpoint_plot, cycle_plot)
    print_final_summary(graph_edges, breakpoint_plot, cycle_plot)
    return graph_edges, cycle_segments, cycles, segment_plot, breakpoint_plot, cycle_plot


if __name__ == "__main__":
    run_workflow()
