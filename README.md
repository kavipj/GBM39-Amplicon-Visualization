Project title: CSE 182 Project 2 GBM39 Amplicon 1 Reconstruction

Sample name: GBM39
Amplicon name: GBM39 Amplicon 1

Files used:
- GBM39_amplicon1_graph.txt
- GBM39_amplicon1_cycles.txt
- GBM39_amplicon1.png
- GBM39_amplicon1.pdf
- GBM39_summary.txt

Main final files:
- GBM39_Project2_final.ipynb
- GBM39_Project2_final.py
- recreated_amplicon_plot.png
- recreated_amplicon_plot.pdf
- recreated_amplicon_plot_labeled.png

Intermediate output tables:
- graph_edges.csv
- cycle_segments.csv
- cycles.csv
- segment_plot_data.csv
- breakpoint_plot_data.csv
- cycle_plot_data.csv

How to run:
1. Open GBM39_Project2_final.ipynb in Jupyter.
2. Run all cells from top to bottom.
3. The notebook parses the graph and cycles files, creates plotting-ready tables, and saves the final figure outputs.

Alternative command-line run:
python3 GBM39_Project2_final.py

What the final figure shows:
The final figure recreates the GBM39 Amplicon 1 structure on chr7 from coordinate 54729879 to 56219661. Blue horizontal lines represent sequence edges placed at their estimated copy-number values. Brown arcs represent valid intrachromosomal discordant breakpoint edges. Colored tracks below the x-axis show which cycle segments are used by Cycle 1 and Cycle 2.

Cycle summary:
- Cycle 1: copy count 99.5734741489, segment order 2+,4+,1+
- Cycle 2: copy count 13.3171554681, segment order 3+,1+

Assumptions and notes:
- The visualization focuses on chr7 because all parsed segments and valid breakpoints are on chr7.
- The 4 source breakpoint rows with pos1 = -1 are kept in the data tables but excluded from arc drawing because they do not have complete genomic coordinates.
- Concordant one-base adjacency edges are not drawn as arcs because they represent neighboring sequence continuity rather than structural breakpoint arcs.
- The provided GBM39_amplicon1.png and GBM39_amplicon1.pdf were used only as visual references to check that the recreated structure is reasonable.
- The main figure is cleaner for a final writeup. The labeled figure is useful if breakpoint copy numbers need to be visible directly on the plot.
