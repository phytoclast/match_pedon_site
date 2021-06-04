From plots database, export to data folder "Observed_Species" table as tab delimited text with headers.
Copy "Sites" table directly to existing text file to ensure proper significant figures.
May need to updata NASISPEDONS.RDS using SoilDB fetchPedons. To do you you need to query NASIS for all pedons in the last 50 years for the Grand Rapids and Flint Offices.

Run "get_spatial_soil.R" to identify map units of new points.

Run "match_pedo_site_NASIS3.R" to get appropriate grouping of plots.

Run  USNVC_compare_specieslists_loop_by_plot.R in USNVC project

Run "aggregate_plot_summary_handpicked.R" to summarize species composition by group.



