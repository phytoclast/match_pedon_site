From plots database, export to data folder "Observed_Species" table as tab delimited text with headers.
Copy "Sites" table directly to existing text file to ensure proper significant figures.
May need to updated NASISPEDONS.RDS using SoilDB fetchPedons. To do you you need to query NASIS for all pedons in the last 50 years for the Grand Rapids and Flint Offices.

Run "01.get_spatial_soil.R" to identify map units of new points.

Run "02.match_pedo_site_NASIS3.R" to get appropriate grouping of plots.

Run  03.USNVC_compare_specieslists_loop_by_plot.R in USNVC project

Run "04.aggregate_plot_summary_handpicked.R" to summarize species composition by group.




NASIS Query:
FROM site
INNER JOIN nasisgroup ON nasisgroup.grpiid=site.grpiidref AND  nasisgroup.group_name LIKE ? "NASIS Group Name, Use % for wildcard" 
INNER JOIN siteobs ON site.siteiid=siteobs.siteiidref AND siteobs.obsdate BETWEEN ? "Begin date" AND  ? "End date"
LEFT OUTER JOIN pedon ON siteobs.siteobsiid=pedon.siteobsiidref
WHERE site.latstddecimaldegrees BETWEEN ? "Minimum Latitude" AND  ? "Maximum Latitude" AND site.longstddecimaldegrees BETWEEN ? "Minimum Longitude" AND  ? "Maximum Longitude"



