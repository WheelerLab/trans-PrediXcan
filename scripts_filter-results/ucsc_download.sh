# The script is intended to demonstrate how to fetch data from UCSC easily
# (without mysql) in a way which will enable reproducibility of the download.

# It is NOT INTENDED as a generic tool for programmatic access to UCSC;
# Please respect UCSC request of avoiding issuing heavy queries as described in mysql usage manual;
# for more information on conditions of use, see: http://genome.ucsc.edu/goldenPath/help/mysql.html

# For advanced programmatic access you may be interested in setting you own local copy of UCSC database
# and using other, existing tools like: PyUCSC (see: https://pyucsc.readthedocs.io/en/latest)

# Licence: MIT

# Example usage:

# You can load functions from this file with:
# source ucsc_download.sh

# Example 1. Download RefSeq mRNA descriptions to a gzipped file:
# get_whole_genome_table refseq_summary.tsv.gz genes refGene hgFixed.refSeqSummary gzip

# Example 2. Fetch data and pass to some other commad through a pipe
# get_whole_genome_table - genes refGene hgFixed.refSeqSummary | grep protease


unset sid
hg_tables_url="http://genome.ucsc.edu/cgi-bin/hgTables"


function fetch_sid {
	sid=$(wget $hg_tables_url -O - | sed -n "s/.*NAME='hgsid' VALUE='\(.*\)'.*/\1/p" | head -n 1)
}


function get_from_ucsc {
	if [ -z $sid ]; then fetch_sid; fi;
 	query=$(sed "s/\(.*\):\(.*\)/\1=\2/" | tr '\n' '&' $2 | sed 's/\t*//g' | tr ' ' + | sed 's/,/%0D%0A/g')
	echo "hgsid=$sid&${query:0:-1}"
	wget $hg_tables_url --post-data "hgsid=$sid&${query:0:-1}" -O $1
}


function get_whole_genome_table {
	filename=${1:-ucsc_table}
	compression=${5:-none}
	db=${6:-hg19}

	get_from_ucsc $1 <<-QUERY
		jsh_pageVertPos:0
		position:chr21:33031597-33041570
		clade:mammal
		org:Human
		db:$db
		hgta_group:$2
		hgta_track:$3
		hgta_table:$4
		hgta_regionType:genome
		hgta_outputType:primaryTable
		boolshad.sendToGalaxy:0
		boolshad.sendToGreat:0
		boolshad.sendToGenomeSpace:0
		hgta_outFileName:output
		hgta_compressType:$compression
		hgta_doTopSubmit:get output
	QUERY
}