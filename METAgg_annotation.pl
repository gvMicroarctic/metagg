#!/usr/bin/perl 

$infile = $ARGV[0];
($outfile) = $infile =~ /(.*)\./;
$database = $ARGV[1];

`mkdir ${outfile}_tmp`;

if (defined $ARGV[2]) {
	$evalue = $ARGV[2];
} else {
	$evalue = 0.0001;
}

if (defined $ARGV[3]) {
	$maxmatch = $ARGV[3];
} else {
	$maxmatch = 50;
}

`diamond blastx -d $database -q $infile -o ${outfile}2.m8 -t ${outfile}_tmp --un unaligned.txt --strand plus -e $evalue --id 40 -f 6`;

`rm -r ${outfile}_tmp`;
