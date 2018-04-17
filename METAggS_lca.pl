#!/usr/bin/perl

use List::Util qw(reduce);

$taxid = $ARGV[0];
$lcainput = $ARGV[1];
$lcaoutput = $ARGV[2];
$database = $ARGV[3];

open(FILE_NOD, "<${taxid}") or die "Coudn't open the file $taxid: $!";
open(FILE_MY, "<${lcainput}") or die "Coudn't open the file $lcainput: $!";
open(FILE_REF, "<${database}") or die "Coudn't open the file $database: $!";

open($file, ">${lcaoutput}") or die "Coudn't open the file $lcaoutput: $!";

while (defined(my $input = <FILE_NOD>)) {
    my @node = split(/\t/, $input);
    chomp($node[7]);
    push @{$offs{$node[0]}}, ($node[1],$node[2],$node[3],$node[4],$node[5],$node[6],$node[7]);
}

while (defined($input = <FILE_MY>)) {
    chomp($input);
    @contig = split(/\t/, $input);
    if (!(exists $taxa{$contig[1]})) {
        push @{$taxa{$contig[1]}}, 'x';
    }
    $taxa1{$contig[0]}{$contig[1]} = $contig[11];
}

while (defined($input = <FILE_REF>)) {
    ($ac0,$taxid) = split(/\t/, $input, 2);
    chomp($taxid);
    $ac = 'UniRef100_' . $ac0;
    if (defined $taxa{$ac}) {
        push @{$taxa{$ac}}, ($offs{$taxid}[0], $offs{$taxid}[1],$offs{$taxid}[2], $offs{$taxid}[3],$offs{$taxid}[4], $offs{$taxid}[5],$offs{$taxid}[6]);
    }
}

undef %offs;

foreach $key (sort keys %taxa1) {
    $flag = 0;
    print $file "$key";
    foreach $key1 (sort {$taxa1{$key}{$b} <=> $taxa1{$key}{$a}} keys %{$taxa1{$key}}) {r
        if (defined $taxa{$key1}[1]) {
            if ($flag == 0) {
                $score = $taxa1{$key}{$key1};
                $flag = 1;
                ${$final{$taxa{$key1}[1]}}[0] += ($taxa1{$key}{$key1}/$score)**20;
                ${$final{$taxa{$key1}[1]}}[1] = ${$final{$taxa{$key1}[1]}}[1] . ' ' . $key1;
            } else {
                ${$final{$taxa{$key1}[1]}}[0] += ($taxa1{$key}{$key1}/$score)**20;
                ${$final{$taxa{$key1}[1]}}[1] = ${$final{$taxa{$key1}[1]}}[1] . ' ' . $key1;
            }
        }
    }
    my $highest = reduce { ${$final{$a}}[0] > ${$final{$b}}[0] ? $a : $b } keys %final;
    print $file "\t$highest";
    foreach $z (2..7) {
        $flag = 0;
        @keep = split / /, ${$final{$highest}}[1];
        shift(@keep);
        undef %final;
        foreach $ke (@keep) {
            if ($flag == 0) {
                $score = $taxa1{$key}{$ke};
                $flag = 1;
                ${$final{$taxa{$ke}[$z]}}[0] += ($taxa1{$key}{$ke}/$score)**20;
                ${$final{$taxa{$ke}[$z]}}[1] = ${$final{$taxa{$ke}[$z]}}[1] . ' ' . $ke;
            } else {
                ${$final{$taxa{$ke}[$z]}}[0] += ($taxa1{$key}{$ke}/$score)**20;
                ${$final{$taxa{$ke}[$z]}}[1] = ${$final{$taxa{$ke}[$z]}}[1] . ' ' . $ke;
            }
        }
        $highest = reduce { ${$final{$a}}[0] > ${$final{$b}}[0] ? $a : $b } keys %final;
        print $file "\t$highest";
    }
    print $file "\n";
    undef %final;
}
