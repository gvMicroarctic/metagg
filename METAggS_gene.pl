#!/usr/bin/perl

open(FILE, "<$ARGV[0]") or die "Coudn't open the file $ARGV[0]: $!";
open($file, ">$ARGV[1]") or die "Coudn't open the file $ARGV[1]: $!";

$a = 0;

while (defined($input = <FILE>)) {
    chomp($input);
    @info = split(/\t/, $input);
    $inside = 0;
    if ($info[6] < $info[7]) {
        $start=$info[6];
        $end=$info[7];
    } else {
        $start=$info[7];
        $end=$info[6];
    }
    if ($a == 0) {
        $a++;
        $final{$a}{'S'}=$start;
        $final{$a}{'E'}=$end;
        $finalprint{$start} = $input;
    } else {
        foreach $key (sort keys %final) {
            if ((($start <= $final{$key}{'S'}) && ($end >= $final{$key}{'S'})) or (($start <= $final{$key}{'E'}) && ($end >= $final{$key}{'E'})) or (($start >= $final{$key}{'S'}) && ($end <= $final{$key}{'E'})) && ($inside == 0)) {
                $inside = 1;
            }
        }
        if ($inside == 0) {
            $a++;
            $final{$a}{'S'}=$start;
            $final{$a}{'E'}=$end;
            $finalprint{$start} = $input;
        }
    }
}


foreach $key (sort {$a <=> $b} keys %finalprint) {
    $b++;
    print $file "$finalprint{$key}\n";
}

print "Total number of genes: $b\n";
