#!/usr/bin/perl

@file = `ls *_bwacut20less1.sam`;

open($file, ">TOT_SNPS.txt") or die "Coudn't open the file $file: $!";
open($file1, ">TOT_INSERTION_SNPS1.txt") or die "Coudn't open the file $file: $!";
open($file2, ">TOT_DELETION_SNPS1.txt") or die "Coudn't open the file $file: $!";
open($file0, ">TOT_base.txt") or die "Coudn't open the file TOT_base.txt: $!";

foreach $fi (@file) {
    chomp($fi);
    ($name_sample) = $fi =~ /^(.*?)_/;
    print "Processing $fi\n";
    print $file "\@${fi}";
    print $file0 "\@${fi}";
    print $file1 "\t$name_sample";
    print $file2 "\t$name_sample";
    push @name, $name_sample;
    foreach $big (1..2) {
        open(FILE_SAM, "<$fi") or die "Coudn't open the file $fi: $!";
        $tot_read = 0;
        print "Cycle $big\n";
        while (defined($input = <FILE_SAM>)) {
            $tot_read++;
            ($con,$read,$start_al,$cigar,$md,$seq) = split(/\t/, $input,6);
            $delins = 0;
            @cigar_snp = split /(?<=\M|I|D|S|H)/, $cigar;
            $match_read = 0;
            foreach $c (@cigar_snp) {
                if ($c =~ /M/) {
                    ($match) = $c =~ /(.*)[A-Z]/;
                    $match_read = $match_read + $match;
                } elsif ($c =~ /I/) {
                    ($insertion) = $c =~ /(.*)[A-Z]/;
                    $triple = substr ($seq, $match_read, $insertion, "");
                    if (($insertion/3) =~ /^\d+$/) {
                        while($triple =~ /([CATGN]{3})/g) {
                            $codon = $1;
                            $final_ins{$codon}{$name_sample}++;
                        }
                    }
                }
            }
            @md_snp = split /(?<=A|T|C|G|\^)/, $md;
            $match_tot = 0;
            $match_contig = $start_al;
            $del = 0;
            $deletion = "";
            foreach $m (@md_snp) {
                if ($m =~ /A|T|C|G/) {
                    if (($del > 0) && ($m !~ /[0-9]/)) {
                        $deletion = $deletion . $m;
                        $del++;
                        if ($big == 1) {
                            $onlysnps{$con}{$match_contig} = '';
                        } else {
                            $snsn{$match_contig} = '';
                        }
                        $delins++;
                        $match_contig++;
                    } else {
                        $del = 0;
                        ($match) = $m =~ /(.*)[A-Z]/;
                        $let = chop($m);
                        $match_tot = $match_tot + $match;
                        $match_contig = $match_contig + $match;
                        $snp = substr($seq, $match_tot, 1);
                        if ($big == 1) {
                            $con_read{$con}{$match_contig}{$let}{$snp}++;
                            $onlysnps{$con}{$match_contig} = '';
                        } else {
                            $snsn{$match_contig} = '';
                        }
                        $match_tot++;
                        $match_contig++;
                    }
                } elsif ($m =~ /\^/) {
                    $del = 1;
                    ($match) = $m =~ /(.*)\^/;
                    $match_tot = $match_tot + $match;
                    $match_contig = $match_contig +$match;
                    if ($big == 1) {
                        $onlysnps{$con}{$match_contig} = '';
                    } else {
                        $snsn{$match_contig} = '';
                    }
                }
            }
            if ((($del-1)/3) =~ /^\d+$/) {
                while($deletion =~ /([CATGN]{3})/g) {
                    $codon = $1;
                    if ($big == 1) {
                        $final_del{$codon}{$name_sample}++;
                    }
                }
            }
            if ($big == 2) {
                $seq_l = length($seq);
                $seq_l1 = $seq_l + $start_al - 2 + $delins;
                foreach $ll ($start_al..$seq_l1) {
                    if (!(exists $snsn{$ll}) && (exists $onlysnps{$con}{$ll})) {
                        $def_onlysnps{$con}{$ll}++;
                    }
                }
                undef %snsn;
            }
        }
        if ($big == 1) {
            $TOT{$name_sample} = $tot_read;
            foreach $contig (sort keys %con_read) {
                print $file "\n\>${contig}\n";
                foreach $pos (sort {$a <=> $b} keys %{$con_read{$contig}}) {
                    foreach $letcon (sort keys %{$con_read{$contig}{$pos}}) {
                        print $file "${pos}${letcon}:";
                        foreach $letread (sort keys %{$con_read{$contig}{$pos}{$letcon}}) {
                            print $file "${letread}$con_read{$contig}{$pos}{$letcon}{$letread};";
                        }
                    }
                }
            }
            undef %con_read;
            close FILE_SAM;
            print $file "\n";
        } else {
            foreach $b1 (sort keys %def_onlysnps) {
                print $file0 "\n\>$b1\n";
                foreach $b2 (sort keys %{$def_onlysnps{$b1}}) {
                    print $file0 "$b2:$def_onlysnps{$b1}{$b2}\;";
                }
            }
            close FILE_SAM;
        }
    }
}

undef %onlysnp;

foreach $codon (sort keys %final_del) {
    print $file1 "\n$codon";
    print $file2 "\n$codon";
    foreach $sample (@name) {
        if (exists $final_ins{$codon}{$sample}) {
            $fraction = ($final_ins{$codon}{$sample}/$TOT{$sample})*100;
            print $file1 "\;";
            printf $file1 "%.4f", $fraction;
        } else {
            print $file1 "\;0";
        }
        if (exists $final_del{$codon}{$sample}) {
            $fraction = ($final_del{$codon}{$sample}/$TOT{$sample})*100;
            print $file2 "\;";
            printf $file2 "%.4f", $fraction;
        } else {
            print $file2 "\;0";
        }
    }
}

close $file;
close $file1;
close $file2;

undef %final_del;
undef %final_ins;
print "Done with creation main file.\n";

open(SNP, "<TOT_SNPS.txt") or die "Coudn't open the file TOT_SNPS.txt: $!";
open(CONTIG, "<combined_trimmed.fna") or die "Coudn't open the file combined.fna: $!";
open(MATCH, "<combined_matches_mod_cut.m8") or die "Coudn't open the file combined_matches_mod_cut.m8: $!";

open($file, ">TOTSUM_SNPS.txt") or die "Coudn't open the file TOTSUM_SNPS.txt: $!";
open($file1, ">AA_SNPs1.txt") or die "Coudn't open the file AA_SNPs.txt: $!";
open($file2, ">CODON_SNPs1.txt") or die "Coudn't open the file CODON_SNPs.txt: $!";

$count = 0;
while (defined($input = <SNP>)) {
    if ($input =~ /^>/) {
        ($contig) = $input =~ />(.*)/;
        if (!(exists $index{$contig})) {
            $count++;
            $index{$contig} = '';
        }
    }
}

$cicle = 200;
$round = int($count/$cicle) + 1;
undef %index;

print "The total number of contig is $count\n";

print "The data have een divided in group of $round contigs\n";

while (defined($input=<MATCH>)) {
    @info = split /\t/, $input;
    $diamond{$info[0]}{$info[1]}{'S'} = $info[6];
    $diamond{$info[0]}{$info[1]}{'E'} = $info[7];
}
close MATCH;

$count1 = 0;
$group = 1;
while (defined($input = <CONTIG>)) {
    if ($input =~ /^>/) {
        ($contig) = $input =~ />(.*)/;
    } else {
        $count1++;
        if (!(exists $diamond{$contig})) {
            $definitive{$group}{$contig}='';
        } else {
            foreach $gene (keys %{$diamond{$contig}}) {
                if ($diamond{$contig}{$gene}{'S'} < $diamond{$contig}{$gene}{'E'}) {
                    $start = $diamond{$contig}{$gene}{'S'};
                    $end = $diamond{$contig}{$gene}{'E'};
                } else {
                    $start = $diamond{$contig}{$gene}{'E'};
                    $end = $diamond{$contig}{$gene}{'S'};
                    $revert{$contig}{$start} = '';
                }
                $start1 = $start-1;
                $seqp = substr $input, $start1, ($end-$start1);
                $definitive{$group}{$contig}{$start}{$end} = $seqp;
            }
        }
    }
    if ($count1 == $round) {
        $count1 = 0;
        $group++;
    }
}

close CONTIG;
undef %diamond;

print "Done with combined file\n";

$equal = 0;
$totcodonaa = 0;
foreach $i (sort {$a <=> $b} keys %definitive) {
    open(SNP, "<TOT_SNPS.txt") or die "Coudn't open the file TOT_SNPS.txt: $!";
    while (defined($input = <SNP>)) {
        chomp($input);
        if ($input =~ /^>/) {
            ($contig) = $input =~ />(.*)/;
            if (exists $definitive{$i}{$contig}) {
                $flag = 1;
            }
        } elsif ($input =~ /^@/) {
            ($name_sample) = $input =~ /^@(.*?)_/;
            $sampleindex{$name_sample} = '';
        } else {
            if ($flag == 1) {
                @SNP = split /;/, $input;
                foreach $s (@SNP) {
                    if ($s =~ /:/) {
                        @original = split /:/, $s;
                        $lcon =chop($original[0]);
                        $lread = substr($original[1], 0, 1,"");
                        $final{$contig}{$original[0]}{$lcon}{$lread}{$name_sample} += $original[1];
                    } else {
                        $lread = substr($s, 0, 1,"");
                        $final{$contig}{$original[0]}{$lcon}{$lread}{$name_sample} += $s;
                    }
                }
                $flag = 0;
            }
        }
    }
    foreach $key (sort keys %final) {
        $a = 0;
        foreach $key1 (sort {$a <=> $b} keys %{$final{$key}}) {
            foreach $key2 (keys %{$final{$key}{$key1}}) {
                $b = 0;
                foreach $key3 (sort keys %{$final{$key}{$key1}{$key2}}) {
                    $sum_single = 0;
                    $c = 0;
                    foreach $key4 (keys %{$final{$key}{$key1}{$key2}{$key3}}) {
                        $sum_single += $final{$key}{$key1}{$key2}{$key3}{$key4};
                    }
                    if (exists $def_onlysnps{$key}{$key1}) {
                        $percentage = ($sum_single / $def_onlysnps{$key}{$key1})*100;
                    } else {
                        $percentage = 0;
                    }
                    
                    if (($sum_single == 1) or ($percentage <= 1)) {
                        delete $final{$key}{$key1}{$key2}{$key3}{$key4};
                    } else {
                        foreach $start (sort {$a <=> $b} keys %{$definitive{$i}{$key}}) {
                            foreach $end (sort {$a <=> $b} keys %{$definitive{$i}{$key}{$start}}) {
                                $flag1 = 0;
                                $flagcodon = 0;
                                if (($start < $key1) && ($end > $key1)) {
                                    $frame = ($key1-$start)%3;
                                    if ($frame == 0) {
                                        $beg = $key1-$start;
                                    } elsif ($frame == 1) {
                                        $beg = $key1-$start-1;
                                    } elsif ($frame == 2) {
                                        $beg = $key1-$start-2;
                                    }
                                    $codon_con = substr $definitive{$i}{$key}{$start}{$end}, $beg, 3;
                                    @co = split //, $codon_con;
                                    $codon_snp = ();
                                    $a=0;
                                    foreach $c (@co) {
                                        if ($a == $frame) {
                                            $codon_snp = $codon_snp . $key3;
                                        } else {
                                            $codon_snp = $codon_snp . $c;
                                        }
                                        $a++;
                                    }
                                    if (exists $revert{$key}{$start}) {
                                        $codon_con =~ tr/CATG|catg/GTAC/;
                                        $codon_con = reverse($codon_con);
                                        $codon_snp =~ tr/CATG|catg/GTAC/;
                                        $codon_snp = reverse($codon_snp);
                                        $codon_con = uc $codon_con;
                                        $codon_snp = uc $codon_snp;
                                        $aa_con = translate($codon_con);
                                        $aa_snp = translate($codon_snp);
                                    } else {
                                        $codon_con = uc $codon_con;
                                        $codon_snp = uc $codon_snp;
                                        $aa_con = translate($codon_con);
                                        $aa_snp = translate($codon_snp);
                                    }
                                    $totcodonaa++;
                                    
                                    if ($aa_con eq $aa_snp) {
                                        $equal++;
                                        $flag1 = 1;
                                    } else {
                                        if ($codon_snp !~ /N/) {
                                            $flagcodon = 1;
                                            $flag1 = 2;
                                        }
                                    }
                                }
                            }
                        }
                        foreach $key5 (sort {$a <=> $b} keys %{$final{$key}{$key1}{$key2}{$key3}}) {
                            if ($flagcodon == 1) {
                                $differaa{$aa_con}{$aa_snp}{$key5} += $final{$key}{$key1}{$key2}{$key3}{$key5};
                                $differcod{$codon_con}{$codon_snp}{$key5} += $final{$key}{$key1}{$key2}{$key3}{$key5};
                            }
                            if ($a == 0) {
                                print $file "\n>$key\n";
                                $a++;
                            }
                            if ($b == 0) {
                                print $file "${key1}${key2}$def_onlysnps{$key}{$key1}:";
                                $b++;
                            }
                            if ($c == 0) {
                                if ($flag1 == 0) {
                                    print $file "${key3}:";
                                } elsif ($flag1 == 1) {
                                    print $file "${key3}\*:";
                                } else {
                                    print $file "${key3}\+:";
                                }
                                $c++;
                            }
                            print $file "$final{$key}{$key1}{$key2}{$key3}{$key5}\-${key5}\;";
                        }
                    }
                }
            }
        }
    }
    undef %final;
    $done = $round*$i;
    print "Processed $done contigs out of $count\n";
    close SNP;
}

close $file;

print "The number of substitution with synonimous SNPs is $equal out of $totcodonaa\n";

foreach $sample (sort {$a <=> $b} keys %sampleindex) {
    print $file1 "\;$sample";
    print $file2 "\;$sample";
    
}

foreach $key (sort keys %differaa) {
    foreach $key1 (sort keys %{$differaa{$key}}) {
        print $file1 "\n${key}\-${key1}";
        foreach $sample (sort {$a <=> $b} keys %sampleindex) {
            if (exists $differaa{$key}{$key1}{$sample}) {
                $fraction = ($differaa{$key}{$key1}{$sample}/$TOT{$sample})*100;
                print $file1 "\;";
                printf $file1 "%.4f", $fraction;
            } else {
                print $file1 "\;0";
            }
        }
    }
}

foreach $key (sort keys %differcod) {
    foreach $key1 (sort keys %{$differcod{$key}}) {
        print $file2 "\n${key}\-${key1}";
        foreach $sample (sort {$a <=> $b} keys %sampleindex) {
            if (exists $differcod{$key}{$key1}{$sample}) {
                $fraction = ($differcod{$key}{$key1}{$sample}/$TOT{$sample})*100;
                print $file2 "\;";
                printf $file2 "%.4f", $fraction;
            } else {
                print $file2 "\;0";
            }
        }
    }
}

close $file1;
close $file2;

sub translate {
    $lookup{AAC}="N";
    $lookup{AAG}="K";
    $lookup{ACA}="T";
    $lookup{ACC}="T";
    $lookup{ACG}="T";
    $lookup{AGA}="R";
    $lookup{AGC}="S";
    $lookup{AGG}="R";
    $lookup{AGT}="S";
    $lookup{ATA}="I";
    $lookup{ATC}="I";
    $lookup{CAC}="H";
    $lookup{CAG}="Q";
    $lookup{CCA}="P";
    $lookup{CCC}="P";
    $lookup{CCG}="P";
    $lookup{CGA}="R";
    $lookup{CGC}="R";
    $lookup{CGG}="R";
    $lookup{CTA}="L";
    $lookup{CTC}="L";
    $lookup{CTG}="L";
    $lookup{GAC}="D";
    $lookup{GAG}="E";
    $lookup{GCA}="A";
    $lookup{GCC}="A";
    $lookup{GCG}="A";
    $lookup{GGA}="G";
    $lookup{GGC}="G";
    $lookup{GGG}="G";
    $lookup{GTA}="V";
    $lookup{GTC}="V";
    $lookup{GTG}="V";
    $lookup{TAC}="Y";
    $lookup{TCA}="S";
    $lookup{TCC}="S";
    $lookup{TCG}="S";
    $lookup{TGC}="C";
    $lookup{TTA}="L";
    $lookup{TTC}="F";
    $lookup{TTG}="L";
    $lookup{AAA}="K";
    $lookup{AAT}="N";
    $lookup{ACT}="T";
    $lookup{ATG}="M";
    $lookup{ATT}="I";
    $lookup{CAA}="Q";
    $lookup{CAT}="H";
    $lookup{CCT}="P";
    $lookup{CGT}="R";
    $lookup{CTT}="L";
    $lookup{GAA}="E";
    $lookup{GAT}="D";
    $lookup{GCT}="A";
    $lookup{GGT}="G";
    $lookup{GTT}="V";
    $lookup{TAT}="Y";
    $lookup{TCT}="S";
    $lookup{TGG}="W";
    $lookup{TGT}="C";
    $lookup{TTT}="F";
    $lookup{TAG}="STOP";
    $lookup{TAA}="STOP";
    $lookup{TGA}="STOP";
    
    my $cod = shift;
    $aa = $lookup{$cod};
    return($aa);
}
