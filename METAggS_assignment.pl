#!/usr/bin/perl

open(LCA, "<$ARGV[0]") or die "Coudn't open the file $ARGV[0]: $!";

open(DIAMOND, "<$ARGV[1]") or die "Coudn't open the file $ARGV[1]: $!";

open(INDEX, "<$ARGV[2]") or die "Coudn't open the file $ARGV[2]: $!";

open(TOTREADS, "<tot_reads.txt") or die "Coudn't open the file tot_reads.txt: $!";

while (defined($input = <TOTREADS>)) {
    chomp($input);
    @reads = split /\t/, $input;
    $readtot{$reads[0]} = ($reads[1]/2);
}

while (defined($input = <LCA>)) {
    chomp($input);
    @taxa = split(/\t/, $input);
    if (!(defined $taxa[1])) {
        $offs{$taxa[0]}{'domain'} = 'Unclassified';
        $offs{$taxa[0]}{'phylum'} = 'Unclassified';
        $offs{$taxa[0]}{'class'} = 'Unclassified';
        $offs{$taxa[0]}{'order'} = 'Unclassified';
        $offs{$taxa[0]}{'family'} = 'Unclassified';
        $offs{$taxa[0]}{'genus'} = 'Unclassified';
        $offs{$taxa[0]}{'species'} = 'Unclassified';
    } else {
        $offs{$taxa[0]}{'domain'} = $taxa[1];
        $offs{$taxa[0]}{'phylum'} = $taxa[2];
        $offs{$taxa[0]}{'class'} = $taxa[3];
        $offs{$taxa[0]}{'order'} = $taxa[4];
        $offs{$taxa[0]}{'family'} = $taxa[5];
        $offs{$taxa[0]}{'genus'} = $taxa[6];
        $offs{$taxa[0]}{'species'} = $taxa[7];
    }
}

while (defined($input = <DIAMOND>)) {
    chomp ($input);
    @matches = split /\t/, $input;
    if (defined $offs{$matches[0]}) {
        if ($matches[6] < $matches[7]) {
            $ha{$matches[1]}{$matches[0]}{'S'} = $matches[6];
            $ha{$matches[1]}{$matches[0]}{'E'} = $matches[7];
            $def{$matches[0]}{$matches[1]}{'S'} = $matches[6];
            $def{$matches[0]}{$matches[1]}{'E'} = $matches[7];
        } else {
            $ha{$matches[1]}{$matches[0]}{'S'} = $matches[7];
            $ha{$matches[1]}{$matches[0]}{'E'} = $matches[6];
            $def{$matches[0]}{$matches[1]}{'S'} = $matches[7];
            $def{$matches[0]}{$matches[1]}{'E'} = $matches[6];
        }
    }
}

open(GO, "<../../Database/index.txt") or die "Coudn't open the file $ARGV[1]: $!";
while (defined($input = <GO>)) {
    chomp ($input);
    $input1 = substr $input, 1;
    if ($input =~ /^>/) {
        @go_accession = split /;/, $input1;
        foreach $g (@go_accession1) {
            $g1 = "UniRef100_" . $g;
            if (exists $ha{$g1}) {
                foreach $k0 (keys %prov) {
                    foreach $k1 (keys %{$prov{$k0}}) {
                        $go{$g1}{$k0}{$k1} = '';
                    }
                }
            }
        }
        undef %prov;
    } else {
        @go_code = split /\t/, $input;
        $prov{$go_code[0]}{$go_code[1]} = '';
        @go_accession1 = @go_accession;
    }
}

print "done1\n";

while (defined($input = <INDEX>)) {
    chomp($input);
    ($acc, $prot) = split /\t/, $input;
    $acc1 = "UniRef100_" . $acc;
    if (exists $ha{$acc1}) {
        $uni{$acc1} = $prot;
    }
}

@cat =('P','F','C');

@phyl = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species');

@file = `ls ../SvalbardSeqData/trimmed_reads/combined/SNPS/*_bwacut20.sam`;
foreach $fi (@file) {
    %match = ();
    chomp($fi);
    ($name_sample) = $fi =~ /SNPS\/(.*?)_/;
    push @name, $name_sample;
    print "Processing $fi\n";
    open(FILE_SAM, "<$fi") or die "Coudn't open the file $fi: $!";
    print "Creating ${name_sample}_bwacut20less1.sam\n";
    $readtrimmed = 0;
    while (defined($input = <FILE_SAM>)) {
        if ($input =~ /^\@SQ/) {
            ($contig_length) = $input =~ /SN:(.*?)\t/;
            ($contig_bp) = $input =~ /LN:(.*?)\n/;
            $length{$contig_length} = $contig_bp;
        } elsif ($input !~ /^\@/) {
            ($score) = $input =~ /AS:i:(.*?)\t/;
            if ($score > 40) {
                $readtrimmed++;
                ($read,$bitflag,$con_read,$start_al,$mapscore,$cigar,$second_cont,$second_start,$inc,$seq,$rest) = split(/\t/, $input,11);
                ($md) = $input =~ /MD:Z:(.*?)\t/;
                @cigar_snp = split /(?<=\M|I|D|S|H)/, $cigar;
                $match_tot = 0;
                print $filesam "$con_read\t$read\t$start_al\t";
                $sum = 0;
                foreach $c (@cigar_snp) {
                    if ($c =~ /M/) {
                        $match_tot++;
                        print $filesam "$c";
                        ($sum0) = $c =~ /(.*?)[A-Z]/;
                        $sum += $sum0;
                    } elsif ($c =~ /S/) {
                        ($soft) = $c =~ /(.*?)[A-Z]/;
                        if ($match_tot == 0) {
                            $seq = substr($seq, $soft);
                        } else {
                            $seq = substr($seq, 0, -($soft));
                        }
                    } elsif ($c =~ /I/){
                        print $filesam "$c";
                        ($sum0) = $c =~ /(.*?)[A-Z]/;
                        $sum += $sum0;
                    } elsif ($c =~ /D/) {
                        print $filesam "$c";
                    }
                }
                print $filesam "\t$md\t$seq\n";
                $taxa{$con_read}++;
                $end_al = $start_al+$sum;
                $inside = 0;
                $gene_tot_length = 0;
                if (exists $def{$con_read}) {
                    foreach $key (keys %{$def{$con_read}}) {
                        $gene_length = $def{$con_read}{$key}{'E'} - $def{$con_read}{$key}{'S'};
                        $gene_tot_length += $gene_length;
                        if ((($start_al > $def{$con_read}{$key}{'S'}) && ($start_al < $def{$con_read}{$key}{'E'})) or (($end_al > $def{$con_read}{$key}{'S'}) && ($end_al < $def{$con_read}{$key}{'E'}))) {
                            $inside = 1;
                            $proportion = (1/$gene_length)*1000;
                            $finalgene{$uni{$key}}{$name_sample} += $proportion;
                            $finalgene0{$name_sample} += $proportion;
                            $combined{$con_read}{$uni{$key}} += $proportion;
                            if (!(defined $go{$key})) {
                                $combined_go{$con_read}{'P'}{'Unknown Biological Process'}+=$proportion;
                                $combined_go{$con_read}{'C'}{'Unknown Cellular Component'}+=$proportion;
                                $finalgene_go0{$name_sample} += $proportion;
                                $finalgene_go{'F'}{'Unknown Function'}{$name_sample} += $proportion;
                                $finalgene_go{'P'}{'Unknown Biological Process'}{$name_sample} += $proportion;
                                $finalgene_go{'C'}{'Unknown Cellular Component'}{$name_sample} += $proportion;
                            } else {
                                foreach $c (@cat) {
                                    if (!(defined $go{$key}{$c})) {
                                        if ($c eq 'F') {
                                            $def = 'Unknown Function';
                                        } elsif ($c eq 'P') {
                                            $def = 'Unknown Biological Process';
                                        } elsif ($c eq 'C') {
                                            $def = 'Unknown Cellular Component';
                                        }
                                        $combined_go{$con_read}{$c}{$def}+=$proportion;
                                        $finalgene_go0{$name_sample} += $proportion;
                                        $finalgene_go{$c}{$def}{$name_sample} += $proportion;
                                    } else {
                                        $num = keys %{$go{$key}{$c}};
                                        $proportion_go = $proportion/$num;
                                        foreach $k1 (keys %{$go{$key}{$c}}) {
                                            $combined_go{$con_read}{$c}{$k1} += $proportion_go;
                                            $finalgene_go0{$name_sample} += $proportion_go;
                                            $finalgene_go{$c}{$k1}{$name_sample} += $proportion_go;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    $noncoding_length = $length{$con_read}-$gene_tot_length;
                } else {
                    $noncoding_length = $length{$con_read};
                }
                if ($inside == 0) {
                    $proportion = (1/$noncoding_length)*1000;
                    $finalgene{'Non_coding'}{$name_sample} += $proportion;
                    $finalgene0{$name_sample} += $proportion;
                    $combined_go{$con_read}{'F'}{'Non_coding'} += $proportion;
                    $combined_go{$con_read}{'P'}{'Non_coding'} += $proportion;
                    $combined_go{$con_read}{'C'}{'Non_coding'} += $proportion;
                    $finalgene_go0{$name_sample} += $proportion;
                    $finalgene_go{'F'}{'Non_coding'}{$name_sample} += $proportion;
                    $finalgene_go{'P'}{'Non_coding'}{$name_sample} += $proportion;
                    $finalgene_go{'C'}{'Non_coding'}{$name_sample} += $proportion;
                }
            }
        }
    }
    $fraction = $readtrimmed/$readtot{$name_sample};
    print "ok $fraction and $readtrimmed and $readtot{$name_sample} and $finalgene0{$name_sample}\n";
    $finalgene0{$name_sample} = $finalgene0{$name_sample}/$fraction;
    print "ok1 $finalgene0{$name_sample} and $fraction\n";
    $flag = 0;
    foreach $key (keys %taxa) {
        foreach $ph (@phyl) {
            $proportion = ($taxa{$key}/$length{$key})*1000;
            $finaltaxa{$ph}{$offs{$key}{$ph}}{$name_sample} += $proportion;
            $finaltaxa0{$name_sample} += $proportion;
            foreach $key1 (keys %{$combined{$key}}) {
                $finalcombined{$ph}{$offs{$key}{$ph}}{$key1}{$name_sample} += $combined{$key}{$key1};
            }
            foreach $key1 (keys %{$combined_go{$key}}) {
                foreach $key2 (keys %{$combined_go{$key}{$key1}}) {
                    $finalcombined_go{$ph}{$offs{$key}{$ph}}{$key1}{$key2}{$name_sample} += $combined_go{$key}{$key1}{$key2};
                }
            }
        }
    }
    $finaltaxa0{$name_sample} = $finaltaxa0{$name_sample}/$fraction;
    close FILE_SAM;
    close $filesam;
    undef %taxa;
    undef %combined;
    undef %combined_go;
}

foreach $ph (@phyl) {
    open($file, ">TAXA_${ph}.txt") or die "Coudn't open the file TAXA_${ph}.txt: $!";
    open($file1, ">TAXA_GENE_${ph}.txt") or die "Coudn't open the file TAXA_GENE_${ph}.txt: $!";
    open($filep, ">TAXA_GOP_GENE_${ph}.txt") or die "Coudn't open the file TAXA_GOP_GENE_${ph}.txt: $!";
    open($filef, ">TAXA_GOF_GENE_${ph}.txt") or die "Coudn't open the file TAXA_GOS_GENE_${ph}.txt: $!";
    open($filec, ">TAXA_GOC_GENE_${ph}.txt") or die "Coudn't open the file TAXA_GOC_GENE_${ph}.txt: $!";
    print $file "\t", join("\t", @name);
    print $file1 "\t", join("\t", @name);
    print $filep "\t", join("\t", @name);
    print $filef "\t", join("\t", @name);
    print $filec "\t", join("\t", @name);
    foreach $key (sort keys %{$finaltaxa{$ph}}) {
        print $file "\n$key";
        foreach $na (@name) {
            if (exists $finaltaxa{$ph}{$key}{$na}) {
                $p = sprintf "%.7f", $finaltaxa{$ph}{$key}{$na};
                print $file "\t$p";
            } else {
                print $file "\t0";
            }
        }
        foreach $key1 (sort keys %{$finalcombined{$ph}{$key}}) {
            print $file1 "\n$key\t$key1";
            foreach $na (@name) {
                if (exists $finalcombined{$ph}{$key}{$key1}{$na}) {
                    $p = sprintf "%.7f", $finalcombined{$ph}{$key}{$key1}{$na};
                    print $file1 "\t$p";
                } else {
                    print $file1 "\t0";
                }
            }
        }
        foreach $go0 (sort keys %{$finalcombined_go{$ph}{$key}}) {
            if ($go0 eq 'F') {
                foreach $ff (sort keys %{$finalcombined_go{$ph}{$key}{$go0}}) {
                    print $filef "\n$key\t$ff";
                    foreach $na (@name) {
                        if (exists $finalcombined_go{$ph}{$key}{$go0}{$ff}{$na}) {
                            $p = sprintf "%.7f", $finalcombined_go{$ph}{$key}{$go0}{$ff}{$na};
                            print $filef "\t$p";
                        } else {
                            print $filef "\t0";
                        }
                    }
                }
            } elsif ($go0 eq 'P') {
                foreach $pp (sort keys %{$finalcombined_go{$ph}{$key}{$go0}}) {
                    print $filep "\n$key\t$pp";
                    foreach $na (@name) {
                        if (exists $finalcombined_go{$ph}{$key}{$go0}{$pp}{$na}) {
                            $p = sprintf "%.7f", $finalcombined_go{$ph}{$key}{$go0}{$pp}{$na};
                            print $filep "\t$p";
                        } else {
                            print $filep "\t0";
                        }
                    }
                }
            } elsif ($go0 eq 'C') {
                foreach $cc (sort keys %{$finalcombined_go{$ph}{$key}{$go0}}) {
                    print $filec "\n$key\t$cc";
                    foreach $na (@name) {
                        if (exists $finalcombined_go{$ph}{$key}{$go0}{$cc}{$na}) {
                            $p = sprintf "%.7f", $finalcombined_go{$ph}{$key}{$go0}{$cc}{$na};
                            print $filec "\t$p";
                        } else {
                            print $filec "\t0";
                        }
                    }
                }
            }
        }
    }
    close $file;
    close $file1;
    close $filep;
    close $filef;
    close $filec;
}

open($file2, ">GENE_TOT.txt") or die "Coudn't open the file GENE_TOT.txt: $!";

print $file2 "gene;", join(";", @name);
foreach $key (sort keys %finalgene) {
    print $file2 "\n$key";
    foreach $na (@name) {
        if (exists $finalgene{$key}{$na}) {
            $p = sprintf "%.7f", $finalgene{$key}{$na};
            print $file2 "\t$p";
        } else {
            print $file2 "\t0";
        }
    }
}

open($filept, ">GENE_TOT_GOP.txt") or die "Coudn't open the file GENE_TOT_GOP.txt: $!";
open($fileft, ">GENE_TOT_GOS.txt") or die "Coudn't open the file GENE_TOT_GOS.txt: $!";
open($filect, ">GENE_TOT_GOC.txt") or die "Coudn't open the file GENE_TOT_GOC.txt: $!";

print $filept "gene;", join(";", @name);
print $fileft "gene;", join(";", @name);
print $filect "gene;", join(";", @name);
foreach $key (sort keys %finalgene_go) {
    if ($key eq 'P') {
        foreach $key1 (sort keys %{$finalgene_go{$key}}) {
            print $filept "\n$key1";
            foreach $na (@name) {
                if (exists $finalgene_go{$key}{$key1}{$na}) {
                    $p = sprintf "%.7f", $finalgene_go{$key}{$key1}{$na};
                    print $filept "\t$p";
                } else {
                    print $filept "\t0";
                }
            }
        }
    } elsif ($key eq 'F') {
        foreach $key1 (sort keys %{$finalgene_go{$key}}) {
            print $fileft "\n$key1";
            foreach $na (@name) {
                if (exists $finalgene_go{$key}{$key1}{$na}) {
                    $p = sprintf "%.7f", $finalgene_go{$key}{$key1}{$na};
                    print $fileft "\t$p";
                } else {
                    print $fileft "\t0";
                }
            }
        }
    } elsif ($key eq 'C') {
        foreach $key1 (sort keys %{$finalgene_go{$key}}) {
            print $filect "\n$key1";
            foreach $na (@name) {
                if (exists $finalgene_go{$key}{$key1}{$na}) {
                    $p = sprintf "%.7f", $finalgene_go{$key}{$key1}{$na};
                    print $filect "\t$p";
                } else {
                    print $filect "\t0";
                }
            }
        }
    }
}

close $fileft;
close $filect;
close $filept;
