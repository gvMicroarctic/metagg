#!/usr/bin/perl

#Mapping step: all the fastq file need to be in the same folder

%args = @ARGV;

$path = $args{-p}; #path to the files
$contig = $args{-c}; #assembly file
$thread = $args{-t}; #number of threads
$seed = $args{-s}; #seed length

if (!(defined $seed)) {
    $seed = 20;
}

if (!(defined $thread)) {
    $thread = 1;
}

if ($path != '') {
    @file_folder = `ls $path/*_R1.fastq`;
} else {
    @file_folder = `ls *_R1.fastq`;
}


foreach $file (@file_folder) {
    ($file1) = $file =~ /(.*)_R/;
    print "Running bwa mem -t $thread -k $seed ${contig} ${file1}_R1.fastq ${file1}_R2.fastq > ${file1}.sam\n";
    `bwa mem -t thread -k $seed $contig ${file1}_R1.fastq ${file1}_R2.fastq > ${file1}.sam`;
}




