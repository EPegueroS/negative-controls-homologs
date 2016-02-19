#!/usr/bin/perl

use strict;
use warnings;
no warnings ('uninitialized', 'substr');

my ($inputDir_BLAST,@inputFiles_BLAST,$inputDir_mast_clusters_fungi, @inputFiles_mast_clusters_fungi,$inputDir_mast_fungi_vs_mosca, @inputFiles_mast_fungi_vs_mosca);
my $n_clusters;


#Directory of fly vs. fungi BLAST
#/home/epeguero/motif_find_meme/MEME_8_8_opt/MOSCA/BLAST_mosca
$inputDir_BLAST = $ARGV[0];
chomp $inputDir_BLAST;
if(substr($inputDir_BLAST,length($inputDir_BLAST)-1,1) eq '/'){
	chop $inputDir_BLAST;
} #if

#Directory of mast performed on fungi clusters
#/home/epeguero/motif_find_meme/MEME_8_8_opt/MAST
$inputDir_mast_clusters_fungi = $ARGV[1];
chomp $inputDir_mast_clusters_fungi;
if(substr($inputDir_mast_clusters_fungi,length($inputDir_mast_clusters_fungi)-1,1) eq '/'){
        chop $inputDir_mast_clusters_fungi;
} #if

#Directory of3 mast fungi vs fly
#/home/epeguero/motif_find_meme/MEME_8_8_opt/MOSCA
$inputDir_mast_fungi_vs_mosca = $ARGV[2];
chomp $inputDir_mast_fungi_vs_mosca;
if(substr($inputDir_mast_fungi_vs_mosca,length($inputDir_mast_fungi_vs_mosca)-1,1) eq '/'){
        chop $inputDir_mast_fungi_vs_mosca;
} #if

$n_clusters = $ARGV[5];
#
my $org_file = $ARGV[3];

my $ptt_fly = $ARGV[4];


##open directory of blasts fly vs. fungi and read files $ARGV[0];
opendir(DIR1,"$inputDir_BLAST") or die "$!";
@inputFiles_BLAST = sort grep {/mosca_vs_todos.BLAST.out$/} readdir(DIR1);
#print "Input files BLAST fly vs. fungi:	@inputFiles_BLAST\n";
closedir(DIR1);
my $inputFile_BLAST = $inputFiles_BLAST[0];

my $path_blast = "$inputDir_BLAST"."/"."$inputFile_BLAST";
#print "$path_blast\n";
open (IN1,"$path_blast") or die "$!";

###Identify pairs of homologs
my %id_pairs_fungi_fly_blast;
while (my $a = <IN1>) {
	chomp $a;
	my ($fly, $fungi) = split(/\t/,$a);
	my (undef,$id_fungi) = split(/\|/,$fungi);
	#print "$id_fungi\n";
	my (undef,$id_fly) = split(/\|/,$fly);
	#print "$fly\n";
	#print "$id_fly\n";
	my $pair_fungi_fly = $id_fungi."\t".$id_fly;
	#print "$pair_fungi_fly\n";
	$id_pairs_fungi_fly_blast{$pair_fungi_fly} = $pair_fungi_fly;
} #while
close IN1;

##open directory of mast performed on fungi clusters and read files $ARGV[1];
opendir(DIR2,"$inputDir_mast_clusters_fungi") or die "$!";
#print "n_clusters: $n_clusters\n";

while (my $file =  readdir(DIR2)) {
	@inputFiles_mast_clusters_fungi = sort grep {/mast_out_/} readdir(DIR2);
	#print "@inputFiles_mast_clusters_fungi\n";
	#print "Input files mast on clusters: $i	$j	 $file_test\n"
	#my $i = "a";
	#my $j = "b";
	#if (length($file) == 17 )  {
	#	$i = substr($file,9,2);
	#	#print "$i\n";
	#	$j = substr($file,15,2);
	#	#print "$j\n";
	#} #if
	#if (length($file) == 15 ) {
	#	$i = substr($file,9,1);
        #        #print "$i\n";
        #        $j = substr($file,14,1);
        #        #print "$j\n";
	#} #if

	#if (($i eq $j) && $file =~ /^mast_out/) {
	#	#print "Input files mast on clusters:	$file\n";
	#	push(@inputFiles_mast_clusters_fungi,$file);
	#	#print "@inputFiles_mast_clusters_fungi\n";
	#} #if
} #while

closedir(DIR2);

open (IN2,"$org_file") or die "$!";
my %gi_gen_description;
while (my $b = <IN2>) {
        chomp $b;
        #print "$b\n";
        my ($gi, $gen,$description) = split(/\t/,$b);
        $gi_gen_description{$gi} = $gen."\t".$description;
} #while
close IN2;

open (IN_fly_description,"$ptt_fly") or die "$!";
my %gi_gen_description_fly;
while (my $f = <IN_fly_description>) {
	chomp $f;
	if ($f =~ /^[0-9]/) {
		my (undef,undef,undef,$gi_f ,$gen_f ,$code,undef,undef,$description_f) = split(/\t/,$f);
        	$gi_gen_description_fly{$gi_f} = $gen_f."\t".$description_f;
	} #if
} #while

close IN_fly_description;

##open directory of mast performed on fungi vs fly $ARGV[2]
opendir(DIR3,"$inputDir_mast_fungi_vs_mosca") or die "$!";
@inputFiles_mast_fungi_vs_mosca = sort grep {/_vs_mosca$/} readdir(DIR3);
#print "@inputFiles_mast_fungi_vs_mosca\n";
closedir(DIR3);

####open each of the files corresponding to mast 
###(only those corresponding to the same clusters)
my %pair_description_out;
foreach my $directory (@inputFiles_mast_clusters_fungi) {
	my %id_found;
	my %id_found_m;
	chomp $directory;
	#print "$directory\n";
	my %result_out;
	my $string = substr($directory,8,6);
	#print "$string\n";
	my @matches = grep(/$string/,@inputFiles_mast_fungi_vs_mosca);
	#print "Matches:	@matches\n";
	#my $match = $matches[0];
	#my $path_mast_fly =$inputDir_mast_fungi_vs_mosca."/".$match."/"."mast.txt";
	#print "$path_mast_fly\n";
	foreach my $match (@matches) {
		#print "$match\n";
	 
		my $path_mast_fly =$inputDir_mast_fungi_vs_mosca."/".$match."/"."mast.txt";
	
		open (INfly,"$path_mast_fly") or die "$!";
        	while (my $m = <INfly>) {
                	chomp $m;
                	#print "$m\n";
                	if ($m =~ /^\d{6,}\|/) {
                        	my ($id_m,undef) = split(/\|/,$m);
                        	$id_found_m{$id_m} = $directory."\t".$match;
                	} #if
        	} #while	
		close INfly;

		my $path_mast_clusters = $inputDir_mast_clusters_fungi."/".$directory;
		#print "$path_mast_clusters\n";
		my $file_mast = $path_mast_clusters."/"."mast.txt"; 
		open (IN3,"$file_mast") or die "$!";
		while (my $c = <IN3>) {
			chomp $c;
			#print "$c\n";
			if ($c =~ /^\d{6,}\|/) {
				my ($id,undef) = split(/\|/,$c);
				$id_found{$id} = $id;
			} #if
		} #while
		close IN3;
		my %pairs;
		my %results;
		while (my($llave,$valor) = each(%id_found)) {
			#print "$llave\n";
			chomp $llave;
			while (my($key,$value) = each(%id_found_m)) {
				chomp $key;
				my $pair = $llave."\t".$key;
				$pairs{$pair} = $pair;
				#print "$pair\n";
			} #while	
			while (my($key_p,$value_p) = each(%pairs)) {
				if (defined($id_pairs_fungi_fly_blast{$key_p})) { 
					#print "$key_p	$directory	$match\n";
					my $key_out = $key_p."\t".$directory."\t".$match;
					$pair_description_out{$key_out} = $directory."\t".$match; 
				} #if
			} #while

		} #while
		#while (my ($k,$v) = each(%results)) {
        	#	$result_out{$k} = $v;
			#print "$k	$v\n";
        	#} # while	

	} #foreach
	#while (my($key_out,$value_out) = each(%result_out)) {
		#print "$key_out	$value_out\n";
	#	$pair_description_out{$key_out} = $value_out;
	#} #while
} #foreach


while (my($key,$value) = each(%pair_description_out)) {
	#print "$key	$value\n";
	my ($identifier_fungi, $identifier_fly) = split(/\t/,$key);
	print "$key	$identifier_fungi\t$gi_gen_description{$identifier_fungi}\t$identifier_fly\t$gi_gen_description_fly{$identifier_fly}\n";

} #while

exit;
