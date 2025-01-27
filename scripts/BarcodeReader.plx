#!/usr/bin/env perl


use strict;
use warnings;
use Data::Dumper;
use File::Find;
use File::Basename;
#use Math::Derivative qw(Derivative1 Derivative2);
use Cwd;
my %qtag_counts = ();

my $BARCODE_MOTIF = "CGA[ACTG]{3}C[ACTG]{4}AATTCGATGG";
my $MCOUNT_MOTIF = "[ATCG]{0,3}C[ACTG]{3}C[ACTG]{3}C[ACTG]{3}GCGCAACGCG";
my %BARCODE_COUNTS = ();


my %qtags = (
        "TCGGCTAGATGT" => "19",
        "AGGAACACCAAG" => "23",
	"TCGCCGAGCAGT" => "22",
        "CGAGCGCGAGGA" => "24",
        "TGGCGAATATGG" => "25",
        "TCTTCTACAACA" => "26",
        "AGCACGCCTTGT" => "27",
	"GCAACTTCTTCA" => "26_2",
	"AAGAAGTCCAAC" => "17",
);

sub log2 {
        my $n = shift;
        return log($n)/log(2);
    }

sub hamming{
     #String length is assumed to be equal
     my ($p,$o) = @_;
     return 0 if $p eq $o; 
     my $len = length ($p);
     my $num_mismatch = 0;
     for (my $i=0; $i<$len; $i++){
     	++$num_mismatch if substr($p, $i, 1) ne substr($o, $i, 1);
     }

     return $num_mismatch;
}

my %FILES = get_files();

sub get_files {
	my %hash = ();
	my @files = glob("*.fastq.gz");
	foreach my $file (@files) {
		my @array = split/\_/, $file;
		my $string = join("_", @array[0,1]);
		push @{$hash{$string}}, $file;
	}
	return %hash;
}


foreach my $reads (keys %FILES) {
	#next unless scalar(@{$FILES{$reads}}) == '2';
	print "Processing $reads\n";
	my %UNKNOWN = ();
	my %data = ();
	my $data_file = $reads . "_" . 'reads.csv';
	my $chimera_file = $reads . "_" . 'chimera_data.txt';
	#my $threshold_file = $reads . "_" . 'threshold_data.txt';
	my $mcount_barcode_file = $reads . "_" . 'mcount_barcode_data.txt';
	my $summary_file = $reads . "_" . 'summary_data.txt';
	my $novel_qtags = $reads . "_" . 'novel_qtags_data.txt';
	my $umi_distribution = $reads . "_" . "umi_distribution.txt";

	open my($of1), ">$chimera_file";
	#open my($of2), ">$threshold_file";
	open my($of3), ">$mcount_barcode_file";
	open my($of4), ">$summary_file";
	open my($of9), ">$novel_qtags";
	open my($of5), ">$umi_distribution";
	open my($of6), ">$data_file";
	print $of1 join("\t", 'run', 'total_reads', 'match_all_features_counts'), "\n";
	#print $of2 join("\t", 'run', 'qtag', 'index', 'barcode', 'counts', 'norm', 'dydx', 'dy2dx2'), "\n";
	print $of3 join("\t", 'id', 'mcount_barcode', 'qtag', 'count'), "\n";
	print $of4 join("\t", 'umi', 'barcode', 'qtag', 'counts'), "\n";
 
	open my($fh1), "seqtk seq -A ${$FILES{$reads}}[0] | ";

	my $counter = '0';
	my %feature_counts = ();
	my %summary = ();
	my %chimera = ();
	my %umi = ();
	my %alt = ();
	while(defined (my $line1 = <$fh1>) ) {
		next if $line1 =~ /^\>/;
		chomp $line1;
		#chomp $line2;
		my @tags = ();
		my @features = ();
		my $qtag = 'undef';
		$counter++;
		if ( $line1 =~ /($MCOUNT_MOTIF)([ATCG]+)($BARCODE_MOTIF)([ATCG]+TGGTGTTCAAGCTT)([ATCG]{12})/) {
			my @test = bin_qtags($5, keys %qtags);
			#next unless scalar(@test) == '1';
			if (scalar(@test) == '1') {
				
				push @tags, $1;
				push @features, 'mc';
				my @bc = split('', $3);
                        	my $string = join('', @bc[3,4,5,7,8,9,10]);
                        	push @tags, $string;
                        	push @features, 'bc';
				push @tags, $test[0];
                        	push @features, 'qtag';
			}else {
				$UNKNOWN{$5}++;
				push @tags, 'MISSING';
                        	push @features, 'na';
                        	push @tags, 'MISSING';
                        	push @features, 'na';
                        	push @tags, 'MISSING';
                        	push @features, 'na';
                        	push @tags, 'MISSING';
                        	push @features, 'na'; 
			}
		}else {
			push @tags, 'MISSING';
			push @features, 'na';
                    	push @tags, 'MISSING';
                        push @features, 'na';
			push @tags, 'MISSING';
                        push @features, 'na';	
                        push @tags, 'MISSING';
                        push @features, 'na';
		}
	

		if (defined $qtags{$tags[$#tags]} ) {
			$qtag = $qtags{$tags[$#tags]};
		}

		if ( grep { $_ eq 'MISSING'} @tags ) {
			$summary{'hits'}++;
		}
		$summary{'total'}++;
		my $string = join(",", @features);
		$feature_counts{$string}++;
		next if grep { $_ eq 'MISSING'} @tags;
		$BARCODE_COUNTS{$tags[1]}{$reads}++;
		my $umi_bc = $tags[0] . '_' . $tags[1];
		$umi{$tags[0]}{$tags[1]}++;
		$alt{$umi_bc}{$qtag}++;## change name of variable.
		print $of6 join("\t", $reads, $qtag, $tags[1], $tags[0], $line1), "\n";
		$data{$qtag}{$tags[1]}{$tags[0]}++;
		$chimera{$tags[1]}{$qtags{$tags[2]}}++;
	}

	my %COUNTS = ();
	foreach my $umi ( keys %umi) {
		foreach my $barcode ( keys %{$umi{$umi}} ) {
			$COUNTS{$umi{$umi}{$barcode}}++;
		}
	}

	foreach my $umi (sort {$a <=> $b }keys %COUNTS) {
		print $of5 join("\t", $reads, $umi, $COUNTS{$umi}), "\n";
		}
	

	foreach my $item (keys %UNKNOWN){
		print $of9 join("\t", $reads, $item, $UNKNOWN{$item}), "\n";
	}

	foreach my $item (keys %feature_counts) {
		my @array = split/\,/, $item;
		my $percent = '100' * ($feature_counts{$item} / $counter);
                my $rounded = sprintf("%.4f", $percent);
		print $of4 join("\t", $reads, @array, $feature_counts{$item}, $counter, $rounded), "\n";
	}

	foreach my $bc (keys %chimera) {
		my @q = keys %{$chimera{$bc}};
		next if scalar(@q) == '1';
		for my $item (@q) {
       			print $of1 join("\t", $reads, $bc, $item, $chimera{$bc}{$item}), "\n";
        	}
	}
	my $cnt = '0';
	foreach my $tag (keys %alt) {
		next unless scalar(keys %{$alt{$tag}})  > 1;
		my @stuff = keys %{$alt{$tag}};
		$cnt++;
		my $flag = '0';
		foreach my $bc ( keys %{$alt{$tag}} ) {
			for my $item (@stuff) { 
				if ($alt{$tag}{$item} >= '10') {
					$flag++;
				}
	
			}
			if ($flag > '0' ) {
                		print $of3 join("\t", $cnt, $tag, $bc, $alt{$tag}{$bc}), "\n";
        		}

		}
	}
	###umi barcode filter
	foreach my $qtag ( keys %data) {
		foreach my $barcode (keys %{$data{$qtag}}) {
			foreach my $umi (keys %{$data{$qtag}{$barcode}} ) {
				if (defined $umi{$umi}{$barcode}) {
					if ($umi{$umi}{$barcode} > '1') {
					#	delete $data{$qtag}{$barcode}{$umi};
					}
				}
			}
		}
	}		

}

#foreach my $barcode ( keys %BARCODE_COUNTS) {
#        foreach my $read ( keys %{$BARCODE_COUNTS{$barcode}}) {
#                print join("\t", $barcode, $read, $BARCODE_COUNTS{$barcode}{$read}), "\n";
#        }
#}

sub get_mcounts {
	my %data = @_;
      	my %clusters = ();
      	foreach my $qtag (keys %data) {
        	foreach my $barcode (keys %{$data{$qtag}}) {
                	$clusters{$qtag}{$barcode} = scalar(keys %{$data{$qtag}{$barcode}});
                }
        }
 	return %clusters;
}

sub return_list {
	my @array = @_;
	my @list = ();
	for my $item (@array) {
		push @list, $item->[7];
	}
	@list = sort {$a <=> $b} @list;
	return @list
}
		

sub remove_barcodes {
	my $run = shift;
	my @barcodes = @_;
	my %cache = ();
	my $out = $run . '.chimera_one_off.txt';
	open my($df), ">$out";
	my @barcodes_filtered = ();
	for (my $i = 0; $i <= $#barcodes; $i++) {
		next if $cache{$i};
		#print $df join("\t", $run, $barcodes[$i][1], $barcodes[$i][0], $barcodes[$j][1], $barcodes[$j][0], '0'), "\n";
                for (my $j = $i +1; $j <= $#barcodes; $j++) {
                	next if $cache{$j};
                        my $mm = hamming($barcodes[$i][1], $barcodes[$j][1]);
                        #if ($mm == '1' or $mm == '0') {
			if ($mm <= '1') {
			#print join("\t", @{$barcodes[$i]}). "\n";	
				print $df join("\t", $run, $barcodes[$i][1], $barcodes[$i][0], $barcodes[$j][1], $barcodes[$j][0], $mm), "\n";
				#$barcodes[$i][0] += $barcodes[$j][0];
				$cache{$j}++;
                        }
                }
        }

	for my $x (0.. $#barcodes) {
        	next if $cache{$x};
                push @barcodes_filtered, $barcodes[$x];
       	}
	close($out);
	return @barcodes_filtered;
}

##############################
sub remove_barcodes_edit {
        my $add = shift;
        my $run = shift;
        my @barcodes = @_;
        my %cache = ();
        my $out = $run . '.chimera_one_off.txt';
        open my($df), ">$out";
        my @barcodes_filtered = ();
        for (my $i = 0; $i <= $#barcodes; $i++) {
                next if $cache{$i};
                for (my $j = $i +1; $j <= $#barcodes; $j++) {
                        next if $cache{$j};
                        my $mm = hamming($barcodes[$i][1], $barcodes[$j][1]);
                        if ($mm <= '1') {
                                print $df join("\t", $run, $barcodes[$i][1], $barcodes[$i][0], $barcodes[$j][1], $barcodes[$j][0], $mm), "\n";
                                $cache{$j}++;
                                #next if $add eq 'FALSE';
                                #next if $mm == '1';#Adding one offs to parent
                                #$barcodes[$i][0] += $barcodes[$j][0];
                        }
                }
        }

        for my $x (0.. $#barcodes) {
                next if $cache{$x};
                push @barcodes_filtered, $barcodes[$x];
        }
        close($out);
        return @barcodes_filtered;
}




sub bin_qtags {
	my $test = shift;
        my @qtags = @_;
        my %cache = ();
        my @array = ();
        for (my $i = 0; $i <= $#qtags; $i++) {
        	my $mm = hamming($qtags[$i], $test);
                if ($mm <= '2') {
                       	push @array, $qtags[$i];
               }
        }
	return @array;
}




exit;
