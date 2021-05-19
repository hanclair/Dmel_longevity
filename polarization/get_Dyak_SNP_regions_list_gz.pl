#!perl
# reads in a bedfile of regions and creates a simple fasta of sequences

#  zcat Dyak_Tai18E2.fasta.gz.Z | perl get_Dyak_SNP_regions_list_gz.pl Dyak_snp_coordinates.bed

# list is a list of SNPs in the following format
# Dmel_2L_499720    no_hit
# Dmel_2L_645186    2L    699109    699208    +
# Dmel_2L_838003    2L    699109    699208    +

if (@ARGV < 1){print "\nusage: zcat reffile.z.Z | get_Dyak_SNP_regionslist_gz.pl list\n\n";exit;}

my $list = shift(@ARGV); chomp $list;
open IN, $list or die "wrong format for infile\n";

#**************** read in sequence & create output sequence
my %data = ();
my @sequence_names = ();
my @sequence_temp =();
my @data_temp =();

    while (<>) {
        my $line = $_;
        chomp $line;
        $line =~ s/\s+/\t/g; # replace all white spaces with a tab
        if ($line =~ />/)
            {
            (my $name, my $other) = split(/\t/, $line); chomp $name;
            push(@sequence_names, $name);
            $sequence_temporary = join ("", @data_temp);
            push(@sequence_temp,$sequence_temporary);
            @data_temp=();
            }
        else {push(@data_temp, $line);}
        }
    close REF;
    $sequence_temporary = join ("", @data_temp);
    push(@sequence_temp,$sequence_temporary);
    for $x (0..(scalar(@sequence_names)-1))    # put sequences into @data
        {
        $data{$sequence_names[$x]}=$sequence_temp[$x+1];
        }
my $seq_count = scalar (keys (%data));
# print "number of sequences: $seq_count\n";
#****************


#**************** read in bed file
my @snps = (); my @starts = (); my @ends = (); my @arms = ();
while (my $line2 = <IN>) {
    chomp $line2;
    $line =~ s/\s+/\t/g; # replace all white spaces with a tab
    if ($line2 !~ /no_hit/){ # skip no hits
        (my $snp, my $species, my $arm, my $start, my $end) = split(/\t/, $line2);
        push (@snps,$snp."_".$species."_".$arm."_".$start."_".$end);
        push (@arms, $arm);
        push (@starts,$start);
        push (@ends,$end);
    } # for each SNP
    else {
        push (@snps,"Dyak_no_hit");
        push (@arms, "no_hit");
        push (@starts,"no_hit");
        push (@ends,"no_hit");
    }
} # for all lines in list
close IN;

my @Dyak_seqs=();
open OUT, ">Dyak_snp_sequences.fa";
for $j (0..scalar(@snps-1)){
    if ($arms[$j] !~ /no_hit/){
        my $key = ">".$arms[$j]; chomp $key;
        # print "$key\t$starts[$j]\t$ends[$j]\n";
        my $seq = substr($data{$key}, $starts[$j], $ends[$j]-$starts[$j]);
        print ">$snps[$j]\n";
        print OUT ">$snps[$j]\n$seq\n";
        push (@Dyak_seqs,$seq);
        } # if there is a hit
    else{push (@Dyak_seqs,"no_hit"); print OUT ">$snps[$j]\nno_hit\n";}
    } # for all snp regions
   
close OUT;


