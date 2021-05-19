#!perl 

# takes a list of SNPs and uses blast to generate a bedfile coordinates in Dsim
# get_Dmel_SNP_regions_in_Dsim_gz.pl

# perl ancestral_states_list.pl list
# list is a list of SNPs in the following format
# CHR    START    END    QVAL LTR    CHR OVERLAP    START GENE    END GENE    GENE ID    GENE NAME    NT OVERLAP
#chr2L    499719    499720    0.038673216    chr2L    476220    540560    FBgn0003963    ush    1
#chr2L    645185    645186    5.79E-07    chr2L    640013    714983    FBgn0284247    ds    1
#chr2L    838002    838003    0.026840972    chr2L    833584    851071    FBgn0020304    drongo    1

system("rm check_blast_results_Dsim");

if (@ARGV < 2){print "\nusage: zcat reffile.z.Z | perl get_Dmel_SNP_regions_in_Dsim_gz.pl list buffer\n\n";exit;}
# buffer is the number of bases to be added on either side of SNP
# i.e. enter 50 for SNP to be at 50th base position

my $list = shift(@ARGV); chomp $list;
open IN, $list or die "wrong format for infile\n";

my $buffer= shift(@ARGV); chomp $buffer;

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


#**************** read in snp file

my @snps = (); my @starts = (); my @arms = ();
while (my $line2 = <IN>) {
    chomp $line2;
    $line =~ s/\s+/\t/g; # replace all white spaces with a tab
    if ($line2 !~ /CHR/){ # skip header
        (my $arm, my $position) = split(/\t/, $line2); $arm=~s/chr//g; $position=$position+1;
        $start=$position-$buffer;
        push (@snps,"Dmel_".$arm."_".$position);
        push (@arms, $arm);
        push (@starts,$start);
    } # for each SNP
} # for all lines in list
close IN;

my @Dmel_seqs=();
open OUT, ">Dmel_snp_sequences.fa";
open OUT2, ">Dsim_snp_coordinates.bed";
for $j (0..scalar(@snps-1)){
    my $key = ">".$arms[$j]; chomp $key;
    # print "$key\t$starts[$j]\t$ends[$j]\n";
    my $seq = substr($data{$key}, $starts[$j], 100);
    print OUT ">$snps[$j]\n$seq\n";
    push (@Dmel_seqs,$seq);
    
    

    ####################################### blast subrountine
    open TEMP1, ">temp.fa"; print TEMP1 ">$snps[$j]\n$Dmel_seqs[$j]\n";
    open TEMP2, "./blast/blastn -query temp.fa -db  Dsim_w501_chromosome_arms_mtDNA_Pilon_w60.fasta -outfmt 6 -max_target_seqs 1 2>/dev/null |";
    system("./blast/blastn -query temp.fa -db  Dsim_w501_chromosome_arms_mtDNA_Pilon_w60.fasta -outfmt 6 -max_target_seqs 1 2>/dev/null >>check_blast_results_Dsim");
    my @lines = <TEMP2>; # print "@lines\n";
    my $tophit_location=();
    if (scalar(@lines)>0){ # if there is a blast hit
        my @positions = (); my $arm = (); my $direction =();
        for $i (0..(scalar(@lines)-1)){ # for all blast hits
            $lines[$i]=~ s/\s+/\t/g;
            my @working_line = split(/\s+/, $lines[$i]);
            my $distance = ();
            $arm = "$working_line[1]";
            if ($i==0){$tophit_location=$working_line[8];}
            else{$distance = abs($tophit_location - $working_line[8]);}
            #print "line: $i\tworking line[8]: $working_line[8]\ttop_hit: $tophit_location\tdistance: $distance\t";
            if ($distance<$buffer*2){ # only keep blast hit locations within $bufferX2 of top hit
                # print "recording positions\n";
                push (@positions,$working_line[8]);
                push (@positions,$working_line[9]);
                if ($i==0){if ($working_line[8]>$working_line[9]){$direction="-";}else{$direction="+";}}
                }
            # else{print "NOT recording positions\n";}
            } # for all blast hits
        # print "positions: ", $positions[0],"\t", $positions[1],"\n";
        my @positions_sorted = sort { $a <=> $b }(@positions);
        my $start = $positions_sorted[0];
        my $end = $positions_sorted[scalar(@positions)-1];
        # print "positions: ", $start,"\t", $end,"\n";
        $key = substr($key, 1);
        print "$snps[$j]\tDsim\t$arm\t$start\t$end\t$direction\t", $end-$start+1, "\n";
        print OUT2 "$snps[$j]\tDsim\t$arm\t$start\t$end\t$direction\t", $end-$start+1, "\n";
        $counter++;
        } # if there is a blast hit
    else{print "$snps[$j]\tDsim\tno_hit\n";; print OUT2 "$snps[$j]\tDsim\tno_hit\n";}
    close TEMP1;close TEMP2;
    ####################################### blast subrountine
    
}
close OUT;close OUT2;
system("rm temp.fa");



