#!perl 

# snp_ancestral_state.pl

# reads in multi fasta for 3 seperate species e.g.
# Dmel_snp_sequences.fa
# Dsim_snp_sequences.fa
# Dsim_snp_sequences.fa
# aligns using muscle, QC and calls ancestral state of focal SNP
# currently set up for sequences for which focal base is position 50 in Dmel sequence, also hard-coded for Dmel, Dsim and Dyak

# requirements: /MUSCLE/muscle3.8.31_i86darwin64

if (@ARGV < 3){print "\nusage: snp_ancestral_state.pl species1.fa species2.fa species3.fa \nNote: input is a fasta multialignments for three species in the following format for each species:\n\n>Dmel_region1\nGACTGACTGACTGACTGACTGGC\n>Dmel_region2\nGACTACATCATACTCATCCATA\n>Dmel_region3\nATTATGTAGATAGATAGGATGAT\n>Dmel_region5\nGATACACCAACCCCATAGGATAG\n>Dmel_region6\nCAGACCCAAATGGGATTAGATAG\n\n";exit;}

#************************ READ IN DATA HERE **************************************
# species 1
my $infile1 = shift(@ARGV); open IN1, $infile1 or die "cannot find fasta file for species 1\n";
my @seq_names1 = (); my @sequences1 = ();
while (my $line1 = <IN1>) {chomp $line1; push (@seq_names1, $line1); $line1 = <IN1>; chomp $line1; push (@sequences1, $line1);} close IN1;
# species 2
my $infile2 = shift(@ARGV); open IN2, $infile2 or die "cannot find fasta file for species 2\n";
my @seq_names2 = (); my @sequences2 = ();
while (my $line2 = <IN2>) {chomp $line2; push (@seq_names2, $line2); $line2 = <IN2>; chomp $line2; push (@sequences2, $line2);} close IN2;
# species 3
my $infile3 = shift(@ARGV); open IN3, $infile3 or die "cannot find fasta file for species 3\n";
my @seq_names3 = (); my @sequences3 = ();
while (my $line3 = <IN3>) {chomp $line3; push (@seq_names3, $line3); $line3 = <IN3>; chomp $line3; push (@sequences3, $line3);} close IN3;
# print "number of sequences: ", scalar(@seq_names1), "\n";

open OUT2, ">snp_ancestral_state_log";
open OUT3, ">snp_ancestral_state_results";
print OUT3 "species\tarm\tsnp_pos\tstates_MSY\tANC\n";

# create MUSCLE alignment and read in alignment **************************************
for $j (0..scalar(@seq_names1)-1){
        open OUT1, ">temp.fa";
        print OUT1 "$seq_names1[$j]\n$sequences1[$j]\n$seq_names2[$j]\n$sequences2[$j]\n$seq_names3[$j]\n$sequences3[$j]\n";
        system("./muscle/muscle3.8.31_i86darwin64 -in temp.fa -out temp2.fa -diags 2>/dev/null");
        
        open ALN, "temp2.fa" or die "can't find muscle output\n";
        my @sequence_names = ();
        my @sequence_temp =();
        my @data_temp =();
        #**************** read in fasta
        while (my $line = <ALN>) {
            chomp $line;
            # print "$line\n\n";
            # $line =~ s/\s+/\t/g; # remove all white spaces
            if ($line =~ />/){
                my $name = $line; chomp $name;
                push(@sequence_names, $name);
                $sequence_temporary = join ("", @data_temp);
                push(@sequence_temp,$sequence_temporary);
                @data_temp=();
                } # if the first line is a name
                else {push(@data_temp, $line);}
                } # while there are lines
        close ALN;
        $sequence_temporary = join ("", @data_temp);
        push(@sequence_temp,$sequence_temporary);
        
        # put sequences into data_hash but reorder so order is mel, sim and yak
        my %data_hash = ();
        for $x (0..(scalar(@sequence_names)-1)){$data_hash{$sequence_names[$x]}=$sequence_temp[$x+1];}
        my @data = (); my @sequence_names_new = ();
        foreach $key (sort {lc($a) cmp lc($b)} keys(%data_hash)){
                if ($key=~ /sim/){$data[1]=$data_hash{$key};$sequence_names_new[1]=$key;}
                elsif ($key=~ /yak/){$data[2]=$data_hash{$key};$sequence_names_new[2]=$key;}
                else {$data[0]=$data_hash{$key};$sequence_names_new[0]=$key;}
                } # for each key
        # print "$sequence_names_new[0]\n$data[0]\n$sequence_names_new[1]\n$data[1]\n$sequence_names_new[2]\n$data[2]\n\n";
        my @parsed_snp_ID = split(/\_/,$sequence_names_new[0]);
        # print substr($parsed_snp_ID[0],1),"\t$parsed_snp_ID[1]\t$parsed_snp_ID[2]\n\n";
    
        if (($sequences2[$j] !~ /no_hit/) && ($sequences3[$j] !~ /no_hit/)){ # if there was no blast hit in either Dsim or Dyak
        
            # remove all insertions relative to Dmel sequence
            my @data_new=();
            for $k (0..length($data[0])){
                if (substr($data[0],$k,1) ne "-"){
                    $data_new[0] = $data_new[0].substr($data[0],$k,1);
                    $data_new[1] = $data_new[1].substr($data[1],$k,1);
                    $data_new[2] = $data_new[2].substr($data[2],$k,1);
                    }
                } # for entire length of sequence
             
            print OUT2 "Dmel insertions removed:\n";
            print OUT2 "$sequence_names_new[0]\n", substr($data_new[0],0,50),"\t",substr($data_new[0],50,50), "\n$sequence_names_new[1]\n", substr($data_new[1],0,50),"\t",substr($data_new[1],50,50), "\n$sequence_names_new[2]\n", substr($data_new[2],0,50),"\t",substr($data_new[2],50,50), "\n";
            print OUT2 "focalSNP(MSY): ", substr($data_new[0],49,1),substr($data_new[1],49,1),substr($data_new[2],49,1),"\t";
            
            # define parsimony rules:
            # GGG = G
            # TGG = G
            # GTG = G
            # GGT = G
            # TTG = T
            # TGT = T etc.
            # assign ancestry states based on above rules **************************************
            my $anc = ""; my $MSY="";
            # print length($data_new[0]), "\n";
            if (length($data_new[0])>50){  # if all three sequences are present print alignment to log, the first segment ends on the focal SNP
                if (substr($data_new[1],49,1) eq substr($data_new[2],49,1)){print OUT2 "ANC:",substr($data_new[1],49,1), "\n\n"; $anc =substr($data_new[1],49,1)}
                elsif (substr($data_new[0],49,1) eq substr($data_new[2],49,1)){print OUT2 "ANC:",substr($data_new[0],49,1), "\n\n";$anc =substr($data_new[0],49,1)}
                elsif (substr($data_new[0],49,1) eq substr($data_new[1],49,1)){print OUT2 "ANC:",substr($data_new[0],49,1), "\n\n"; $anc =substr($data_new[0],49,1)}
                else {print OUT2 "ANC:MULTI\n\n"; $anc= "AMBIG"}
                
                print OUT3 substr($parsed_snp_ID[0],1),"\t$parsed_snp_ID[1]\t$parsed_snp_ID[2]\t",substr($data_new[0],49,1),substr($data_new[1],49,1),substr($data_new[2],49,1),"\t$anc\n";
                print substr($parsed_snp_ID[0],1),"\t$parsed_snp_ID[1]\t$parsed_snp_ID[2]\t",substr($data_new[0],49,1),substr($data_new[1],49,1),substr($data_new[2],49,1),"\t$anc\n";
                }
            } # if all three sequences are present
            else{ # if no blast hit in either Dsim or Dyak
                print OUT2 "$sequence_names_new[0]\tNo blast hit in either Dsim or Dyak\n",substr($data[0],0,50),"\t",substr($data[0],50,50),"\n\n";  # print Dmel sequence to log
                print OUT3 substr($parsed_snp_ID[0],1),"\t$parsed_snp_ID[1]\t$parsed_snp_ID[2]\tno_hit\tno_hit\n";
                print substr($parsed_snp_ID[0],1),"\t$parsed_snp_ID[1]\t$parsed_snp_ID[2]\tno_hit\t\tno_hit\n";
                }
    
} # for all sequences
system("rm temp.fa temp2.fa");
