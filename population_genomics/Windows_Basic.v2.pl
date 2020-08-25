#!/bin/perl
use strict;
use warnings;


# Takes input files (comma-separated list) of samples in "simplified" vcf form (output of SimplifyVCF_Basic.pl),
# along with a comma-separated list of names for the samples and a few input parameters (minimum depths for
# some analyses, etc.).
#
# Key assumptions about input formats:

# 1. The "simplified" vcf file should have columns of: CHROM POS REF IND1 IND2 (etc) NUM_A NUM_T NUM_G NUM_C.
#    As long as there is information about A,G,T,C calls for each individual at each position, the program can
#    be modified relatively easily to accomodate alternative simplified vcf formats.
# 2. It is assumed that the vcf files contain all the sites in the genome, with any desired masking of sites done
#    at this stage.  Modifications restricting attention to variant sites only should be straightforward.
# 3. The masking method implemented here is done by implementing a gff-format list of masked regions ($RM_gff_file),
#    as outputted by the program RepeatMasker.

# In this implementation, all masked sites/regions are discarded from the analysis
# completely (i.e., we are not simply rejecting variants in masked regions).  Currently, the program will
# also discard any sites for which at least one sample has zero depth of called sites, but other
# screening methods (sites containing indels, etc.) can be easily implemented.
#
# Pi, S, and Tajima's D are computed separately in each sample for each observed sampling depth.  The values
# are then combined using a weighted average for Pi and S, and a modification of the weighting scheme of Langley et al (2012) for
# Tajima's D, for each window of size $window_size. The modification of the Tajima's D calculation is as follows:
# instead of assuming that each depth-specific calculation of Tajima's D is an *independent* realization of the
# usual N(0,1) random variable, we assume that the depth-specific Tajima's D calculations are *perfectly correlated*,
# reflecting the likely genealogical correlation of sites within a window. This results in division of the summed
# Tajima's D estimates by the number of estimates (rather than the sqrt of the number of estimates, as in Langley et al.)
# to produce the approximate N(0,1) final Tajima's D values for a window.
#
# Only sites meeting some criteria $min_TajD_depth will be included in the Tajima's D calculation, but all sites
# are used for the final S and Pi estimates.
#
# For pairwise Fst calculations, the Weir-Cockerham estimator is used, and estimates are combined across SNPs
# by averaging the numerator (across SNPs) and the denominator (across SNPs), and taking the ratio of those
# averages. There is a minimum depth (min_Fst_depth in both samples considered) for inclusion in Fst
# calculations, and only SNPs that are diallelic across the entire dataset are considered.
#
# The Fst minimum sample depth is also used for the pairwise dxy calculations, but here we do not restrict to diallelic sites.
#
# Usage:  Windows_Basic.pl in.vcf_list in.sample_list in.repeat_gff min_TajD_depth min_Fst_depth max_depth window_size outfile




my $sample_vcf_list = $ARGV[0]; #comma-separated list of (simplified) vcfs for each sample
my $sample_name_list = $ARGV[1]; #comma-separated list of names for each sample
my $RM_gff_file = $ARGV[2]; # (bare) gff-formatted list of masked sites; excluded *completely* from analysis
my $min_TajD_depth = $ARGV[3]; # minimum sample depth required to include a variant in Tajima's D calculations (required for focal sample)
my $min_Fst_depth = $ARGV[4]; # minimum sample depth required to include a site in Fst calculations (required for all samples)
my $max_depth = $ARGV[5]; # maximum possible sample depth across all samples
my $window_size = $ARGV[6]; # size (in bp) of non-overlapping windows to consider
my $outfile = $ARGV[7]; # name of output file


print "Starting program now, time is ";
system('echo `date`');


#####################################################################################################################################
#                                                                                                                                   #
#       Assignment of variables for window description and calculations.  Reinitialize each window.                                 #
#                                                                                                                                   #
#####################################################################################################################################

# descriptive and place-counting variables
my $XA = 'NA';
my $chrom_state = 'NA';

my $window_count = 0;
my $window_start = 0;
my $window_end = 0;
my $window_mid = 0;

my $discarded_masked = 0;
my $discarded_ND = 0;
my $discarded_depth = 0;

# basic elements of the pi, S, Fst, and TajD calculations

# Arrays of temporary holders for sample x depth data as we go through the window
# One element per possible sample depth within a sample
my @pi_tmp = ();
my @S_tmp = ();
my @site_count_included = ();
my @site_count_all = ();


# Final calculated quantities, one element per sample; max $num_samples (pulled from input data; last element reserved for the population total)
my @TajD = ();
my @final_pi = ();
my @final_S_count = ();

# vectors to hold running totals for pairwise Fst and dxy quantities (note:  these will be sparsely filled)
my @pw_Fst_num = ();
my @pw_Fst_den = ();

my @pw_dxy_numsites = ();
my @pw_dxy = ();

# scalars to hold running total for the overall Fst
my $Fst_num = 0;
my $Fst_den = 0;

# Final calculated quantities, one element per sample pair
my @final_pw_Fst = ();
my $final_Fst = 'NA';

# vectors to store running totals for depths and overall counts of sites
my @final_site_count_included = ();
my @final_site_count_all = ();
my @mean_depth_included = ();
my @mean_depth_all = ();

# Special cases:  must be re-initialized for each site (not just each window)
my @site_depth_samp = (); # sampling depth by sample
my @num_alleles = (); # number of alleles by sample (simple array)
my @samp_atgc = (); # allele count numbers by sample (including total)


# we'll use this PassFail variable later for tests on whether or not a site should be included. Three tests currently implemented.
my $PassFail = 'PASS';


#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#                                                                                                                                   #
#       Read in and do basic manipulations of input files to make the information available to the analysis                         #
#                                                                                                                                   #
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################


my @sample_vcfs = split(",",$sample_vcf_list);
my @sample_names = split(",",$sample_name_list);
my @sample_sizes = ();
my $obs_max_depth = 0;
my $num_samples = 0;
my %all_samples_hash = (); # giant genome hash.  Keys correspond to sample_Chr:pos pointing to the relevant line of the sample vcf.


# List of chromosomes to analyze, and their lengths. (Could read from the vcfs, but usually we want the major arms only, so do this by hand)
# Later in the script, as part of the window description, euchromatic-heterochromatic regions are defined by coordinates.  Adjust as needed.
#my @chrom_list = ('ref_sim_3L:12790000..13275000', '2L','2R','3L','3R','4','X');
#my %chrom_length = (
#'ref_sim_3L:12790000..13275000' => 484999,
#'2L' => 23513712,
#'2R' => 25286936,
#'3L' => 28110227,
#'3R' => 32079331,
#'4' => 1348131,
#'X' => 23542271,
#);


my @chrom_list = ('1','2');
my %chrom_length = (
'1' => 23513712,
'2' => 25286936,
);

$num_samples = @sample_names;
unless ($num_samples == @sample_vcfs)  {
    die "Number of samples and number of sample vcfs does not match! Program died ";
}
push (@sample_names,"overall");




#####################################################################################################################################
#                                                                                                                                   #
#       Read in all the data from the simplified vcfs.  While we're here, we'll pick up information about                           #
#       the max sampling depth of each sample (i.e., the sample size, recorded in @sample_sizes()).                                 #
#                                                                                                                                   #
#       All raw information from the vcfs is stored in the genome hash %all_samples_hash.                                           #
#       Keys are: sample_chrom:pos, and contents are the corresponding line in the corresponding vcf.                               #
#                                                                                                                                   #
#####################################################################################################################################


for (my $i = 0; $i < $num_samples; $i++){
    
    open(FILE, "$sample_vcfs[$i]")|| print ("Can't open $sample_vcfs[$i]!\n");
    foreach my $line(<FILE>)
    {
        chomp($line);
        if (substr($line,0,5) eq 'CHROM'){
            my @vcf_info = split("\t", $line);
            my $sample_size = 0;
            $sample_size = @vcf_info - 7;
            push @sample_sizes, $sample_size;
            $obs_max_depth += $sample_size;
        }
        else{
            if (substr($line,0,1) ne "#"){
                my @vcf_info = split("\t", $line);
                $all_samples_hash{$sample_names[$i]."_".$vcf_info[0].":".$vcf_info[1]} = $line;
            }
        }
    }
    close(FILE);
    
    unless ($obs_max_depth < $max_depth)  {
        die "Given max depth is too low! Program died ";
    }
    
    print "finished loading ", $sample_names[$i], " sample, sample size is ", $sample_sizes[$i], ", time is ";
    system('echo `date`');
}


#set up hash of masked sites (removed from analysis).  Here, assumes gff format, but easily modified to take in a list of sites
my %repeat_hash = ();
open(FILE, "$RM_gff_file")|| print ("Can't open $RM_gff_file!\n");
foreach my $line(<FILE>)
{
    chomp($line);
    if (substr($line, 0, 1) ne '#'){
        my @repeat_info = split("\t", $line);
        for (my $i = $repeat_info[3]; $i <= $repeat_info[4]; $i++){
            $repeat_hash{$repeat_info[0].":".$i} = $repeat_info[8]; # repeat name; not used currently
        }
    }
}
close(FILE);

print "finished loading masked sites, about to begin main program, time is ";
system('echo `date`');



# now we know what we're dealing with, time to set the column headings for the output
open(OUTFILE, ">$outfile") || (print "Can't open $outfile!\n");
print OUTFILE "chrom,chrom_state,window_start,window_end,window_midpt,window_count,discarded_masked,discarded_ND,discarded_depth,XA";
for (my $i = 0; $i <= $num_samples; $i++) {
    print OUTFILE ",avg_depth_inc_",$sample_names[$i],",avg_depth_all_",$sample_names[$i],",S_",$sample_names[$i],",pi_",$sample_names[$i],",TajD_",$sample_names[$i];
}
for (my $i = 0; $i < $num_samples; $i++){
    for (my $j = $i + 1; $j < $num_samples; $j++){
        print OUTFILE ",Fst_",$sample_names[$i],".",$sample_names[$j];
        print OUTFILE ",dxy_",$sample_names[$i],".",$sample_names[$j];
        print OUTFILE ",dxy_numsites_",$sample_names[$i],".",$sample_names[$j];
    }
}
print OUTFILE ",Fst_overall";
print OUTFILE "\n";
close(OUTFILE);





#initialization-- repeated for every new window
# last element reserved for the population total
for (my $i = 0; $i <= $num_samples; $i++){
    
    for (my $j = 0; $j < $max_depth; $j++){
        $pi_tmp[$i * $max_depth + $j] = 0;
        $S_tmp[$i * $max_depth + $j] = 0;
        $site_count_included[$i * $max_depth + $j] = 0;
        $site_count_all[$i * $max_depth + $j] = 0;
        
    }
    for (my $j = 0; $j < 4; $j++){
        $samp_atgc[$i * 4 + $j] = 0;
    }
    
    $TajD[$i] = 'NA';
    $final_pi[$i] = 0;
    $final_S_count[$i] = 0;
    $final_site_count_included[$i] = 0;
    $final_site_count_all[$i] = 0;
    $mean_depth_included[$i] = 0;
    $mean_depth_all[$i] = 0;
    
    $site_depth_samp[$i] = 0;
    $num_alleles[$i] = 0;
}


for (my $i = 0; $i < $num_samples; $i++){
    for (my $j = $i + 1; $j < $num_samples; $j++){
        $pw_Fst_num[$i * $num_samples + $j] = 0;
        $pw_Fst_den[$i * $num_samples + $j] = 0;
        $final_pw_Fst[$i * $num_samples + $j] = 'NA';
        $pw_dxy_numsites[$i * $num_samples + $j] = 0;
        $pw_dxy[$i * $num_samples + $j] = 0;
        
    }
}


$Fst_num = 0;
$Fst_den = 0;
$final_Fst = 'NA';

# Tajima's D elements-- will calculate the final constants (e1, e2) for all possible sampling depths. Note that
# the vector indices directly correspond to the sample depth (ignoring the unused 0,1 indices).

my @a1 = ();
my @a2 = ();
my @b1 = ();
my @b2 = ();
my @c1 = ();
my @c2 = ();
my @e1 = ();
my @e2 = ();

for (my $i = 2; $i <= $max_depth; $i++){
    
    $a1[$i] = 0;
    $a2[$i] = 0;
    $b1[$i] = 0;
    $b2[$i] = 0;
    $c1[$i] = 0;
    $c2[$i] = 0;
    $e1[$i] = 0;
    $e2[$i] = 0;
    
    for (my $j = 1; $j < $i; $j++){
        $a1[$i] += 1 / $j;
        $a2[$i] += 1 / ($j * $j);
    }
    
    $b1[$i] = ($i + 1) / (3 * ($i - 1));
    $b2[$i] = 2 * ($i * $i + $i + 3) / (9 * $i * ($i - 1));
    
    $c1[$i] = $b1[$i] - (1 / $a1[$i]);
    $c2[$i] = $b2[$i] - (($i + 2) / ($a1[$i] * $i)) + ($a2[$i] / ($a1[$i] * $a1[$i]));
    
    $e1[$i] = $c1[$i] / $a1[$i];
    $e2[$i] = $c2[$i] / ($a1[$i] * $a1[$i] + $a2[$i]);
    
}



for (my $chrom_index = 0; $chrom_index < @chrom_list; $chrom_index++) {
#for (my $chrom_index = 1; $chrom_index < 2; $chrom_index++) {
    
    # set the new chromosome
    my $chrom = $chrom_list[$chrom_index];
    print "Starting chromosome ", $chrom, ", time is ";
    system("echo `date`");
    
    #####################################################################################################################################
    #####################################################################################################################################
    #####################################################################################################################################
    #                                                                                                                                   #
    #   Initialize everything (repeated for each new window and when crossing into a new chromosome-- last, incomplete window of the    #
    #   previous chromosome will be ditched).                                                                                           #
    #                                                                                                                                   #
    #####################################################################################################################################
    #####################################################################################################################################
    #####################################################################################################################################
    $window_start = 1; # only 1 at the beginning of a chromosome, (otherwise current position + 1 for new windows within a chromosome)
    $window_count = 0;
    $window_end = 'NA';
    
    $discarded_ND = 0;
    $discarded_masked = 0;
    $discarded_depth = 0;
    
    $final_Fst = 'NA';
    $Fst_num = 0;
    $Fst_den = 0;
    
    
    for (my $i = 0; $i < $num_samples; $i++){
        for (my $j = $i + 1; $j < $num_samples; $j++){
            $pw_Fst_num[$i * $num_samples + $j] = 0;
            $pw_Fst_den[$i * $num_samples + $j] = 0;
            $final_pw_Fst[$i * $num_samples + $j] = 'NA';
        }
    }
    
    for (my $i = 0; $i <= $num_samples; $i++){
        
        for (my $j = 0; $j < $max_depth; $j++){
            $pi_tmp[$i * $max_depth + $j] = 0;
            $S_tmp[$i * $max_depth + $j] = 0;
            $site_count_included[$i * $max_depth + $j] = 0;
            $site_count_all[$i * $max_depth + $j] = 0;
            
        }
        for (my $j = 0; $j < 4; $j++){
            $samp_atgc[$i * 4 + $j] = 0;
        }
        
        $TajD[$i] = 'NA';
        $final_pi[$i] = 0;
        $final_S_count[$i] = 0;
        $final_site_count_included[$i] = 0;
        $final_site_count_all[$i] = 0;
        $mean_depth_included[$i] = 0;
        $mean_depth_all[$i] = 0;
        
        $site_depth_samp[$i] = 0;
        $num_alleles[$i] = 0;
    }
    
    for (my $pos = 1; $pos <= $chrom_length{$chrom}; $pos++){
        $PassFail = 'PASS';
        #make it slightly easier to look up the genome results
        my $posName = $chrom.":".$pos;
        
        #####################################################################################################################################
        #                                                                                                                                   #
        #   First two screening tests (third is internal, after assessment of the samples). Note that if a site fails one of these first    #
        #   tests, it will NOT be included in the "all" counts of sites (when looking at the overall sampling depth in a region for both    #
        #   included and excluded sites.                                                                                                    #
        #                                                                                                                                   #
        #####################################################################################################################################
        
        # First screening test-- do we even have data for this site in all the samples?
        for (my $sample = 0; $sample < $num_samples; $sample++){
            if (!(exists $all_samples_hash{$sample_names[$sample]."_".$chrom.":".$pos})){
                $PassFail = 'FAIL';
                $discarded_ND++;
            }
        }
        
        # Second screening test-- is this site in a masked region?  If so, don't even start the analysis
        if (exists $repeat_hash{$posName}){
            $PassFail = 'FAIL';
            $discarded_masked++;
        }
        
        if ($PassFail eq 'PASS'){
            
            #initialize the site-specific vectors for number of alleles, sample depth, and atgc counts
            for (my $sample = 0; $sample <= $num_samples; $sample++){
                $site_depth_samp[$sample] = 0;
                $num_alleles[$sample] = 0;
                for (my $allele = 0; $allele < 4; $allele++){
                    $samp_atgc[$sample * 4 + $allele] = 0;
                }
            }
 
            #####################################################################################################################################
            #####################################################################################################################################
            #####################################################################################################################################
            #                                                                                                                                   #
            #   Load data for samples in preparation for calculations.  Note that the total population slot ($num_samples) has to be filled     #
            #   in as a running total for some things (e.g. depth), or else assessed fully after all samples have loaded (e.g. num_alleles).    #
            #                                                                                                                                   #
            #####################################################################################################################################
            #####################################################################################################################################
            #####################################################################################################################################

             for (my $sample = 0; $sample < $num_samples; $sample++){
                
                #########################################################################################################################
                #                                                                                                                       #
                #   For the "simplfied" vcf format produced by SimplifyVCF_Basic.pl, for each chrom:pos location, the counts            #
                #   of A,T, G, C alleles are going to be in locations (3+samp_size) + 0,1,2,3.  Adjust as needed for other              #
                #   vcf formats (e.g., ATGC totals can be easily calculated on-the-fly if not in the VCF).                              #
                #                                                                                                                       #
                #########################################################################################################################
                
                my $cur_samp = $all_samples_hash{$sample_names[$sample]."_".$chrom.":".$pos};
                my @samp_info = split("\t", $cur_samp);
                
                #get the observed sampling depth and number of alleles present in the sample
                for (my $j = $sample_sizes[$sample] + 3; $j < $sample_sizes[$sample] + 7; $j++){
                    if ($samp_info[$j] > 0){
                        $num_alleles[$sample]++;
                        $site_depth_samp[$sample] += $samp_info[$j];
                        
                        my $local_allele = $j - $sample_sizes[$sample] - 3;
                        #record the agct info for this sample in compact, easily retrievable form
                        $samp_atgc[$sample * 4 + $local_allele] = $samp_info[$j];
                        
                        #total pop (in the $num_samples slot) is a running total, so add these in there
                        $site_depth_samp[$num_samples] += $samp_info[$j];
                        $samp_atgc[$num_samples * 4 + $local_allele] += $samp_info[$j];
                    }
                }
                
                
                # Add this site's info to the site counts ("all") for avg depth calculations
                $site_count_all[$sample * $max_depth + $site_depth_samp[$sample]]++;
                $mean_depth_all[$sample] += $site_depth_samp[$sample];
                
                
                #####################################################################################################################################
                #                                                                                                                                   #
                #   Third screening test.  This one is an example of a mid-stream test of data quality.                                             #
                #   Here, we're just making sure we have some data for each sample. If any fail, the flag is fipped                                 #
                #   and the whole site fails and is excluded from analysis (site is retained in "all" counts, however)                              #
                #                                                                                                                                   #
                #####################################################################################################################################
                 
                if ($site_depth_samp[$sample] == 0) {
                    $PassFail = 'FAIL'
                }
                
                
            } #end of site data loading for samples
            
            # fill in the total number of alleles
            for (my $j = 0; $j < 4; $j++){
                if ($samp_atgc[$num_samples * 4 + $j] > 0){
                    $num_alleles[$num_samples]++;
                }
            }
            
            # Total population ($num_samples) running totals
            $site_count_all[$num_samples * $max_depth + $site_depth_samp[$num_samples]] = $site_count_all[$num_samples * $max_depth + $site_depth_samp[$num_samples]] + 1;
            $mean_depth_all[$num_samples] += $site_depth_samp[$num_samples];
            
            
            if ($PassFail eq 'PASS'){
                
                $window_count++; #add a site to the current count
                
                
                #####################################################################################################################################
                #####################################################################################################################################
                #####################################################################################################################################
                #                                                                                                                                   #
                #   Pi and S running totals.  Note that these running counts are always tied to their sampling depths.                              #
                #                                                                                                                                   #
                #####################################################################################################################################
                #####################################################################################################################################
                #####################################################################################################################################
                
                
                for (my $cur_sample = 0; $cur_sample <= $num_samples; $cur_sample++) {
                    my $cur_sample_depth = $site_depth_samp[$cur_sample];
                    
                    $mean_depth_included[$cur_sample] += $cur_sample_depth;
                    $site_count_included[$cur_sample * $max_depth + $cur_sample_depth]++;
                    $final_site_count_included[$cur_sample]++;
                    
                    # S and pi for this sample-depth combination added to the totals
                    if ($num_alleles[$cur_sample] > 1) {
                        $S_tmp[$cur_sample * $max_depth + $cur_sample_depth]++;
                        my $tmp_pi = 1;
                        for (my $allele = 0; $allele < 4; $allele++) {
                            if ($samp_atgc[$cur_sample * 4 + $allele] > 0) {
                                $tmp_pi *= $samp_atgc[$cur_sample * 4 + $allele];
                            }
                        }
                        $pi_tmp[$cur_sample * $max_depth + $cur_sample_depth] += $tmp_pi;
                        
                    }
                }
                
                
                #####################################################################################################################################
                #                                                                                                                                   #
                #                                               Pairwise Dxy calculations                                                           #
                #                                                                                                                                   #
                #####################################################################################################################################
                
                for (my $first_sample = 0; $first_sample < $num_samples; $first_sample++){
                    for (my $second_sample = $first_sample + 1; $second_sample < $num_samples; $second_sample++){
                        
                        my $n1 = $site_depth_samp[$first_sample]; # size of first sample at this site
                        my $n2 = $site_depth_samp[$second_sample]; # size of second sample
                        
                        # requirements for a "good" site for Dxy calculations:
                        # both samples meeting the minimum depth requirement
                        if (($n1 >= $min_Fst_depth) && ($n2 >= $min_Fst_depth)){
                            
                            # we've found a usable site, now go through to find the polymorphic alleles and do the calculations
   
                            
                            my $num_comparisons = $n1 * $n2;
                            my $num_diffs = 0;
                            
                            for (my $allele = 0; $allele < 4; $allele++){
                                if ($samp_atgc[$first_sample * 4 + $allele] != 0){
                                    $num_diffs += ($samp_atgc[$first_sample * 4 + $allele]) * ($n2 - $samp_atgc[$second_sample * 4 + $allele]);
                                }
                            }
                            
                            $pw_dxy_numsites[$first_sample * $num_samples + $second_sample] += 1;
                            $pw_dxy[$first_sample * $num_samples + $second_sample] += $num_diffs / $num_comparisons;
                            
                            
                        }
                    }
                }
                

                
                
                
                #####################################################################################################################################
                #####################################################################################################################################
                #####################################################################################################################################
                #                                                                                                                                   #
                #                                               Fst calculations (Weir-Cockerham estimator)                                         #
                #                                                                                                                                   #
                #       All pairwise combinations plus an overall Fst.  Observed heterozygosity of sample is assumed to be zero.                    #
                #                               No attempt is made to weight by sampling depth of an individual SNP.                                #
                #                                                                                                                                   #
                #                                           ONLY DIALLELIC SNPS ARE CONSIDERED                                                      #
                #                                                                                                                                   #
                #####################################################################################################################################
                #####################################################################################################################################
                #####################################################################################################################################
                
                
                #first do the overall Fst estimate
                if ($num_alleles[$num_samples] == 2){
                    my $r = 0; # number of samples included in estimate
                    my $n_sum = 0; # total sample size (depth)
                    my $n2_sum = 0;
                    
                    my $p_bar = 0; # average frequency of focal allele
                    my $n_bar = 0; # average sample size (depth)
                    my @p_vec = (); # vector of focal allele frequencies by sample
                    my @n_vec = (); # vector of sample sizes (depths)
                    
                    #pick the focal allele
                    my $allele = 0;
                    while ($samp_atgc[$num_samples * 4 + $allele] == 0){
                        $allele++;
                    }
                    
                    #only include those samples for which depth is high enough to justify inclusion, add their info to the temporary nvec and pvec data
                    for (my $i = 0; $i < $num_samples; $i++){
                        if ($site_depth_samp[$i] >= $min_Fst_depth){
                            $r++;
                            push (@n_vec, $site_depth_samp[$i]);
                            push (@p_vec, $samp_atgc[$i * 4 + $allele] / $site_depth_samp[$i]);
                            $p_bar += $samp_atgc[$i * 4 + $allele];
                            $n_sum += $site_depth_samp[$i];
                            $n2_sum += $site_depth_samp[$i] * $site_depth_samp[$i];
                        }
                    }
                    
                    #don't bother with further calculations unless you have at least two samples included
                    if ($r > 1){
                        
                        my $nc = 0;
                        my $S2 = 0;
                        my $T1 = 0;
                        my $T2 = 0;
                        
                        $p_bar /= $n_sum;
                        $n_bar = $n_sum / $r;
                        $nc = (1.0 / ($r - 1.0))*($n_sum - ($n2_sum / $n_sum));
                        
                        for (my $i = 0; $i < $r; $i++){
                            $S2 += $n_vec[$i] * ($p_vec[$i] - $p_bar) * ($p_vec[$i] - $p_bar);
                        }
                        
                        $S2 *= (1.0/($n_bar * ($r - 1.0)));
                        
                        $T1 = $S2 - (1.0 / (2.0 * $n_bar - 1.0)) * ($p_bar * (1.0 - $p_bar) - ($r - 1.0) * $S2/ $r);
                        $T2 = $p_bar * (1.0 - $p_bar) * (2.0 * $nc - 1.0) / (2.0 * $n_bar - 1.0) + ($S2 / $r) * (1.0 + (2.0 * ($r - 1.0) * ($n_bar - $nc))/(2.0 * $n_bar - 1.0));
                        
                        $Fst_num += $T1;
                        $Fst_den += $T2;
   
                        undef @p_vec;
                        undef @n_vec;
                    }
                }
                
                
                #####################################################################################################################################
                #                                                                                                                                   #
                #                                               Pairwise Fst calculations                                                           #
                #                                                                                                                                   #
                #####################################################################################################################################
               
                for (my $first_sample = 0; $first_sample < $num_samples; $first_sample++){
                    for (my $second_sample = $first_sample + 1; $second_sample < $num_samples; $second_sample++){
                        
                        my $n1 = $site_depth_samp[$first_sample]; # size of first sample at this site
                        my $n2 = $site_depth_samp[$second_sample]; # size of second sample
                        
                        # requirements for a "good" site for Fst calculations:  diallelic overall,
                        # with both samples meeting the minimum depth requirement
                        if (($num_alleles[$num_samples] == 2) && ($n1 >= $min_Fst_depth) && ($n2 >= $min_Fst_depth)){
                            
                            # we've found a usable site, now go through to find the polymorphic allele and do the calculations
                            
                            #pick the focal allele
                            my $allele = 0;
                            while ($samp_atgc[$first_sample * 4 + $allele] == 0){
                                $allele++;
                            }
                            
                            
                            my $n_sum = $n1 + $n2;
                            my $n2_sum = $n1 * $n1 + $n2 * $n2;
                            
                            my $p1 = $samp_atgc[$first_sample * 4 + $allele];
                            my $p2 = $samp_atgc[$second_sample * 4 + $allele];
                            
                            my $p_bar = ($p1 + $p2) / $n_sum;
                            
                            $p1 /= $n1;
                            $p2 /= $n2;
                            
                            my $n_bar = $n_sum / 2.0;
                            my $nc = ($n_sum - ($n2_sum / $n_sum));
                            
                            
                            my $S2 = ($n1 * ($p1 - $p_bar) * ($p1 - $p_bar) + $n2 * ($p2 - $p_bar) * ($p2 - $p_bar)) / $n_bar;
                            
                            
                            my $T1 = $S2 - (1.0 / (2.0 * $n_bar - 1.0)) * ($p_bar * (1.0 - $p_bar) - $S2 / 2);
                            my $T2 = $p_bar * (1.0 - $p_bar) * (2.0 * $nc - 1.0) / (2.0 * $n_bar - 1.0) + ($S2 / 2) * (1.0 + (2.0 * ($n_bar - $nc))/(2.0 * $n_bar - 1.0));
                            
                            
                            if ($T2 > 0){
                                $pw_Fst_num[$first_sample * $num_samples + $second_sample] += $T1;
                                $pw_Fst_den[$first_sample * $num_samples + $second_sample] += $T2;
                            }
                            
                            
                        }
                    }
                }
                
   
                
                
                
                
                
                #####################################################################################################################################
                #####################################################################################################################################
                #####################################################################################################################################
                #                                                                                                                                   #
                #   Check if we've completed a window.  If so, stop, assess where we are, do final pi, Tajima's D, and Fst calculations,            #
                #   and print the results.                                                                                                          #
                #                                                                                                                                   #
                #####################################################################################################################################
                #####################################################################################################################################
                #####################################################################################################################################
                
                
                if ($window_count == $window_size) {
                    
                    # collect remaining descriptive information about the window
                    $window_end = $pos;
                    $window_mid = 0.5 * ($window_end + $window_start);
                    
                    $XA = 'A'; # default to simple X vs autosome designation, add other categories if needed
                    if ($chrom eq 'X') {
                        $XA = 'X';
                    }
                    
                    if ((($chrom eq '2L') && ($window_mid>=1000000) && ($window_mid<=20000000)) ||
                        (($chrom eq '2R') && ($window_mid>=6500000) && ($window_mid<=20000000)) ||
                        (($chrom eq '3L') && ($window_mid>=1000000) && ($window_mid<=18000000)) ||
                        (($chrom eq '3R') && ($window_mid>=8500000) && ($window_mid<=26500000)) ||
                        (($chrom eq 'X') && ($window_mid>=2500000) && ($window_mid<=19500000))) {
                            $chrom_state = 'EU';
                        }
                    else {
                        $chrom_state = "HET";
                    }
                    
                    
                    #####################################################################################################################################
                    #                                                                                                                                   #
                    #   Tajima's D calculations.  Uses a modification of the Langley et al (2012) method of averaging across different sampling depths. #
                    #   D values are summed for all sites with at least $min_TajD_depth across all depths with at least 4 observed sites, divided by    #
                    #   the number of depths retained.                                                                                                  #
                    #                                                                                                                                   #
                    #####################################################################################################################################
                    
                    
                    for (my $cur_sample = 0; $cur_sample <= $num_samples; $cur_sample++){
                        
                        
                        #Tajima's D. Only use the subset of data for which sampling depth exceeds some cutoff ($min_sample_depth).
                        #temporarily set TajD to zero for calculation purposes
                        $TajD[$cur_sample] = 0;
                        
                        my $num_depths_used = 0;
                        # Calculate a TajD for each depth (greater than min_depth), as long as there are at least 4 observations at that depth
                        for (my $cur_depth = $min_TajD_depth; $cur_depth < $max_depth; $cur_depth++) {
                            if ($S_tmp[$cur_sample * $max_depth + $cur_depth] > 3){
                                
                                $num_depths_used++;
                                #retrieve the "local" pi and S values for this sample and this depth
                                my $local_pi = $pi_tmp[$cur_sample * $max_depth + $cur_depth] * 2 / ($cur_depth * ($cur_depth - 1));
                                my $local_S = $S_tmp[$cur_sample * $max_depth + $cur_depth];
                            
                                
                                # we need this count to make the final adjustment to the sum of Tajima's D's
                                $final_S_count[$cur_sample] += $local_S;
                                
                                 #TajD calculation for this depth
                                $TajD[$cur_sample] += (($local_pi - ($local_S / $a1[$cur_depth])) / sqrt($e1[$cur_depth] * $local_S + $e2[$cur_depth] * $local_S * ($local_S - 1)));
                            }
                        }
                        
                        if ($num_depths_used > 0){
                            $TajD[$cur_sample] = $TajD[$cur_sample] / $num_depths_used;
                        }
                        else{
                            $TajD[$cur_sample] = 'NA';
                        }
                        
                        #####################################################################################################################################
                        #                                                                                                                                   #
                        #          Fully recalculate pi, S values to include low-sample-size elements skipped in the TajD calculations                      #
                        #                                                                                                                                   #
                        #####################################################################################################################################
                        
                        $final_pi[$cur_sample] = 0;
                        $final_S_count[$cur_sample] = 0;
                        
                        $final_site_count_included[$cur_sample] = $site_count_included[$cur_sample * $max_depth + 1]; # depth 1 sites (can't go through loop below)
                        $final_site_count_all[$cur_sample] = $site_count_all[$cur_sample * $max_depth + 1]; # depth 1 sites
                        
                        
                        for (my $cur_depth = 2; $cur_depth < $max_depth; $cur_depth++) {
                            
                            $final_pi[$cur_sample] += $pi_tmp[$cur_sample * $max_depth + $cur_depth] * 2 / ($cur_depth * ($cur_depth - 1));
                            $final_S_count[$cur_sample] += $S_tmp[$cur_sample * $max_depth + $cur_depth];
                            $final_site_count_included[$cur_sample] += $site_count_included[$cur_sample * $max_depth + $cur_depth];
                            $final_site_count_all[$cur_sample] += $site_count_all[$cur_sample * $max_depth + $cur_depth];
                            
                        }
    
                        #####################################################################################################################################
                        #                                                                                                                                   #
                        #                               Final mean depths, and check for "no data" windows                                                  #
                        #                                                                                                                                   #
                        #####################################################################################################################################
                      
                        $mean_depth_all[$cur_sample] /= $final_site_count_all[$cur_sample];

                        #these seem like silly checks, but may be necessary at small window sizes depending on how sites are designated for inclusion/exclusion
                        if ($final_site_count_included[$cur_sample] > 0){
                            $mean_depth_included[$cur_sample] /= $final_site_count_included[$cur_sample];
                            $final_pi[$cur_sample] /= $final_site_count_included[$cur_sample];
                        }
                        else{
                            $mean_depth_included[$cur_sample] = 'NA';
                            $final_pi[$cur_sample] = 'NA';
                            $final_S_count[$cur_sample] = 'NA';
                        }
                        
                        
                        
                    }
                    
                    #####################################################################################################################################
                    #                                                                                                                                   #
                    #                                                   Final Fst calculations                                                          #
                    #                                                                                                                                   #
                    #####################################################################################################################################
                    
                    if ($Fst_den > 0) {
                        $final_Fst = $Fst_num / $Fst_den;
                     }
                    else{
                        $final_Fst = 'NA';
                    }
                    
                    for (my $i = 0; $i < $num_samples; $i++){
                        for (my $j = $i + 1; $j < $num_samples; $j++){
                            
                            if ($pw_Fst_den[$i * $num_samples + $j] > 0){
                                $final_pw_Fst[$i * $num_samples + $j] = ($pw_Fst_num[$i * $num_samples + $j] / $pw_Fst_den[$i * $num_samples + $j]);
                            }
                            else {
                                $final_pw_Fst[$i * $num_samples + $j] = 'NA';
                            }
                        }
                    }
                    
                    
                    
                    
                    #####################################################################################################################################
                    #####################################################################################################################################
                    #####################################################################################################################################
                    #                                                                                                                                   #
                    #                                               Print out the results                                                               #
                    #                                                                                                                                   #
                    #####################################################################################################################################
                    #####################################################################################################################################
                    #####################################################################################################################################
                    
                    open(OUTFILE, ">>$outfile")|| print ("Can't open $outfile!\n");
                    print OUTFILE $chrom, ",",$chrom_state,",",$window_start,",",$window_end,",",$window_mid,",",$window_count,",";
                    print OUTFILE $discarded_masked,",",$discarded_ND,",",$discarded_depth,",",$XA;
                    for (my $i = 0; $i <= $num_samples; $i++) {
                        print OUTFILE ",",$mean_depth_included[$i],",",$mean_depth_all[$i],",",$final_S_count[$i],",",$final_pi[$i],",",$TajD[$i];
                    }
                    for (my $i = 0; $i < $num_samples; $i++){
                        for (my $j = $i + 1; $j < $num_samples; $j++){
                            print OUTFILE ",",$final_pw_Fst[$i * $num_samples + $j];
                            print OUTFILE ",",$pw_dxy[$i * $num_samples + $j];
                            print OUTFILE ",",$pw_dxy_numsites[$i * $num_samples + $j];
                        }
                        
                    }
                    print OUTFILE ",",$final_Fst;
                    print OUTFILE "\n";
                    close(OUTFILE);
                    
                    
                    
                    #####################################################################################################################################
                    #####################################################################################################################################
                    #####################################################################################################################################
                    #                                                                                                                                   #
                    #                                           Initialize everything for the next window                                               #
                    #                                                                                                                                   #
                    #####################################################################################################################################
                    #####################################################################################################################################
                    #####################################################################################################################################
                    $window_start = $pos+1;
                    $window_count = 0;
                    $window_end = 'NA';
                    
                    $discarded_ND = 0;
                    $discarded_masked = 0;
                    $discarded_depth = 0;
                    
                    $final_Fst = 'NA';
                    $Fst_num = 0;
                    $Fst_den = 0;
                    
                    
                    for (my $i = 0; $i < $num_samples; $i++){
                        for (my $j = $i + 1; $j < $num_samples; $j++){
                            $pw_Fst_num[$i * $num_samples + $j] = 0;
                            $pw_Fst_den[$i * $num_samples + $j] = 0;
                            $final_pw_Fst[$i * $num_samples + $j] = 'NA';
                            $pw_dxy_numsites[$i * $num_samples + $j] = 0;
                            $pw_dxy[$i * $num_samples + $j] = 0;

                        }
                    }
                    
                    for (my $i = 0; $i <= $num_samples; $i++){
                        for (my $j = 0; $j < 4; $j++){
                            $samp_atgc[$i * 4 + $j] = 0;
                        }
                        
                        for (my $j = 0; $j < $max_depth; $j++){
                            $pi_tmp[$i * $max_depth + $j] = 0;
                            $S_tmp[$i * $max_depth + $j] = 0;
                            $site_count_included[$i * $max_depth + $j] = 0;
                            $site_count_all[$i * $max_depth + $j] = 0;
                            
                        }
                        $TajD[$i] = 'NA';
                        $final_pi[$i] = 0;
                        $final_S_count[$i] = 0;
                        $final_site_count_included[$i] = 0;
                        $final_site_count_all[$i] = 0;
                        $mean_depth_included[$i] = 0;
                        $mean_depth_all[$i] = 0;
                        
                        $site_depth_samp[$i] = 0;
                        $num_alleles[$i] = 0;
                    }
                    
                    
                } #end recording window
            } # end internal "passing site" routine
            else {
                $discarded_depth++;
            }
        } #end external "passing site" routine
    } #end loop over positions
} #end loop over chromosomes


print "All done!\n";
system("echo Time is `date`");
