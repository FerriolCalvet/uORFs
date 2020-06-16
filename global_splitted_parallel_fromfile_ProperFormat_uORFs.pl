#!/usr/bin/perl

use vars qw($USAGE);

use strict;
use warnings;
use Bio::EnsEMBL::Registry; 
use Bio::SeqIO;
use Getopt::Long;
use Parallel::ForkManager;# https://metacpan.org/pod/Parallel::ForkManager

my $time_start0 = time;


$USAGE = "perl file.pl [--help] [--output_file] [--output_directory] [--max_processes] [--verbose] [--extension label] [--species] --input_file ensemblIDs.txt\n";

my ($input_filename, $output_filename, $output_directory, $max_processes, $label_file, $species, $verb, $help) =
   (undef,           undef,            undef,             undef,          undef,       "Human",  undef, undef);

&GetOptions('input_file|i=s'          => \$input_filename,
            'output_file|of=s'        => \$output_filename,
            'output_directory|od=s'   => \$output_directory,
            'max_processes|m=s'       => \$max_processes,
            'extension|e=s'           => \$label_file,
            'species|sp=s'            => \$species,
            'verbose|v'               => \$verb,
            'help|h'                  => \$help,
            );

# We could use sort -R to randomly sort the rows of a file so that the IDs of protein
#  coding genes are randomly assigned among processes in a more proportional way
#### system "sort -R filename > new_file";

if ($help) {
   print $USAGE;
   die;
}

if (!defined $input_filename) {
   die($USAGE . "\nPlease specify an input filename.\n");
}

my $default_output = "/home/fcalvet/Desktop/uORFs/DATA/raw_data/raw_uORFs/uORFs_defaultOUT/";

if (!defined $output_filename) {
   my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
   $mon += 1;
   $year += 1900;
   if (!defined $label_file) {
      $output_filename = "${default_output}output_uORFs_$mday-$mon-${year}_at_$hour-$min-$sec.tsv";
   }
   else {
      system "mkdir -p ${default_output}${label_file}";
      $output_filename = "${default_output}${label_file}/${label_file}.output_uORFs_$mday-$mon-${year}_at_$hour-$min-$sec.tsv";
   }
}

### START define stats of uORFs per gene output

my @spl = split('/', reverse($output_filename), 2);
my $stats_file = reverse($spl[0]);
my $stats_path = reverse($spl[1]);
system "mkdir -p $stats_path/stats\n";
my $output_stats = "$stats_path/stats/$stats_file.stats";

### END define stats of uORFs per gene output




### START defining start and stop codons

my $selected_start_codons="atg";
#my $selected_start_codons="atg|ttg|gtg|ctg|aag|agg|acg|ata|att|atc";

my $selected_stop_codons="taa|tga|tag";
# my $selected_stop_codons="taa|tga|tag";

### END defining start and stop codons




# my $input_filename = '/home/fcalvet/Desktop/uORFs/input_data/ppp.txt';
#                                                              more_ppp.txt';
#                                                              allEnsemblGeneIDs_withoutPAR.txt';
#                                                              allEnsemblGeneIDs_noVersion.txt';




### START of the function

# This function receives a transcript object as input, and return all the uORFs that can be found
# in it as output.
# It is based on regular expressions to find start and end codons and then find the matches by 
# finding the first match in frame
sub uORFs_explorer {

   ### START initialize several variables

   my $s;            # start of the transcript
   my $e;            # end of the transcript
   my $updated_e;    # use this end coordinate to update the end points of the transcript
   my $orf="";       # open reading frames
   my $prot_seq="";  # open reading frames translated_sequence
   my $uORFtype;     # type of ORF
   my $length;       # ORF nucleotide length
   my $aalength;     # ORF aminoacid length
   my $uORF_in_trans; # ID to identify the uORFs individually
   $uORF_in_trans = 0;# we initialize it at 1

   my $output="";    # this is the text that will be returned as output of the function

   ### END initialize several variables




   ###******* TO AVOID UNNECESSARY COPIES, I WILL NOT CREATE THIS VARIABLE
   ###******* I have made this change and it keeps working

#   ### START receive the input of the function

#   my $transcript_object=$_[0];

#   ### END receive the input of the function




   ### START distribute the input into different variables the will be used afterwards

   # transcript genomic sequence without introns
   my $dna=Bio::Seq->new(-seq => $_[0]->spliced_seq);
   my $sequence = $dna->seq; # store the sequence here

   # strand that is transcribed
   my $strand_input=$_[0]->strand();

   # get the transcript id
   my $transcript_id=$_[0]->stable_id();

   # get the gene id
   my $gene_id=$_[0]->get_Gene()->stable_id();

   # get the name of the chromosome
   my $name_tran = $_[0]->seq_region_name();

   # CDS start and end coordinates relative to the cDNA sequence (after splicing)
   my $cds_start_int=$_[0]->cdna_coding_start();
   my $cds_end_int=$_[0]->cdna_coding_end();


   # build transcript mapper with the information
   # of the current transcript
   # this will allow us to change from cDNA to genomic coordinates
   my $trmapper = Bio::EnsEMBL::TranscriptMapper->new($_[0]);

   # we store the start and end coordinates of the transcript. start < end.
   # this will be useful for removing the partial CDS transcripts and ORFs
   # identical to those from the transcript at which we are looking
   my $transcript_start = $_[0]->start();
   my $transcript_end = $_[0]->end();


   ### END distribute the input into different variables the will be used afterwards




   ### START adjusting data to the corresponding strand

   # we define the variable strand to use the + and - signs.

   # adjusting the direction of the sequence is not needed because it
   # it is already provided in the proper direction of transcription

   # transcripts in the reverse strand use the
   # reverse complement of the forward sequence
   # but when retrieving the sequence from the
   # transcript it is already taken into account


   # adjust the value of strand accordingly
   my $strand;
   if ($strand_input == 1){
      $strand = '+';
   } else {
      $strand = '-';
   };

   ### END adjusting data to the corresponding strand




   ### START print all the information to check that it is OK

   # print "Sequence:\n", $sequence,"\nStrand:$strand_input\nTranscriptID:$transcript_id\nGeneID:$gene_id\nChromosomeName:$name_tran\nCDSstartPosition:$cds_start_int\n";

   ### END print all the information to check that it is OK




   ### START look for start and end ORF candidates

   # considering only ATG codon for start
   # considering only TAA|TGA|TAG STOP codons

   # store all the start candidates
   my @starts=();

   # store all the end candidates
   my @ends=();


   # !!! ORFs using these starts will be tagged as '.hypothetical_start.'
   # Introduce the first three positions as start candidates
   # accounting for possible ORFs initiated upstream of the TSS
   # generating ORFs not starting with an ATG (or M)
   for (my $frame=0;$frame<3;$frame++) {
      # prevent the inclusion of a STOP codon as a ORF "start" candidate
      # unless ($sequence=~m/^.{$frame}(taa|tga|tag)/i) {
      # also prevent the inclusion of a start codon, as it would be repeated twice
      unless ($sequence=~m/^.{$frame}($selected_stop_codons|atg)/i) {
         push @starts,$frame+1;
#         print "hypo start candidate: ", $frame+1, "\n";
      }
   }


   # Look for ATG matches in any frame and store them into the list
#   while ($sequence=~m/(atg)/gi) {
   while ($sequence=~m/($selected_start_codons)/gi) {
      # only add the ATG if it is downstream the TIS
      # check whether the ATG is downstream compared to
      # the ATG initiating the main ORF
      # we must make sure that the initial position of the codon is 
      # smaller than the initial position of the annotated
      # CDS of that transcript

      # pos(sequence) is the last position of the match
      # so to get the position of the first letter, we substract 2
      if ((pos($sequence)-2) < $cds_start_int){
         push @starts,(pos($sequence)-2);
#         if ($transcript_id eq "ENST00000612152"){
#            print "start candidate: ", pos($sequence)-2, "\n";
#         }
      }
      # if we knew that the matches are obtained in order, we could stop this loop earlier, but it seems they are not
   }

   # Look for STOP codon matches in any frame and store them into the list
   while ($sequence=~m/($selected_stop_codons)/gi) {
      # again
      # pos(sequence) is the last position of the match
      # so to get the position of the first letter, we substract 2
      push @ends,pos($sequence)-2;
   }

   # !!! ORFs using these ends will be tagged as '.hypothetical_end.'
   # introduce the last three positions of the sequence as candidate ends
   # this is because in some cases, a STOP codon is not present
   # within the transcript so the ORF needs to be recorded, and this is the way to do it
   push @ends,( (length $sequence)-2, (length $sequence)-1, (length $sequence));

   ### END look for start and end ORF candidates




#   ### START sort starts and ends increasingly

#   The way we add the elements into the list should ensure that
#   the order is the correct one, this is not needed
#   if you think that these can solve your problems,
#   uncomment them.

#?   Why is the order important?
#?   It is important because the ORFs end at the first inframe STOP
#?   codon after the ATG, so the ends must be iterated in increasing order

#   # not sure this is necessary but it is safer having it
#   my @sorted_starts= sort { $a <=> $b }  @starts;
#   my @sorted_ends= sort { $a <=> $b } @ends;

#   ### END sort starts and ends increasingly




   ### START matching start and ends

   for my $s (@starts) {
      # we iterate only the meaningful ends
      # we use grep to iterate only the ends that are greater
      # than the starting position
      # for my $e (@ends) { 
      #    if ( $e%3==$s%3 and $e>$s) {

      for my $e (grep {$_ > $s} @ends) { 
         if ( $e%3==$s%3 ) {
      #       print "$s start and $e end selected\n";

   ### END matching start and ends

            $uORF_in_trans = $uORF_in_trans + 1;


            ### START compute the lenght of the uORF

            # compute the length by making a substruction
            # the position of the first nucleotide of the STOP minus
            # the position of the first nucleotide of the ATG
            $length=$e-$s;
            # we compute the aminoacid length by dividing by 3
            # (the STOP codon is not included)
            $aalength = $length/3;

##            $e=$e; # end coordinate not included
#            $e=$e-1; # end coordinate included
#            if ($e < $dna->length-2){
#               $orf=Bio::Seq->new(-seq=>$sequence)->subseq($s,$e+2); # include STOP codon
#            } else {
#               $orf=Bio::Seq->new(-seq=>$sequence)->subseq($s,$e-1); # without STOP codon
#            }

            ### END compute the lenght of the uORF




            ### START update the end of the transcript

            # correct the end of the transcript to include the STOP codon
            # Until here, $e stores the first coordinate of the STOP codon
            # in all the cases, except the two last positions of the transcript,
            # we can include the next two positions as they will be the
            # ones completing the STOP codon

            # The two last positions, imply that there is a hypothetical end so that
            # the protein does not finish with an stop codon.


            # !!!
            # when $e == $dna->length-2, we also add a new codon, but in
            # this case it is not necessarily a STOP codon.
            # This is the case of all the CDSs tagged as EndNotFound.
            my $updated_e;
            if ($e+2 <= (length $sequence) ){
#               print "before $e \n";
               # include STOP codon
               $updated_e=$e+2;
#               print "after $updated_e \n";

               #!!!!!!! add +3 to the nucleotide sequence length
               #!!!!!!! add +1 to the protein sequence length

            } else {
               # until the end of the transcript following the preferred frame
               # (in some cases, a codon may not appear complete, !!! not sure about this)
               $updated_e=$e; 
            }

            ### END update the end of the transcript




            ### START get DNA and protein sequence

            # retrieve the sequence of the ORF using the start and end coordinates defined
            # (remember this is nucleotide sequence)
            # don't know if it is worth to declare it as a new variable
            $orf=Bio::Seq->new(-seq=>$sequence)->subseq($s,$updated_e);
            $prot_seq=Bio::Seq->new(-seq=>$orf)->translate->seq;

            ### END get DNA and protein sequence




            ### START remove uORFs that are exactly the same as partial CDSs or start inside the CDS of the annotated transcript

            # if we want to use this, we have to move it below the update of s_new and e_new
#            if ($e_new == $transcript_end && $s_new == $transcript_start){
#               last;
#            } els

            if ($s >= $cds_start_int || $updated_e == $cds_end_int){ # I dont know if this is removing transcripts that end
                                                                     # the same but are not the same, is this possible, or
                                                                     # if they end in the same way they should be the same??
               last;
            }

            ### END remove uORFs that are exactly the same as partial CDSs or start inside the CDS of the annotated transcript




            ### START label the transcript using the start and end codons (real or hypothetical)

            my $first_codon = substr($orf, 0, 3);
            my $last_codon = substr($orf, (length $orf)-3, (length $orf) );

            my $extra = ".";
            if ($first_codon!~m/($selected_start_codons)/i){
               $extra = "${extra}hypothetical_start.";

               # call to a function that looks for an ATG above that is in frame with that hypothetical uORF

            };

            if ($last_codon!~m/($selected_stop_codons)/i ){
               $extra = "${extra}hypothetical_stop.";
               $length = $length + 3; # I have not tested this
               $aalength = $aalength + 1; # I have not tested this
            };
#            print "$first_codon\t$last_codon\t$extra\n";

            ### END label the transcript using the start and end codons (real or hypothetical)




            ### START label the uORF overlapping or upstream comparing the end of the uORF and the beginning of the CDS

            if ($updated_e > $cds_start_int){
               $uORFtype = "overlapping_uORF";
            } else {
               $uORFtype = "upstream_uORF";
            };

            ### END label the uORF overlapping or upstream comparing the end of the uORF and the beginning of the CDS




            ### START transform cDNA coordinates to genomic coordinates uORF start and uORF end

            # call to a function that transforms start and end
            # coordinates to genomic coordinates

            # these two variables will store the positions in the genomic coodinates
            my $s_general = $s;
            my $e_general = $updated_e;

            # we call the method cdna2genomic using the start and end in cDNA coordinates
            my @coords_gen = $trmapper->cdna2genomic($s, $updated_e);

            if (scalar(@coords_gen)){
               if ($strand_input == 1){
                  $s_general = $coords_gen[0]->start();
                  $e_general = $coords_gen[-1]->end();
               } else {
                  $s_general = $coords_gen[-1]->start();
                  $e_general = $coords_gen[0]->end();
               };
            };
                        # # negative strand CDS coordinates using cdna2genomic method
                        # X:100635558-100635569:-1
                        # X:100635178-100635252:-1
                        # X:100633931-100634029:-1
                        # X:100633405-100633539:-1
                        # X:100632485-100632568:-1
                        # X:100632063-100632068:-1

                        # # positive strand CDS coordinates using cdna2genomic method
                        # 20:417820-417930:1
                        # 20:419347-419468:1
                        # 20:419558-419731:1
                        # 20:422127-422238:1
                        # 20:427313-427492:1
                        # 20:428491-428589:1
                        # 20:428951-429094:1
                        # 20:430350-430430:1

            ### END transform cDNA coordinates to genomic coordinates to compute uORF start and uORF end




            ### START label the split start and/or stop codons and print them

            my $uORFoutput = "";


            # these two variables will store the positions in the genomic coodinates
            my $sc_s_new = $s;
            my $sc_e_new = $s+2;

            # we call the method cdna2genomic using the start and end in cDNA coordinates
            my @sc_coords = $trmapper->cdna2genomic($sc_s_new, $sc_e_new);



            # these two variables will store the positions in the genomic coodinates
            my $st_s_new = $updated_e-2;
            my $st_e_new = $updated_e;

            # we call the method cdna2genomic using the start and end in cDNA coordinates
            my @st_coords = $trmapper->cdna2genomic($st_s_new, $st_e_new);



            if ($extra!~m/(hypothetical_start)/i ){
               # not sure what will happen when the start codon appears splitted
               if (scalar(@sc_coords) > 1){
                  $extra = "${extra}split_start.";
#                  print "$extra\n";
               };
            };


            if ($extra!~m/(hypothetical_stop)/i ){
               # not sure what will happen when the stop codon appears splitted
               if (scalar(@st_coords) > 1){
                  $extra = "${extra}split_stop.";
#                  print "$extra\n";
               };
            };



#            if ($extra!~m/(hypothetical_start)/i ){ # this way all the start codons will be printed
               if (scalar(@sc_coords)){
                  foreach my $sc_cooo (@sc_coords){
                     $sc_s_new = $sc_cooo->start();
                     $sc_e_new = $sc_cooo->end();
#                     $uORFoutput = "${uORFoutput}chr$name_tran\t$uORFtype\tstart_codon\t$sc_s_new\t$sc_e_new\t$strand\tuORF-${transcript_id}-${uORF_in_trans};chr$name_tran:$s_general-$e_general:$strand\t$gene_id\t$transcript_id\t$length\t$aalength\t$extra\t$orf\t$prot_seq\n";
                     $uORFoutput = "${uORFoutput}chr$name_tran\t$uORFtype\tstart_codon\t$sc_s_new\t$sc_e_new\t$strand\tuORF-$first_codon-${transcript_id}-${uORF_in_trans};chr$name_tran:$s_general-$e_general:$strand\t$gene_id\t$transcript_id\t$length\t$aalength\t$extra\t$orf\t$prot_seq\n";

                  };
               };
#            };



#            if ($extra!~m/(hypothetical_stop)/i ){ # this way all the start codons will be printed
               if (scalar(@st_coords)){
                  foreach my $st_cooo (@st_coords){
                     $st_s_new = $st_cooo->start();
                     $st_e_new = $st_cooo->end();
                     $uORFoutput = "${uORFoutput}chr$name_tran\t$uORFtype\tstop_codon\t$st_s_new\t$st_e_new\t$strand\tuORF-$first_codon-${transcript_id}-${uORF_in_trans};chr$name_tran:$s_general-$e_general:$strand\t$gene_id\t$transcript_id\t$length\t$aalength\t$extra\t$orf\t$prot_seq\n";
                  };
               };
#            };



            ### END label the split start and/or stop codons and print them




            ### START transform cDNA coordinates to genomic coordinates

            # call to a function that transforms start and end
            # coordinates to genomic coordinates

            # these two variables will store the positions in the genomic coodinates
            my $s_new = $s;
            my $e_new = $updated_e;

            # we call the method cdna2genomic using the start and end in cDNA coordinates
            my @coords = $trmapper->cdna2genomic($s, $updated_e);
#            print "scalar(@coords)\n";


            if (scalar(@coords)){
#               my $uORFoutput = "";
               foreach my $cooo (@coords){
                     $s_new = $cooo->start();
                     $e_new = $cooo->end();
                     $uORFoutput = "${uORFoutput}chr$name_tran\t$uORFtype\texon\t$s_new\t$e_new\t$strand\tuORF-$first_codon-${transcript_id}-${uORF_in_trans};chr$name_tran:$s_general-$e_general:$strand\t$gene_id\t$transcript_id\t$length\t$aalength\t$extra\t$orf\t$prot_seq\n";
               };
#               $output = "${output}\n${uORFoutput}";
               $output = "${output}${uORFoutput}";
            };
                        # # negative strand CDS coordinates using cdna2genomic method
                        # X:100635558-100635569:-1
                        # X:100635178-100635252:-1
                        # X:100633931-100634029:-1
                        # X:100633405-100633539:-1
                        # X:100632485-100632568:-1
                        # X:100632063-100632068:-1

                        # # positive strand CDS coordinates using cdna2genomic method
                        # 20:417820-417930:1
                        # 20:419347-419468:1
                        # 20:419558-419731:1
                        # 20:422127-422238:1
                        # 20:427313-427492:1
                        # 20:428491-428589:1
                        # 20:428951-429094:1
                        # 20:430350-430430:1

            ### END transform cDNA coordinates to genomic coordinates




            ### START print the output of the current uORF to the variable we use to store it

            # print {$output_file} "chr$name_tran\t$uORFtype\t$s_new\t$e_new\t$strand\t$gene_id\t$transcript_id\t$length\t$aalength\t$extra\t$orf\t",Bio::Seq->new(-seq=>$orf)->translate->seq,"\n";


#            my $uORFoutput = "chr$name_tran\t$uORFtype\t$s_new\t$e_new\t$strand\t$gene_id\t$transcript_id\t$length\t$aalength\t$extra\t$orf\t$prot_seq";

#            $output = "${output}${uORFoutput}\n";

            ### END print the output of the current uORF to the variable we use to store it




            ### START other printing options

#            if ($strand_input == 1){
#               print {$output_file} "chr$name_tran\t$uORFtype\t$s:$s_new\t$updated_e:$e_new\t$strand\t$gene_id\t$transcript_id\t$length\t$aalength\t$extra\t$orf\t$prot_seq\n";
#            } else {
#               print {$output_file} "chr$name_tran\t$uORFtype\t$updated_e:$s_new\t$s:$e_new\t$strand\t$gene_id\t$transcript_id\t$length\t$aalength\t$extra\t$orf\t$prot_seq\n";
#            };

#~           If you want to convert a file with this output format to the regular one
#~           without cDNA:genomic coordinates, you can use the following command
#~                awk -F '[\t:]' '{print $1"\t"$2"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}'

            ### END other printing options




            ### START stop looking for other ends, break the loop

            last;

            ### END stop looking for other ends, break the loop
         }
      }
   }

   return ($output, $uORF_in_trans); # return the whole output of that transcript
}

### END of the function











# This subroutine will be called once by each process we generate
# It receives input and output filenames as input, and generates a file 
# containing all the uORFs found using the IDs in the input file.

sub manage_input_output {

   ### START storing input and output files from the function attributes and opening them

#   my $internal_input_filename = $_[0];
#   my $internal_output_filename = $_[1];

#   open(my $internal_input_file, '<', $internal_input_filename) or die "Could not open file '$internal_input_filename' $!";
#   open(my $internal_output_file, '>', $internal_output_filename) or die "Could not open file '$internal_output_filename' $!";

   open(my $internal_input_file, '<', $_[0]) or die "Could not open file '$_[0]' $!";
   open(my $internal_output_file, '>', $_[1]) or die "Could not open file '$_[1]' $!";
   open(my $internal_output_file_stats, '>', $_[2]) or die "Could not open file '$_[2]' $!";


   ### END storing input and output files from the function attributes and opening them



   ### START loading registry of ensembl database

   #******# START CONNECTION

   my $registry = 'Bio::EnsEMBL::Registry';

   $registry->load_registry_from_db(
       -host => 'useastdb.ensembl.org', # alternatively 'ensembldb.ensembl.org'
       -user => 'anonymous',
   );

   ### END loading registry of ensembl database




   ### START create a gene adaptor to fetch the gene IDs
   print "Searching in $species species.\n";
   my $gene_adaptor = $registry->get_adaptor($species, "Core", "Gene");

   # my $gen = $gene_adaptor->fetch_by_stable_id('ENSG00000125826');

   ### END create a gene adaptor to fetch the gene IDs




   # A possible improvement would be to read IDs in chunks and fetch multiple genes at the same time
   # https://www.perlmonks.org/?node_id=263605

   ### START iterate over the list of all gene IDs, for each, fetch it, and get the transcripts

   # we create a counter to keep track of the amount of genes processed per each process
   # and be able to estimate how long it will take to finish
   my $counter = 0;
   while (my $row = <$internal_input_file>) {

      # Have the gene ID as a string and without unwanted characters
      chomp $row;

      # fetch the gene using the ID obtained from the row
      # and the gene adaptor we built outside
      my $gen = $gene_adaptor->fetch_by_stable_id($row);

      # get all transcripts
      my @gen_transcripts = @{ $gen->get_all_Transcripts };


      # get the number of transcripts in the gene
      # print "Number of transcripts:" , scalar(@gen_transcripts) , "\n";


   ### END iterate over the list of all gene IDs, for each, fetch it, and get the transcripts




   ### START loop that iterates the transcripts, and calls the uORFs explorer function once per each transcript
      my $total_uORFs_gene = 0;
      my $chunk_write="";
      foreach my $tran (@gen_transcripts){
         # here we filter only protein coding genes, then we can also filter the
         #input to have a better parallelization
         if ( $tran->translation()) {
   #         if ($tran->cdna_coding_start()){

               ### START call to the function

               (my $small_chunk, my $total_uORFs_trans) = &uORFs_explorer($tran);

               ### END call to the function



               ### START increase uORFs per gene

               $total_uORFs_gene = $total_uORFs_gene + $total_uORFs_trans;

               ### END increase uORFs per gene



               ### START packing the output of the transcript with that of the gene

               # this condition is to avoid those transcripts that have no uORF
               unless ($small_chunk eq ""){
                  $chunk_write = "${chunk_write}${small_chunk}";
               }
               # probably it makes no sense, if there is nothing we don't care to attach it at the end

               ### END packing the output of the transcript with that of the gene


   #         print "\n\n--\n\n";
   #         print "\n";
            }
         }
      print {$internal_output_file} $chunk_write;
      $counter += 1;
      print {$internal_output_file_stats} $gen->display_id,"\t$total_uORFs_gene\n";
      print $gen->display_id,"\t$counter\n";

   } # end while

   ### END loop that iterates the transcripts, and calls the uORFs explorer function once per each transcript




   ### START closing input and output files and connections

   # Remove the connection to the Ensembl database server
   $registry->clear();

   #******# END CONNECTION



   # close the input and output files
   close $internal_input_file;
   close $internal_output_file;
   close $internal_output_file_stats;

   ### END closing input and output files and connections




   return;
}









### START run the program using "max_processes" number of processes

if (defined $max_processes){
#   open(my $input_file, '<', $input_filename) or die "Could not open file '$input_filename' $!";
#   open(my $output_file, '>', $output_filename) or die "Could not open file '$output_filename' $!";

   # Output header
#   print {$output_file} "chr\tuORFtype\tstart\tend\tstrand\tgene_id\ttranscript_id\tnuc_len\taa_len\textra\tcdna_seq\tprot_seq\n";
#   close $output_file;



   ### START work with files randomization, count lines and compute chunk size

   system "mv $input_filename $input_filename.ref; sort -R $input_filename.ref > $input_filename";

   my $total_lines = `wc -l < $input_filename`;
   chomp $total_lines;
   my $chunk_size = int($total_lines/$max_processes);
#   print "$total_lines lines divided in pieces of $chunk_size\n";
#   print "This is the total number of processes $max_processes\n";

   ### END work with files randomization, count lines and compute chunk size




   ### START generate the process manager and an array with the number of processes to generate

   my $manager = Parallel::ForkManager->new($max_processes);
   my @processes = (1..$max_processes);

   ### END generate the process manager and an array with the number of processes to generate




   ### START generate all the processes that will split the work

   foreach my $entry (@processes) {
      $manager->start and next;
      print "Process $entry says Hi!\n";

   ### END generate all the processes that will split the work


      ### START preparing the input file for each process

      my $temp_input_filename = "${input_filename}${entry}";

      unless ($entry == $max_processes){
         my $head_number = ($entry) * $chunk_size;
         my $tail_number = $chunk_size;
         system "head -n $head_number $input_filename | tail -n $tail_number > $temp_input_filename";
      } else {
         my $Dtail_number = $total_lines - (($entry-1) * $chunk_size);
         system "tail -n $Dtail_number $input_filename > $temp_input_filename";
      }

      ### END preparing the input file for each process



      ### START define the output file for each process and run the subroutine

      my $temp_output_filename = "${output_filename}${entry}";
      my $temp_output_stats = "${output_stats}${entry}";

      manage_input_output($temp_input_filename,
                          $temp_output_filename,
                          $temp_output_stats); 

      ### END define the output file for each process and run the subroutine



      ### START remove the input file for that process, and kill the process

      system "rm $temp_input_filename";
      print "Process $entry says Bye!\n";

      $manager->finish;

      ### END remove the input file for that process, and kill the process

   }

   # wait until all the processes have generated the output
   $manager->wait_all_children;
   print "All processes finished\n";


   ### START merge all the output files into a single one, and remove them once used

   # here we could try to merge all the generated files into a single one
   foreach my $child (@processes) {
      my $output_child = "${output_filename}${child}";
      my $output_child_stats = "${output_stats}${child}";
      system "cat $output_child >> $output_filename; rm $output_child";
      system "cat $output_child_stats >> $output_stats; rm $output_child_stats";
   }
   print "Output merged into $output_filename\n";
   print "Count stats merged into $output_stats\n";
   system "mv $input_filename.ref $input_filename"; # return the input file to its original order, not random

   ### END merge all the output files into a single one, and remove them once used


### END run the program using "max_processes" number of processes




} else {

   ### START run the program with only one process, call the corresponding function once

   print "$output_filename file CREATED!\n";
   manage_input_output($input_filename,
                       $output_filename,
                       $output_stats);

   ### END run the program with only one process, call the corresponding function once
}




### START keep track of the execution time

my $duration0 = time - $time_start0;

my $time_filename="/home/fcalvet/Desktop/uORFs/time/time";
open(my $time_file, '>>', $time_filename) or die "Could not open file $time_filename !";
print {$time_file} "Execution time ($input_filename) : $duration0 s\n";
close $time_file;

### END keep track of the execution time




__END__


Author: Ferriol Calvet
ferriol[at]ebi.ac.uk
Vertebrate Genomics team

# Core ORF finder extracted from:
# longorf.pl v0208020920
# (c) Dan Kortschak 2002
https://metacpan.org/release/BioPerl/source/examples/longorf.pl

https://stackoverflow.com/questions/43463831/how-to-run-a-perl-script-in-parallel


**** USAGE ****
# Run in parallel, recommended for lots of IDs
# take into account the memory limits of your computer, this could kill some processes and make that some genes are not processed
# always check your error file

# it is important to take memory restrictions into account

perl global_splitted_parallel_fromfile_ProperFormat_uORFs.pl -i ../../DATA/input_data/all_human_EnsemblGeneIDs_v34.txt -of ../../DATA/raw_data/raw_uORFs/allENSG_15-5-2020_at_17-12.tsv -m 100 -sp Human 2> error

OR

# Run sequentially when trying if the changes in the code work.
perl global_splitted_parallel_fromfile_ProperFormat_uORFs.pl -i ../DATA/input_data/allEnsemblGeneIDs_withoutPAR.txt -of ../DATA/raw_data/raw_uORFs/allENSG_27-2-2020_at_11-21.tsv 2> error

**** CREATE INPUT ****
1.- Obtain the list fo genes from BioMart (you need to remove the first line as it is a header

2.- Obtain the list of genes from the GENCODE gff3
gencode_file="/home/fcalvet/Desktop/uORFs/REFERENCE_DATA/gencode.v34.annotation.gff3";
output="/home/fcalvet/Desktop/uORFs/DATA/input_data/all_human_EnsemblGeneIDs_v34.txt";
grep -w gene $gencode_file | cut -d ';' -f2 | cut -d '=' -f2 | sort -u > $output;
# grep -w gene gencode.vM24.annotation.gff3 | cut -d ';' -f2 | cut -d '=' -f2 | sort -u > uORFs/DATA/input_data/all_mouse_EnsemblGeneIDs.txt





**** TOOLS FOR COMPARISONS ****

# whether two files are different or not
cmp --silent p1 p2 || echo "files are different"

# obtain the lines that are different between two files
diff --new-line-format="" --unchanged-line-format=""  file1 file2

# get the duplicated lines of a file
uniq -d file.txt




# NOTE

Whenever an error occurs and some process has died before finishing all the work, the best way to proceed is
wait until all the processes have finished, and then compare the genes in the stats file with the input file.
Then, we can see which genes have not been processed, and run again the script with those genes, concatenate
the input and follow with the other steps of the pipeline.





# When one connection dies, that process cannot finish the work.
# I can see this because the input file is not deleted, and once the error is found, no other ENSG IDs are processed.

# to check which are the missing 
diff <(sort allENSG_25-2-2020_at_15-57.tsv) <(sort allENSG_27-2-2020_at_11-21.tsv) | cut -f 6 | grep ENSG | sort -u > ../../input_data/not_present

grep -f ./input_data/not_present process_input_file.txt.pid | wc -l
wc -l ./input_data/not_present 
# compare the output of both commands


