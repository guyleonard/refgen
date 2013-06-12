#!/usr/bin/perl -w
# Author: Guy Leonard, Copyright MMVIII-MMXI
# Date: 2011

# Complete re-write!

# 1.beta.5 - trying to rewrite it so that each recognised set of seqs (e.g. NCBI or UNK) are output into their own table and can be hidden/toggled for viewing!
#		so far I have the javascript sorted for one case (i.e. UNK) i just need it to be selective
#		then I will need each dataset to be output in to tables but that will probably require ordering of arrays
#		what if people want to preserve order(will they? alignment destroys order anyway!)

# Import packages
use strict;
use warnings;
set_message(
"Please send any errors messages with a description of what you were trying to do and any relevant files to guy.leonard \@ gmail.com"
);
use diagnostics;
use Cwd;
use File::Path;
use File::Basename;
use CGI qw(:cgi-lib);
use CGI::Carp qw(fatalsToBrowser set_message);
use ExtUtils::Installed;
use feature qw{ switch };
use Socket;
use Bio::AlignIO;
use Bio::SeqIO;
use Digest::MD5 'md5_hex';

### User Editable Values
our $UPLOAD_DIR  = "\/var\/www";
our $REFGEN_PATH = "$UPLOAD_DIR\/refgen";
$CGI::POST_MAX = 1024 * 51200;    # Max Upload (1024 * 51200 = 50MB)

# Bunch of GLOBAL Variables
our $VERSION       = "1.beta.5";
our $EMPTY         = q{};
our $ID            = $EMPTY;
our $FILE          = $EMPTY;
our $EXT           = $EMPTY;
our $DIR           = $EMPTY;
our $CONTENT       = $EMPTY;
our @DEFLINE_ARRAY = $EMPTY;
our @DNA_AA_ARRAY  = $EMPTY;
our $DUPLICATES    = 0;

####
# Workflow - MAIN
our $BIOPERL = &check_bioperl();

my $query = CGI->new;

# Figure out which cgi web form we are retrieving
my $stage = $query->param("stage");
if ( $stage eq "one" ) {

    &session_id;
    &upload;
    &process_initial_form;

}
elsif ( $stage eq "two" ) {
    &process_refgen_form;

    # Usage goes here, as until the user clicks through to the final page
    # I don't consider the tool to have been used.
    &usage;
}

&tidy_up;

# End
####

# Check for BIOPERL on local server
# Was useful for legacy support - BioPerl is now REQUIRED
sub check_bioperl {
    my $Inst    = ExtUtils::Installed->new();
    my @Modules = $Inst->modules();
    my $flag    = "FALSE";
    foreach my $mod ( @Modules ) {
        if ( $mod eq "Bio" ) {
            $flag = "TRUE";
        }
    }

    return "$flag";
}

# Assign unique session ID to user each run
sub session_id {

# Get current epoch time - used as a lazy and "unique" identifier
# Perhaps I could hash? it in future - unlikely two people click exactly at the same time.
    $ID = time();
    mkdir( "$REFGEN_PATH\/$ID", 0755 )
      or die "Cannot create directory: $REFGEN_PATH\/$ID\n";
    return $ID;
}

sub upload {

    my $query = CGI->new;

    # Filehandle for uploaded file
    my $upload_file_name = $query->param( "file" );
    if ( !$upload_file_name ) {

        die( $query->header()
              . "No file selected, please choose a file to upload.\n" );
    }

    # Get the MIME type for the upload file
    my $type = $query->uploadInfo( $upload_file_name )->{ 'Content-Type' };
    ## This was causing far too many problems - disabled for now
#
# All uploads should be of text/plain, however without a .txt extension
# files are not recognised as plain text and so are assigned application/octet-stream.
# Other files are application/octet-stream and so we need to only allow a subset of
# these based on other acceptable extensions e.g. *.fas and *.fasta
# unless ( $type eq 'text/plain'
#     or $type eq 'application/octet-stream'
#     or $type eq 'application/x-extension-fas' )
# {
#     die "Only text/plain files are allowed\nYour type = $type\n";
# }

    ( $FILE, $DIR, $EXT ) = fileparse( $upload_file_name, qr{\..*} );

    # Check for tainting...
    # Convert any spaces to underscores "_"
    $FILE =~ tr/ /_/;
    $FILE =~ tr/./_/;
    $EXT  =~ tr/ /_/;

    # Check for illegal characters
    my $allowed_filename_characters = 'a-zA-Z0-9_.-';
    $FILE =~ s/[^$allowed_filename_characters]//g;
    $EXT  =~ s/[^$allowed_filename_characters]//g;
    if ( $FILE =~ /^([$allowed_filename_characters]+)$/ ) {
        $FILE = $1;
    }
    else {
        die(
"The filename is not valid. Filenames can only contain these characters: $allowed_filename_characters\n"
        );
    }

    # This isn't perfect as you could rename any file .fas or .fasta
    # and it would become application/octet-stream hmmmmmm
    if (   $EXT eq ".fa"
        or $EXT eq ".faa"
        or $EXT eq ".fna"
        or $EXT eq ".fas"
        or $EXT eq ".fasta"
        or $EXT eq ".fst"
        or $EXT eq ".txt" )
    {

        # Open output file
        open
          my $input_file, '>', "$REFGEN_PATH\/$ID\/$FILE$EXT"
          or die(
"Error 1: Unable to open file $ID\/$FILE$EXT: $!\nPlease change your file to have a valid fasta extension, e.g. .fasta or .fna\n"
          );
        binmode $input_file;

        # Print the contents of the uploaded file
        while ( <$upload_file_name> ) {
            my $count = length $_;
            $CONTENT = $CONTENT + $count;
            print $input_file $_;
        }
        close( $input_file );
    }
    else {
        die(
"File type (extension) must be .fasta or .fas etc, please try again.\n"
        );
    }

    &convertEOL( "$REFGEN_PATH\/$ID\/$FILE$EXT" );
    &remove_duplicate_sequences( "$REFGEN_PATH\/$ID\/$FILE$EXT\_EOLN" );
    unlink( "$REFGEN_PATH\/$ID\/$FILE$EXT\_EOLN" );

    # Split here to KEEP duplicates - WHY!?
    &file2array( "$REFGEN_PATH\/$ID\/$FILE\_removed_duplicates$EXT" );

}

# Conversion of EOL from input to native form.
sub convertEOL {
    my $file_name = shift;

    open
      my $input_file, '<', "$file_name"
      or die( "Error 2: Unable to read file $file_name: $!\n" );
    open
      my $output_file, '>', "$file_name\_EOLN"
      or die( "Error 3: Unable to save file $file_name\_EOLN: $!\n" );

    #Convert EOL from \r to \n
    while ( <$input_file> ) {
        my ( $line ) = $_;

        $line =~ s/(\015?\012|\015)/\n/gs;

        print $output_file "$line";
    }
    close( $input_file );
    close( $output_file );
}

sub remove_duplicate_sequences {

    my $file_name = shift;

    my %digests;
    my $in = Bio::SeqIO->new(
        - file   => "$file_name",
        - format => 'fasta'
    );

    my $out = Bio::SeqIO->new(
        - file   => ">$REFGEN_PATH\/$ID\/$FILE\_removed_duplicates$EXT",
        - format => 'fasta'
    );
    my $out2 = Bio::SeqIO->new(
        - file   => ">$REFGEN_PATH\/$ID\/duplicates$EXT",
        - format => 'fasta'
    );

    while ( my $seq = $in->next_seq ) {

        $digests{ md5_hex( $seq->seq ) }++;

        if ( $digests{ md5_hex( $seq->seq ) } == 1 ) {

            #print $seq->display_id . "\n" . $seq->seq . "\n";
            $out->write_seq( $seq );
        }
        else {

         #print "\tDUPE " . $seq->display_id . "\n\tDUPE " . $seq->seq . "\n\n";
            $out2->write_seq( $seq );
            $DUPLICATES++;
        }
    }

}

# Read the data file into an array, allows for more efficient handling
sub file2array {
    my $input_file = shift;

    ###
    # SPLIT METHOD HERE TO USE BIOPERL
    ###

    ## Read in the deflines
    open
      my $in_defline, '<', "$input_file"
      or die( "Error 4: Unable to read file $input_file: $!\n" );

    while ( my $line = <$in_defline> ) {
        chomp( $line );
        if ( $line =~ m/^>/ ) {
            push( @DEFLINE_ARRAY, $line );
        }
    }
    close( $in_defline );

    # Read in the sequences, line by line concatening the strings to one
    open
      my $in_sequence, '<', "$input_file"
      or die( "Error 5: Unable to read file $input_file: $!\n" );

    #
    my $temp = $EMPTY;
    while ( my $line = <$in_sequence> ) {
        chomp( $line );
        if ( $line !~ m/^>/ ) {
            $temp = "$temp$line";
        }
        else {
            push( @DNA_AA_ARRAY, $temp );
            $temp = "";
        }
    }
    push( @DNA_AA_ARRAY, $temp );
    close( $in_sequence );

    # Remove empty first element
    @DNA_AA_ARRAY  = splice( @DNA_AA_ARRAY,  2 );
    @DEFLINE_ARRAY = splice( @DEFLINE_ARRAY, 1 );

    #unlink("$input_file");
}

sub process_refgen_form {

    my $query = CGI->new;

    my $ID       = $query->param( "id" );
    my $num_seqs = $query->param( "num_seqs" );
    my $filename = $query->param( "filename" );

    my @taxa        = $EMPTY;
    my @accessions  = $EMPTY;
    my @refgen_ids  = $EMPTY;
    my @dnaaa_array = $EMPTY;

    # Extract the hidden TAXA from the form
    for ( my $i = 0 ; $i <= $num_seqs ; $i++ ) {

        my $taxon = $query->param( "hidden_taxa$i" );
        push( @taxa, $taxon );
    }

    # Extract the hidden ACCESSIONS from the form
    for ( my $i = 0 ; $i <= $num_seqs ; $i++ ) {

        my $accession = $query->param( "hidden_acce$i" );
        push( @accessions, $accession );
    }

    # Extract the REFGENIDS from the form
    for ( my $i = 0 ; $i <= $num_seqs ; $i++ ) {

        my $refgen_id = $query->param( "refgenid$i" );
        push( @refgen_ids, $refgen_id );
    }

    # Extract the hiddeen dna/aa sequences
    for ( my $k = 0 ; $k <= $num_seqs ; $k++ ) {

        my $dnaaa = $query->param( "dnaaa$k" );
        chomp( $dnaaa );
        push( @dnaaa_array, $dnaaa );
    }

    &refgen_html_out( \@refgen_ids, \@dnaaa_array, "$ID", "$filename" );
    &taxa_list_out( \@refgen_ids, \@accessions, \@taxa, "$ID", "$filename" );
    &tree_graph_out( \@refgen_ids, \@taxa, "$ID", "$filename" );

    if ( $BIOPERL == "TRUE" ) {

        &refgen_file_out( \@refgen_ids, \@dnaaa_array, "$ID", "$filename" );
        my ( $FILE, $DIR, $EXT ) = fileparse( $filename, qr{\..*} );

        # NEXUS
        my $in = Bio::AlignIO->new(
            - file   => "$REFGEN_PATH\/$ID\/$FILE\_refgen$EXT",
            - format => 'fasta'
        );
        my $out = Bio::AlignIO->new(
            - file   => ">$REFGEN_PATH\/$ID\/$FILE\_refgen\.nex",
            - format => 'nexus'
        );
        while ( my $aln = $in->next_aln() ) {
            $out->write_aln( $aln );
        }

        # PHYLIP
        my $in = Bio::AlignIO->new(
            - file   => "$REFGEN_PATH\/$ID\/$FILE\_refgen$EXT",
            - format => 'fasta'
        );
        my $out = Bio::AlignIO->new(
            - file   => ">$REFGEN_PATH\/$ID\/$FILE\_refgen\.phy",
            - format => 'phylip'
        );
        while ( my $aln = $in->next_aln() ) {
            $out->write_aln( $aln );
        }

# Maybe bioperl will fix this!?
# MASE
# my $in  = Bio::AlignIO->new(-file => "$REFGEN_PATH\/$ID\/$FILE\_refgen$EXT" , -format => 'fasta');
# my $out = Bio::AlignIO->new(-file => ">$REFGEN_PATH\/$ID\/$FILE\_refgen\.mase", -format => 'mase');
# while ( my $aln = $in->next_aln() ) { $out->write_aln($aln); }

    }
    else {

        # LEGACY - Remove
        &refgen_file_out( \@refgen_ids, \@dnaaa_array, "$ID", "$filename" );
        &refgen_nexus_out( \@refgen_ids, \@dnaaa_array, "$ID", "$filename" );
    }

    return "$ID";
}

sub taxa_list_out {

    my @refgen_ids = @{ $_[ 0 ] };
    my @accessions = @{ $_[ 1 ] };
    my @taxa       = @{ $_[ 2 ] };
    my $ID         = $_[ 3 ];
    my $filename   = $_[ 4 ];
    ( $FILE, $DIR, $EXT ) = fileparse( $filename, qr{\..*} );

    open
      my $taxa_out, '>', "$REFGEN_PATH\/$ID\/$FILE\_species_list.csv"
      or die
"Error 6: Unable to save file $REFGEN_PATH\/$ID\/$FILE\_species_list.csv: $!\n";

    for ( my $j = 1 ; $j <= $#refgen_ids ; $j++ ) {

        print $taxa_out "$accessions[$j],$refgen_ids[$j],$taxa[$j],\n";
    }
    close( $taxa_out );
}

sub tree_graph_out {

    my @refgen_ids = @{ $_[ 0 ] };
    my @taxa       = @{ $_[ 1 ] };
    my $ID         = $_[ 2 ];
    my $filename   = $_[ 3 ];
    ( $FILE, $DIR, $EXT ) = fileparse( $filename, qr{\..*} );

    open
      my $tg_out, '>', "$REFGEN_PATH\/$ID\/$FILE\_treegraph_name_table.txt"
      or die
"Error 7: Unable to save file $REFGEN_PATH\/$ID\/$FILE\_treegraph_name_table.txt: $!\n";

    for ( my $j = 1 ; $j <= $#refgen_ids ; $j++ ) {

        print $tg_out "$refgen_ids[$j]\t$taxa[$j]\n";
    }
    close( $tg_out );
}

sub refgen_html_out {

    my @refgen_ids  = @{ $_[ 0 ] };
    my @dnaaa_array = @{ $_[ 1 ] };
    my $ID          = $_[ 2 ];
    my $filename    = $_[ 3 ];
    ( $FILE, $DIR, $EXT ) = fileparse( $filename, qr{\..*} );

    open
      my $html_out, '>', "$REFGEN_PATH\/$ID\/$ID.html"
      or die "Error 8: Unable to save file $REFGEN_PATH\/$ID\/$ID.html: $!\n";

    print &PrintHeader;
    print &redirect_header( $ID, $#refgen_ids );    ## REDIRECT!

    print $html_out &normal_header( "REFGEN \| Results for $filename" );
    print $html_out <<"ENDOFTEXT";
	<body>
	<div class=\"title\">
    <span class=\"ref\">REF</span><span class=\"gen\">GEN</span> <img src=
    \"/refgen/css/refgen_title2.png\" class=\"title\" alt=\"REFGEN Logo\" /><span class=
    \"subtitle\">REFGEN - REFormat GEne Names for Phylogenies</span>
  </div>
  <div class="container">
        <form  id="form1"
  name="form1">
	<fieldset>
                <legend>REFGEN ID</legend>
                <h1><span></span>$ID</h1>
                <p>Please bookmark <a href="/refgen/$ID/$ID.html">this</a> URL if you wish to review your results at another time.<br /><br />
                Your results will be stored and available for around one week.</p>
        </fieldset>
        <fieldset>
                <legend>Original Files</legend>
                <ol>
                <li><a href="/refgen/$ID/$filename">$filename</a></li>
                <li><a href="/refgen/$ID/$FILE\_removed_duplicates$EXT">$FILE\_removed_duplicates$EXT</a></li>
                <li><a href="/refgen/$ID/duplicates$EXT">duplicates$EXT</a> &raquo; Based on sequence identity only.</li>
                </ol>
        </fieldset>
        <fieldset>
                <legend>REFGEN Files</legend>
                <ol>
                <li><a href="/refgen/$ID/$FILE\_refgen$EXT">$FILE\_refgen$EXT</a></li>
                <li><a href="/refgen/$ID/$FILE\_species_list.csv">$FILE\_species_list.csv</a></li>
                </ol>
        </fieldset>
        <fieldset>
                <legend>Other Files</legend>
                <ol>
                <li><a href="http://en.wikipedia.org/wiki/Nexus_file" target="_blank">Nexus</a> &raquo; <a href="/refgen/$ID\/$FILE\_refgen.nex">$FILE\_refgen.nex</a></li>
                <li><a href="http://evolution.genetics.washington.edu/phylip/doc/sequence.html" target="_blank">Phylip</a> &raquo; <a href="/refgen/$ID\/$FILE\_refgen.phy">$FILE\_refgen.phy</a></li>
                <!--<li>MASE &raquo; <a href="/refgen/$ID\/$FILE\_refgen.mase">$FILE\_refgen.mase</a></li>-->
                <li><a href="http://treegraph.bioinfweb.info/" target="_blank">TreeGraph 2</a> &raquo; <a href="/refgen/$ID\/$FILE\_treegraph_name_table.txt">$FILE\_treegraph_name_table.txt</a></li>
                </ol>
        </fieldset>
        <fieldset>
                <legend>Information</legend>
                <ol>
                        <li>Thank you for using REFGEN $VERSION</li><li>&nbsp;</li>
                        <li>Other Tools: <a href="../../treenamer.html">TREENAMER</a> - A tool to change the unique identifiers made by REFGEN to a more 'human' readable format (accession & species name) in your tree files.</li><li>&nbsp;</li>
                        <li>Please come back again and don't forget to suggest these tools to your friends!</li><li>&nbsp;</li><li>Citation: <a href=
        "http://www.la-press.com/refgen-and-treenamer-automated-sequence-data-handling-for-phylogenetic-a1451">
        Leonard, G., et al. (2009). REFGEN and TREENAMER: Automated Sequence Data Handling for
        Phylogenetic Analysis in the Genomic Era. <em>Evolutionary Bioinformatics Online</em>,
        5:1-4.</a></li><li>&nbsp;</li>
                </ol>
        </fieldset>
        </form>
</div>
ENDOFTEXT
    print $html_out &HtmlBot;
    close( $html_out );
}

sub refgen_file_out {

    my @refgen_ids  = @{ $_[ 0 ] };
    my @dnaaa_array = @{ $_[ 1 ] };
    my $ID          = $_[ 2 ];
    my $filename    = $_[ 3 ];

    ( $FILE, $DIR, $EXT ) = fileparse( $filename, qr{\..*} );

    open
      my $refgen_file_out, '>', "$REFGEN_PATH\/$ID\/$FILE\_refgen$EXT"
      or die
      "Error 9: Unable to save file $REFGEN_PATH\/$ID\/$FILE\_refgen$EXT: $!\n";

    for ( my $j = 1 ; $j <= $#refgen_ids ; $j++ ) {
        my @pieces = $dnaaa_array[ $j ] =~ m/(.{1,70})/gs;
        print $refgen_file_out ">$refgen_ids[$j]\n";

        foreach ( @pieces ) {
            print $refgen_file_out "$_\n";
        }
        print $refgen_file_out "\n";
    }
    close( $refgen_file_out );
}

sub refgen_nexus_out {

    my @refgen_ids  = @{ $_[ 0 ] };
    my @dnaaa_array = @{ $_[ 1 ] };
    my $ID          = $_[ 2 ];
    my $filename    = $_[ 3 ];

    ( $FILE, $DIR, $EXT ) = fileparse( $filename, qr{\..*} );

    ### NEXUS File
    my $maxLength = 0;
    my $datatype  = "DNA";

    ## Get the largest sequence length for NEXUS nchar value and check for PROTEIN or DNA
    ## Idealy this should only be done if all sequence lengths are the same - i.e. when I can be bothered to
    ## add a checkbox or some form of checking system, for now it just happens.
    for ( my $i = 0 ; $i <= $#refgen_ids ; $i++ ) {
        my $seq_length = length( $dnaaa_array[ $i ] );
        if ( $seq_length ge $maxLength ) {
            $maxLength = $seq_length;
        }
    }

    for ( my $i = 0 ; $i <= $#refgen_ids ; $i++ ) {
        my $code = $dnaaa_array[ $i ];
        if ( $code =~ m/[^BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvxxyz]/i ) {
            $datatype = "PROTEIN";
        }
        my $seq_length = length( $dnaaa_array[ $i ] );
        if ( $seq_length ge $maxLength ) {
            $maxLength = $seq_length;
        }

    }

    open
      my $nexus_output, '>', "$REFGEN_PATH\/$ID\/$FILE\_refgen.nex"
      or die(
"Error 10: Unable to open file $REFGEN_PATH\/$ID\/$FILE\_refgen.nex: $!\n"
      );

    my $ntax = $#refgen_ids;

    print $nexus_output "#NEXUS\n";
    print $nexus_output "[File created by REFGEN Version $VERSION on "
      . localtime() . "]\n";
    print $nexus_output "BEGIN DATA;\n";
    print $nexus_output "  DIMENSIONS NTAX=$ntax NCHAR=$maxLength;\n";
    print $nexus_output "  FORMAT DATATYPE=$datatype\n  ;\n";
    print $nexus_output "MATRIX\n";
    my $count = 1;

    for ( my $i = 1 ; $i <= $#refgen_ids ; $i++ ) {
        print $nexus_output "\[$count\] $refgen_ids[$i]\n";
        my @pieces = $dnaaa_array[ $i ] =~ m/(.{1,70})/gs;
        foreach ( @pieces ) {
            print $nexus_output "$_\n";
        }
        print $nexus_output "\n";
        $count++;
    }
    print $nexus_output "END;\n";
    close( $nexus_output );
}

sub process_initial_form {

    my $query = CGI->new;

    # Get accession length from cgi web form
    my $accession_length = $query->param( "acclen" );

    # Get direction to chop string?
    my $acc_head_tail = $query->param( "head_tail" );

    # Get Species Name length from cgi web form
    my $sp_length_0 = $query->param( "species" );
    my $sp_length_1 = $query->param( "name" );

    # Get separator from cgi webform
    my $separator = $query->param( "sep" );

    my $unknown_db = $EMPTY;

    # Print out form for user-validation
    print &PrintHeader;
    print &normal_header( "REFGEN \| Semi-Automatic Mode" );
    my $number_of_seqs  = @DNA_AA_ARRAY;
    my $number_of_seqs2 = @DEFLINE_ARRAY;
    print <<"ENDOFTEXT";
	<body>
	<script language="JavaScript" type="text/javascript">
<!--
// Disable and Enable Form Elements from a CheckBox
// onclick="f40_Disable(null,this,true);"
// true  = if the checkbox is checked the elements will be 'readonly';
// false = if the checkbox is checked the elements will be 'enabled';

function f40_Disable(f40_par,f40_obj,f40_state){
 if (f40_par){f40_clds=f40_AllElements(document.getElementById(f40_par)); }
 else { f40_clds=f40_AllElements(f40_obj.parentNode); }
 if (!f40_obj.ary){
  f40_obj.ary=new Array();
  for (f40_0=0;f40_0<f40_clds.length;f40_0++){
   if (f40_clds[f40_0].tagName=='INPUT'||f40_clds[f40_0].tagName=='SELECT'||f40_clds[f40_0].tagName=='TEXTAREA'){
    f40_obj.ary[f40_obj.ary.length]=f40_clds[f40_0];
   }
  }
 }
 for (f40_1=0;f40_1<f40_obj.ary.length;f40_1++){
  f40_obj.ary[f40_1].removeAttribute('readonly');
 }
 if (f40_obj.checked==f40_state){
  for (f40_2=0;f40_2<f40_obj.ary.length;f40_2++){
    f40_obj.ary[f40_2].setAttribute('readonly','readonly');
  }
 }
 f40_obj.removeAttribute('readonly');
}

function f40_AllElements(f40_){
  if (f40_.all){ return f40_.all; }
  return f40_.getElementsByTagName('*');
}

function change_refgen_id(z,x,y) {
	var refgenid_new = z
	var acta = x
	var taac = y

	var re1='(acce)'
	var p1 = new RegExp(re1,["i"])
	var m1 = p1.exec(acta)

	var re2='(taxa)'
	var p2 = new RegExp(re2,["i"])
	var m2 = p2.exec(acta)

	// Need to add an if to differentiate between taxa and acce
	if (m1 != null) {
	
		var acce = document.getElementById(acta).value
		var taxa = document.getElementById(taac).value
		document.getElementById(refgenid_new).value=acce + taxa
	}
	else if (m2 != null) {
		var acce = document.getElementById(taac).value
		var taxa = document.getElementById(acta).value
		document.getElementById(refgenid_new).value=acce + taxa
	}	
}
//-->
</script>

<script type="text/javascript" charset="utf-8">
		function findElementsWithClass(tagName, className) {
		    if (document.querySelectorAll) {
		        return document.querySelectorAll(tagName + "." + className);
		    } else {
		        var results = [];
		        var all = document.getElementsByTagName(tagName);
		        var regex = new Regexp("(?:^|\\s)" + tagName + "(?:\\s|$)");
		        for (var i=0, len=all.length; i<len; ++i) {
		            if (regex.test(all[i].className)) {
		                results.push(all[i]);
		            }
		        }
		        return results;
		    }
		}

		window.onload = function() {
		    var visible = true;
		    document.getElementById('toggle').onclick = function() {
		        visible = !visible;
		        var tds = findElementsWithClass('tr', 'UNK');
		        for (var i=0, len=tds.length; i<len; ++i) {
		            tds[i].parentNode.style.display = visible ? '' : 'none';
		        }
		    };
		}
	</script>


          <div class=\"title\">
    <span class=\"ref\">REF</span><span class=\"gen\">GEN</span> <img src=
    \"../refgen/css/refgen_title2.png\" class=\"title\" alt=\"REFGEN Logo\" /><span class=
    \"subtitle\">REFGEN - REFormat GEne Names for Phylogenies</span>
  </div>
	<table class=\"center\">
	        <tr>
	            <td>Filename:</td>
	            <td>$FILE$EXT</td>
	        </tr>
	        <tr>
	            <td>REFGEN ID:</td>
	            <td><h1><span></span>$ID</h1></td>
	        </tr>
	        <tr>
	            <td>Bioperl:</td>
	            <td>$BIOPERL</td>
	            <td>Version:</td>
	            <td>$VERSION</td>
	        </tr>
	        
ENDOFTEXT

    # Time wasting error count
    for ( my $x = 0 ; $x <= $#DNA_AA_ARRAY ; $x++ ) {
        my ( $accession, $species_name, $accession_modified, $type,
            $species_name_modified, $errors )
          = &process_deflines(
            $accession_length, $acc_head_tail,       $sp_length_0,
            $sp_length_1,      $DEFLINE_ARRAY[ $x ], $unknown_db
          );
        if ( $errors eq 1 ) {
            $unknown_db++;
        }
    }
    print <<"ENDOFTEXT";
    	        <tr>
	            <td>Unknowns:</td>
	            <td>$unknown_db</td>
	            <td>Duplicates:</td>
	            <td>$DUPLICATES</td>
	        </tr> 
	</table>
	<br />
	<table class=\"center\">
	        <tr>
	                <td>Key:&nbsp;&nbsp;</td>
	                <td class=\"NCBI\">NCBI</td>
	                <td class=\"JGI\">JGI</td>
	                <td class=\"GEN\">GenBank</td>
	                <td class=\"SILV\">SILVA</td>
	                <td class=\"UNK\">Unknown</td>
	        </tr>
	</table>
ENDOFTEXT

    print
"<form enctype=\"multipart\/form-data\" method=\"post\" action=\"/cgi-bin/refgen_1_beta.pl\"><table class=\"center\"><tr><th>#</th><th>Original Input</th><th>Accession</th><th>Taxa</th><th>REFGEN ID</th></tr><tr><td></td><td></td><td></td><td>Only click here if you are sure there are no errors!</td><td><div class=\"right\"><input name=\"submit\" class=\"btn green\" value=\"REFGEN!\" type=\"submit\" /></div></td></tr>";

    # Reset error counter
    $unknown_db = 1;

    # Actual bit where I loop to get all values...
    for ( my $x = 0 ; $x <= $#DNA_AA_ARRAY ; $x++ ) {

        my ( $accession, $species_name, $accession_modified, $type,
            $species_name_modified, $errors )
          = &process_deflines(
            $accession_length, $acc_head_tail,       $sp_length_0,
            $sp_length_1,      $DEFLINE_ARRAY[ $x ], $unknown_db
          );
        if ( $type eq "UNK" ) {
            $unknown_db++;
            print <<"ENDOFTEXT";
                <input type="checkbox" id="toggle" /> Toggle<br />
                <tr class=\"$type\">
                <td class=\"number\">$x</td>
                <td><textarea name=\"seq$x\" cols=40 rows=4 readonly>$DEFLINE_ARRAY[$x]
		$DNA_AA_ARRAY[$x]</textarea>
        	</td>
		<td>$accession<br />
		
		<input name=\"hidden_acce$x\" type=\"text\" id=\"hidden_acce$x\" value=\"$accession\" readonly />
		<br />
		<input name=\"acce$x\" type=\"text\" id=\"acce$x\" value=\"$accession_modified\" onkeyup=\"change_refgen_id('refgenid$x', this.id, 'taxa$x')\" readonly />
		<input type=\"checkbox\" name=\"box$x\" onclick=\"f40_Disable(null,this,false);\"></td>

		<td>$species_name<br />
		
		<input name=\"hidden_taxa$x\" type=\"text\" id=\"hidden_taxa$x\" value=\"$species_name\" readonly />
		<br />
		<input name=\"taxa$x\" type=\"text\" id=\"taxa$x\" value=\"$species_name_modified\" onkeyup=\"change_refgen_id('refgenid$x', this.id, 'acce$x')\" readonly />
		<input type=\"checkbox\" name=\"box$x\" onclick=\"f40_Disable(null,this,false);\"></td>

		<td class=\"refgenid\">$species_name_modified$accession_modified<br />
		<input name=\"refgenid$x\" type=\"text\" id=\"refgenid$x\" value=\"$species_name_modified$accession_modified\" readonly />
		<input type=\"checkbox\" name=\"box$x\" onclick=\"f40_Disable(null,this,false);\">
		
		<!-- Hidden DNA/AA Sequence -->
		<input name=\"dnaaa$x\" type=\"hidden\" readonly value=\"$DNA_AA_ARRAY[$x]\"/>
		
		</td></tr>
ENDOFTEXT
        }
        else {

            print <<"ENDOFTEXT";
                <tr class=\"$type\">
                <td class=\"number\">$x</td>
                <td><textarea name=\"seq$x\" cols=40 rows=4 readonly>$DEFLINE_ARRAY[$x]
		$DNA_AA_ARRAY[$x]</textarea>
        	</td>
		<td>$accession<br />
		
		<input name=\"hidden_acce$x\" type=\"hidden\" id=\"hidden_acce$x\" value=\"$accession\" readonly />
		
		<input name=\"acce$x\" type=\"text\" id=\"acce$x\" value=\"$accession_modified\" onkeyup=\"change_refgen_id('refgenid$x', this.id, 'taxa$x')\" readonly />
		<input type=\"checkbox\" name=\"box$x\" onclick=\"f40_Disable(null,this,false);\"></td>

		<td>$species_name<br />
		
		<input name=\"hidden_taxa$x\" type=\"hidden\" id=\"hidden_taxa$x\" value=\"$species_name\" readonly />
		
		<input name=\"taxa$x\" type=\"text\" id=\"taxa$x\" value=\"$species_name_modified\" onkeyup=\"change_refgen_id('refgenid$x', this.id, 'acce$x')\" readonly />
		<input type=\"checkbox\" name=\"box$x\" onclick=\"f40_Disable(null,this,false);\"></td>

		<td class=\"refgenid\">$accession_modified$species_name_modified<br />
		<input name=\"refgenid$x\" type=\"text\" id=\"refgenid$x\" value=\"$accession_modified$species_name_modified\" readonly />
		<input type=\"checkbox\" name=\"box$x\" onclick=\"f40_Disable(null,this,false);\">
		
		<!-- Hidden DNA/AA Sequence -->
		<input name=\"dnaaa$x\" type=\"hidden\" readonly value=\"$DNA_AA_ARRAY[$x]\"/>
		
		</td></tr>
ENDOFTEXT
        }
    }

    print <<"ENDOFTEXT";
		        <tr>
			<td><input name=\"stage\" type=\"hidden\" readonly value=\"two\" /></td>
			<td><input name=\"id\" type=\"hidden\" readonly value=\"$ID\"/></td>
			<td><input name=\"filename\" type=\"hidden\" readonly value=\"$FILE$EXT\" /></td>
			<td><input name=\"num_seqs\" type=\"hidden\" readonly value=\"$#DNA_AA_ARRAY\" /></td>
			<td><div class=\"right\"><input name=\"submit\" class=\"btn green\" value=\"REFGEN!\" type=\"submit\" /></div></td>
			</tr>
		</table>
		</form>

ENDOFTEXT
    print &HtmlBot;
}

# Normal HTML Header
sub normal_header {

    my $title = shift;

    #
    return <<"END_OF_TEXT";
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
  <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
  <title>$title</title>
  <link rel="stylesheet" type="text/css" media="screen" href="/refgen/css/refgen.css" />
</head>
END_OF_TEXT
}

# Redirect HTML Header
sub redirect_header {

    my $ID = shift;
    $CONTENT = shift;

    # Do we do this here? Make it log?
    $CONTENT = $CONTENT / 1000;

    #
    return <<"END_OF_TEXT";
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
  <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
  <title>Your REFGEN Results are being generated...</title>
  <META HTTP-EQUIV="Refresh" CONTENT="$CONTENT; URL=../refgen/$ID/$ID.html">
  <link rel="stylesheet" type="text/css" media="screen" href="../refgen/css/refgen.css" />
</head>
END_OF_TEXT
}

# Update January 2011
# Takes one defline at a time and returns
# Updated November 2009
# More explict matching on deflines with expansion to include new types
# any that do not fit exact matches are dumped into a non-match array
sub process_deflines {

    my $accession_length = $_[ 0 ];
    my $acc_head_tail    = $_[ 1 ];
    my $sp_length_0      = $_[ 2 ];
    my $sp_length_1      = $_[ 3 ];
    my $defline          = $_[ 4 ];
    my $unk_count        = $_[ 5 ];

    #
    my $accession             = $EMPTY;
    my $species_name          = $EMPTY;
    my $non_gb_jgi            = $EMPTY;
    my $accession_modified    = $EMPTY;
    my $species_name_modified = $EMPTY;
    my $type                  = $EMPTY;
    my $errors                = $EMPTY;

    # Switch was deprecated - trying given/when
    # switch ($line)
    given ( $defline ) {

        ## Long gi defintion - Standard NCBI Output from BLASTp search
    # >gi|71649003|ref|XP_813263.1| hypothetical protein [Trypanosoma cruzi]
    # Accession is 3rd element
    # Species name is generally in square braces [] of 4th element for proteins.
        when ( /(>gi)(\|)(.*?)(\|)(.*?)(\|)(.*?)(\|)(.*?)/ ) {

            # Split the line by the pipe symbol, 5 in length
            my @current_defline =
              split( /\|/, $defline );

            # Accession Number to Array
            $accession = $current_defline[ 3 ];
            $accession_modified =
              &modify_accessions( $accession, $acc_head_tail,
                $accession_length );

            # Species Name to Array
            $species_name = "gil$current_defline[4]";
            ( $species_name, $species_name_modified ) =
              &modify_taxa_names( $species_name, $sp_length_0, $sp_length_1 );

            $type = "NCBI";
        }

        ## JGI non-scaffold
      # >jgi|Phyra1_1|71587
      # Accession is 2nd element after 2nd pipe|
      # Species name is already concatenated and is 1st element, after 1st pipe|
        when ( /(>jgi)(\|)(.*?)(\|)(.*?)/ ) {

            # Split the line by the pipe symbol, 5 in length
            my @current_defline =
              split( /\|/, $defline );

            # Accession Number to Array
            $accession = "$current_defline[2]";
            $accession_modified =
              &modify_accessions( $accession, $acc_head_tail,
                $accession_length );

            # Species Name to Array
            # Automatically change to Genus species from list
            my $jgi_proper = &get_jgi_species_name( "$current_defline[1]" );
            $species_name = "jgi$jgi_proper";
            ( $species_name, $species_name_modified ) =
              &modify_taxa_names( $species_name, $sp_length_0, $sp_length_1 );

            $type = "JGI";
        }

        ## JGI scaffold
        # >fgenesh1_pg.C_scaffold_1000025 [Phyra1_1:72583]
        # Accession is element after colon
        # Species name is element before colon
        when ( /(>fgene)(.*?)(\\[.*?\\])/ ) {

            # Split the line by colon
            my @current_defline =
              split( /\:/, $defline );

            # Accession Number to Array
            $accession = "$current_defline[1]";
            $accession_modified =
              &modify_accessions( $accession, $acc_head_tail,
                $accession_length );

            # Species Name to Array
            # Automatically change to Genus species from list
            my $jgi_proper = &get_jgi_species_name( "$current_defline[0]" );
            $species_name = "jgi$jgi_proper";
            ( $species_name, $species_name_modified ) =
              &modify_taxa_names( $species_name, $sp_length_0, $sp_length_1 );

            $type = "JGI";
        }

	## Yet another JGI format
	# >16474027_peptide|Bdistachyon|Bradi1g21050|Bradi1g21050.1
	# >16404070_peptide|Cpapaya|evm.TU.supercontig_1.153|evm.model.supercontig_1.153
	# >18152817_peptide|Acoerulea|AcoGoldSmith_v1.003707m.g|AcoGoldSmith_v1.003707m
	# Accession is first number before underscore
	# Species name is after 1st pipe - annoyingly in a different format. Not 3/2 but 1/Full
	when (/(>\d+\_)(.*?)(\|)(?:[a-z][a-z]+)(\|)(.*?)/) {
	
  	    # Split the line by the pipe symbol, ~4 in length
            my @current_defline = split( /\|/, $defline );

	    # Accession Number to Array
            $accession = "$current_defline[0]";
            $accession_modified = &modify_accessions( $accession, $acc_head_tail, $accession_length );

            # Species Name to Array
            # Automatically change to Genus species from list
            my $jgi_proper = &get_jgi_species_name( "$current_defline[1]" );
            $species_name = "jgi$jgi_proper";
            ( $species_name, $species_name_modified ) =
              &modify_taxa_names( $species_name, $sp_length_0, $sp_length_1 );

            $type = "JGI";
	
	}

        ## General database identifier
        # This covers gb, emb, dbk, sp and ref or any similar formats
        # gb|accession|locus
        # Accession is 1st element
        # Species name is potentially absent, use 2nd element
        when ( /(>.{2,3})(\|)(.*?)(\|)(.*?)/ ) {

            # Split the line by the pipe symbol, 3 in length
            my @current_defline =
              split( /\|/, $defline );

            # Accession Number to Array
            $accession = "$current_defline[1]";
            $accession_modified =
              &modify_accessions( $accession, $acc_head_tail,
                $accession_length );

            # Species Name to Array
            $species_name = "gba$current_defline[2]";
            ( $species_name, $species_name_modified ) =
              &modify_taxa_names( $species_name, $sp_length_0, $sp_length_1 );

            $type = "GEN";

        }

        ## SILVA DB DNA
        # >AB513180.1.1495 | Burkholderia acidipaludis
        # >AADB02002333.14401.16248 | Homo sapiens (human)
        when ( /(>)(.*?)(\s)(\|)(\s)(.*?)/ ) {

            # Split the line by the pipe symbol
            my @current_defline =
              split( /\|/, $defline );

            # Accession Number to Array
            my @current_accession =
              split( /\./, "$current_defline[0]" );
            $accession = "$current_accession[0]";
            $accession_modified =
              &modify_accessions( $accession, $acc_head_tail,
                $accession_length );

            # Species Name to Array - can use same method for gil
            $species_name = "gil\[$current_defline[1]\]";
            ( $species_name, $species_name_modified ) =
              &modify_taxa_names( $species_name, $sp_length_0, $sp_length_1 );

            $type = "SILV";
        }

        default {

   # I should probably properly catch the errors - we'll think abotu it for now.
            $accession             = "0000$unk_count";
            $accession_modified    = "$unk_count";
            $species_name          = "Unknown Sequence";
            $species_name_modified = "Unk_seq";
            $type                  = "UNK";
            $errors                = 1;
            $unk_count++;

        }
    }
    return (
        "$accession", "$species_name",          "$accession_modified",
        "$type",      "$species_name_modified", "$errors"
    );
}

sub get_jgi_species_name {

    my $jgi_accession = shift;

    my @jgi_name_array = $EMPTY;

    open
      my $in_jgi, '<', "$REFGEN_PATH\/doe_jgi_list.csv"
      or die(
"Error 11: Unable to read DOE JGI file $REFGEN_PATH\/doe_jgi_list.csv: $!\n"
      );
    while ( <$in_jgi> ) {
        my ( $line ) = $_;
        push( @jgi_name_array, $line );
    }
    close( $in_jgi );
    my $jgi_name_array_length = @jgi_name_array;

    for ( my $i = 0 ; $i < $jgi_name_array_length ; $i++ ) {

        my @temp =
          split( /,/, $jgi_name_array[ $i ] );
        my $jgi_acc   = $temp[ 0 ];
        my $jgi_locus = $temp[ 1 ];
        chomp( $jgi_locus );
        if ( $jgi_accession eq $jgi_acc ) {

            #my $jgi_proper = $jgi_locus;
            return "$jgi_locus";
        }

#else {
# If no match is found (occurs when JGI have new genomes that are not in my list)
# Return the initial scalar variable
#return "$jgi_accession";
#}
    }
}

sub modify_accessions {

    my $accession        = $_[ 0 ];
    my $acc_head_tail    = $_[ 1 ];
    my $accession_length = $_[ 2 ];

    # Assign value from accessionArray to short_accArray
    my $short_accession = $accession;

    # This should remove all non-alphanumeric chars
    $short_accession =~ tr/\W*\.*//d;

    # Create sub string of value...
    if ( $acc_head_tail eq "head" ) {
        $short_accession = substr( $short_accession, 0, $accession_length );
    }
    else {
        $short_accession = substr( $short_accession, -$accession_length );
    }

    return "$short_accession";
}

sub modify_taxa_names {

    my $species_name          = $_[ 0 ];
    my $sp_length_0           = $_[ 1 ];
    my $sp_length_1           = $_[ 2 ];
    my $species_name_modified = $EMPTY;

    # Remove any non-alphanumerics
    $species_name =~ tr/\W*\.*\:*//d;

    given ( $species_name ) {

        # DOE JGI Format
        when ( /^jgi/i ) {
            my $length = length $species_name;

            # Remove jgi tag...
            $species_name = substr( $species_name, 3, $length );

            my @species =
              split( / /, $species_name );
            my $sp   = substr( $species[ 0 ], 0, $sp_length_0 );
            my $name = substr( $species[ 1 ], 0, $sp_length_1 );
            $species_name_modified = "$sp$name";

        }

        # GenBank / EMBL Data Lib / DDBJ / Swiss-Prot / NCBI Reference Sequence
        when ( /^gba/i ) {
            my $length = length $species_name;

            # Remove gb tag...
            $species_name = substr( $species_name, 3, $length );

            my @species =
              split( / /, $species_name );
            my $sp   = substr( $species[ 0 ], 0, $sp_length_0 );
            my $name = substr( $species[ 1 ], 0, $sp_length_1 );
            $species_name_modified = "$sp$name";
        }

        # gi identifiers from NCBI
        when ( /^gil/i ) {

# There is an assumed preference for text in [] followed by () then whatever is left
# if $type contains square braces []
            if ( $species_name =~ m/(.*?)(\[.*?\])/is ) {
                $species_name = $2;

# if within the [] there are round braces (), ignore e.g. [Oryza sativa (japonica cultivar)]
                if ( $species_name =~ m/(.*?)(\(.*\))/is ) {
                    $species_name = $1;

                    # Trim all braces and spaces
                    $species_name =~ tr/\(\)\[\]/ /d;
                    $species_name =~ s/^\s+//;
                    $species_name =~ s/\s+$//;

                    # $species_name = $species;
                    my @species =
                      split( / /, $species_name );
                    my $sp   = substr( $species[ 0 ], 0, $sp_length_0 );
                    my $name = substr( $species[ 1 ], 0, $sp_length_1 );
                    $species_name_modified = "$sp$name";
                }
                else {
                    $species_name =~ tr/\[\]/ /d;
                    $species_name =~ s/^\s+//;
                    $species_name =~ s/\s+$//;

                    # $species_name = $species;
                    my @species =
                      split( / /, $species_name );
                    my $sp   = substr( $species[ 0 ], 0, $sp_length_0 );
                    my $name = substr( $species[ 1 ], 0, $sp_length_1 );
                    $species_name_modified = "$sp$name";
                }
            }

            # Else if it contains round braces
            elsif ( $species_name =~ m/(\(.*\))/is ) {
                $species_name = $1;

                # Remove braces and spaces
                $species_name =~ tr/\(\)/ /d;
                $species_name =~ s/^\s+//;
                $species_name =~ s/\s+$//;

                # $species_name = $species;
                my @species =
                  split( / /, $species_name );
                my $sp   = substr( $species[ 0 ], 0, $sp_length_0 );
                my $name = substr( $species[ 1 ], 0, $sp_length_1 );
                $species_name_modified = "$sp$name";
            }
            else {
                my $length = length $species_name;

                # Remove gil tag...
                $species_name = substr( $species_name, 3, $length );

                my @species =
                  split( / /, $species_name );
                my $sp   = substr( $species[ 0 ], 0, $sp_length_0 );
                my $name = substr( $species[ 1 ], 0, $sp_length_1 );
                $species_name_modified = "$sp$name";
            }
        }
    }
    return ( $species_name, $species_name_modified );
}

sub usage {

    my $user_IP      = $ENV{REMOTE_ADDR};
    my $user_browser = $ENV{HTTP_USER_AGENT};
    my $time         = localtime time;

    my $iaddr = inet_aton( "$user_IP" );
    my $dns_name = gethostbyaddr( $iaddr, AF_INET );

    open
      my $usage_out, '>>', "$REFGEN_PATH\/refgen_usage.txt"
      or die
      "Error 12: Unable to write to file $REFGEN_PATH\/refgen_usage.txt: $!\n";

    print $usage_out "$time\t$user_IP\t$dns_name\t$ID\t$user_browser\n";
    close( $usage_out );
}

sub tidy_up {

    # Clear up old file results
    my $time_now = time();
    my $d        = "$REFGEN_PATH";
    opendir( DIR, $d )
      or die( "Error 13: Cannot open directory $d: $!\n" );
    my @dirs =
      grep { !/^\./ && -d "$d/$_" } readdir( DIR );
    closedir DIR;
    for ( my $i = 0 ; $i <= $#dirs ; $i++ ) {
        my $diff = $time_now - $dirs[ $i ];

        # Delete after four weeks (in seconds)
        if ( $diff >= 2419200 ) {
            rmtree( "$d\/$dirs[$i]", 0, 0 );
        }
    }
}
