#!/usr/bin/perl
# Author: Guy Leonard, Copyright MMVIII-MMX
# Date: 2010-04-29

## Updates * Change/Fix, + Addition, - Deletion
# + 2010-01-19 DOE JGI List of species names lookup from CSV file.
# + 2010-01-06 Hardened upload security - check for tainting, MIME type and extension.
# * 17-11-2009 Fixed and rearranged regexs for ID of NCBI/JGI and new gb/emb/ddbj deflines
# * 17-11-2009 As above fixing an error where unidentified deflines caused DNA sequences to go out of sync
# + XX-XX-2009 Experimental output of NEXUS file for pre-aligned files!
# - 26-02-2009 Removal of 3 items from the usage stats file as they were of no use.
# * 26-02-2009 HTML Output has been altered to be easier to read and to reflect main program design.
# * 20-10-2008 Fixed Oryza sativa (japonica cultivar) problem, now parses [], () and then anything whilst ignoring [()]
# * 20-10-2008 Updated ex.ac.uk to exeter.ac.uk
# * 20-10-2008 Removed dependency on cgi-lib.pl and changed to CGI.pm which is installed by default
# + 03-06-2008 Head / Tail selection for Accession
# * 03-06-2008 Usage.txt is now tab delimited and includes document_root information.
# * 29-04-2008 Fixed logo image resize and tidied html
# * 20-02-2008 Increased duplicate addition of an identifying letter from 26 to 18278
# * XX-XX-2007 Made changes to populate() to remove chars such as :;._| from species name... Better supportd with \W*
# * XX-XX-2007 Fixed bug where () and [] and text were not being properly treated.

# Import packages
use Cwd;
use File::Path;
use File::Basename;
use Switch;
use CGI qw(:cgi-lib);
use CGI::Carp qw(fatalsToBrowser set_message);
set_message(
"Please send any errors messages with a description of what you were trying to do and any relevant files to cs02gl\@exeter.ac.uk"
);

our $VERSION = "0.6.8.2";

# Max Upload (1024 * 51200 = 50MB)
$CGI::POST_MAX = 1024 * 51200;
my $upload_dir = "\/var\/www";
my $domain = "http://localhost";

# Main
&session_id;
&upload;
&convertEOL;
&file2array;
&populate;
&output_files;
&tidy_up;
&usage;
&process_form;
# End

# Assign unique session ID to user each run
sub session_id {

    # Get current time - used as "unique" identifier
    $id = time();
    mkdir( "$upload_dir\/results\/rg\/$id", 0755 );
    return $id;
}

sub upload {

    #my $query    = new CGI;
    $query = CGI->new;

    # Filehandle for uploaded file
    $fileName = $query->param("file");
    if ( !$fileName ) {

        die( $query->header()
              . "No file selected, please choose a file to upload.\n" );
    }

#    # Get the MIME type for the upload file
#    $type = $query->uploadInfo($fileName)->{'Content-Type'};

# All uploads should be of text/plain, however without a .txt extension
# files are not recognised as plain text and so are assigned application/octet-stream.
# Other files are application/octet-stream and so we need to only allow a subset of
# these based on other acceptable extensions e.g. *.fas and *.fasta
#    unless ( $type eq 'text/plain' or $type eq 'application/octet-stream' ) {
#        die "Only text/plain files are allowed\nYour type = $type\n";
#    }

    ( $file, $dir, $ext ) = fileparse( $fileName, qr{\..*} );

    # Check for tainting...
    # Convert any spaces to underscores "_"
    $file =~ tr/ /_/;
    $ext  =~ tr/ /_/;

    # Remove illegal characters
    $allowed_filename_characters = 'a-zA-Z0-9_.-';
    $file =~ s/[^$allowed_filename_characters]//g;
    $ext  =~ s/[^$allowed_filename_characters]//g;
    if ( $file =~ /^([$allowed_filename_characters]+)$/ ) {
        $file = $1;
    }
    else {
        die(
"The filename is not valid. Filenames can only contain these characters: $allowed_filename_characters\n Please explain this to the rubber duck...\n"
        );
    }

    # This isn't perfect as you could rename any file .fas or .fasta
    # and it would become application/octet-stream but it's the best
    # method available for now.
    if ( $ext eq ".fas" or $ext eq ".fasta" or $ext eq ".txt" ) {

        # Open output file
        open( $input_file, '>', "$upload_dir\/results\/rg\/$id\/$fileName" )
          or die("Error 1: Unable to open file $fileName: $!\nPlease explain this to the rubber duck...\n");
        binmode $input_file;
        $upload_filehandle = $query->upload("file");

        # Print the contents of the uploaded file
        while (<$upload_filehandle>) {
            print $input_file $_;
        }
        close($input_file);
    }
    else {
        die("File type(extension) must be FASTA, please try again.\nPlease explain this to the rubber duck...\n");
    }
    
    return "Upload Succesful!\n";
}

# Conversion of EOL from input to native form.
sub convertEOL {
    open( $in_file, '<', "$upload_dir\/results\/rg\/$id\/$fileName" )
      or die("Error 2: Unable to read file $fileName: $!\n");
    open( $out_file, '>', "$upload_dir\/results\/rg\/$id\/$fileName\_EOLN" )
      or die("Error 3: Unable to save file $fileName\_EOLN: $!\n");

    #Convert EOL from \r to \n
    while (<$in_file>) {
        my ($line) = $_;

        $line =~ s/>/\015>/gs;
        $line =~ s/(\015?\012|\015)/\n/gs;

        print $out_file "$line";
    }
    close($in_file);
    close($out_file);
    return "EOL Conversion Succesful\n";
}

# Read the data file into an array, allows for more efficient handling
sub file2array {

    # Read in the deflines
    open( $in_defline, '<', "$upload_dir\/results\/rg\/$id\/$fileName\_EOLN" )
      or die("Error 4: Unable to read file $fileName\_EOLN: $!\n");

    while (<$in_defline>) {
        my ($line) = $_;
        if ( $line =~ m/>/ ) {
            chomp($line);
            push( @deflineArray, $line );
        }
    }
    close($in_defline);

    # Read in the sequences, line by line concatening the strings to one
    open( $in_sequence, '<', "$upload_dir\/results\/rg\/$id\/$fileName\_EOLN" )
      or die("Error 5: Unable to read file $fileName\_EOLN: $!\n");
    while (<$in_sequence>) {
        my ($line) = $_;
        if ( $line !~ m/>/ ) {
            chomp($line);
            $temp = "$temp$line";
        }
        else {
            push( @dnaaaArray, $temp );
            $temp = "";
        }
    }
    close($in_sequence);

    # Print last sequence not printed in loop and remove all white space
    push( @dnaaaArray, $temp );
    @dnaaaArray = grep { /\S/ } @dnaaaArray;

    $deflineSize = @deflineArray;
    $dnaaaSize   = @dnaaaArray;

    unlink("$upload_dir\/results\/rg\/$id\/$fileName\_EOLN");
    return;
}

sub populate {

    # Get accession length from cgi web form
    $accession_length = $query->param("acclen");

    # Get direction to chop string?
    $acc_head_tail = $query->param("head_tail");

    # Get Species Name length from cgi web form
    @sp_length[0] = $query->param("species");
    @sp_length[1] = $query->param("name");

    # Get separator from cgi webform
    $separator = $query->param("sep");

    &process_deflines;
    &remove_duplicates;
    &modify_accessions;
    &modify_taxa_names;
    &check_for_duplicate_IDs;
    return;
}

# Updated November 2009
# More explict matching on deflines with expansion to include new types
# any that do not fit exact matches are dumped in the non-match array
sub process_deflines {

    for ( $i = 0 ; $i < $deflineSize ; $i++ ) {

        # Read array line by line
        $line = $deflineArray[$i];

        switch ($line) {
            ## Long gi defintion - Standard NCBI Output from BLASTp search
    # >gi|71649003|ref|XP_813263.1| hypothetical protein [Trypanosoma cruzi]
    # Accession is 3rd element
    # Species name is generally in square braces [] of 4th element for proteins.
            case m/(>gi)(\|)(.*?)(\|)(.*?)(\|)(.*?)(\|)(.*?)/ {

                # Split the line by the pipe symbol, 5 in length
                @current_defline = split( /\|/, $line );

                # Accession Number to Array
                push( @accessionArray, "$current_defline[3]" );

                # Species Name to Array
                push( @species_nameArray, "gil$current_defline[4]" );
            }

            ## Short gi defintion - conceptual translations?
   # >gi|451623       (U04987) env [Simian immunodeficiency virus]
   # Accession is usually within round braces () 3rd else use gi after pipe| 2nd
   # Species name is not present, default to letters of words 4th
            case m/(>gi)(\|)(.*?)(\(.*\))(.*?)/ {

                # Split the line by the pipe symbol, 5 in length
                @current_defline = split( /\|/, $line );

                if ( $current_defline[3] ne "" ) {

                    # Accession Number to Array
                    push( @accessionArray, "$current_defline[3]" );
                }
                else {

                    # GI Number to Array
                    push( @accessionArray, "$current_defline[2]" );

                }

                # Species Name to Array
                push( @species_nameArray, "gis$current_defline[4]" );
            }

            ## JGI non-scaffold
      # >jgi|Phyra1_1|71587
      # Accession is 2nd element after 2nd pipe|
      # Species name is already concatenated and is 1st element, after 1st pipe|
            case m/(>jgi)(\|)(.*?)(\|)(.*?)/ {

                # Split the line by the pipe symbol, 5 in length
                @current_defline = split( /\|/, $line );

                # Accession Number to Array
                push( @accessionArray, "$current_defline[2]" );

                # Species Name to Array
		# Automatically change to Genus species from list
		&get_jgi_species_name("$current_defline[1]");
                push( @species_nameArray, "jgi$jgi_proper" );
            }

            ## JGI scaffold
      # >fgenesh1_pg.C_scaffold_1000025 [Phyra1_1:72583]
      # Accession is element after colon
      # Species name is element before colon
            case m/(>fgene)(.*?)(\\[.*?\\])/ {

                # Split the line by colon
                @current_defline = split( /\:/, $line );

                # Accession Number to Array
                push( @accessionArray, "$current_defline[1]" );

                # Species Name to Array
		# Automatically change to Genus species from list
		&get_jgi_species_name("$current_defline[0]");
                push( @species_nameArray, "jgis$jgi_proper" );
            }

            ## General database identifier
            # This covers gb, emb, dbk, sp and ref or any similar formats
            # gb|accession|locus
            # Accession is 1st element
            # Species name is potentially absent, use 2nd element
            case m/(>.{2,3})(\|)(.*?)(\|)(.*?)/ {

                # Split the line by the pipe symbol, 3 in length
                @current_defline = split( /\|/, $line );

                # Accession Number to Array
                push( @accessionArray, "$current_defline[1]" );

                # Species Name to Array
                push( @species_nameArray, "gb$current_defline[2]" );

            }

            else {
                push( @non_gb_jgi, $line );
                push( @non,        $dnaaaArray[$i] );
                delete $deflineArray[$i];
                delete $dnaaaArray[$i];

            }
        }
    }
    return;
}

sub get_jgi_species_name {
	
	my $jgi_accession = shift;
	open( $in_jgi, '<', "$upload_dir\/results\/doe_jgi_list.csv" ) or die("Error: Unable to read DOE JGI file $in_jgi: $!\n");
	while (<$in_jgi>) {
		my ($line) = $_;
		push( @jgi_name_array, $line );
	}
	close($in_jgi);
	$jgi_name_array_length = @jgi_name_array;

	for (my $i = 0; $i < $jgi_name_array_length; $i++) {
		
		my @temp = split (/,/, $jgi_name_array[$i]);
		$jgi_acc = $temp[0];
		$jgi_locus = $temp[1];
		chomp($jgi_locus);
		if ($jgi_accession eq $jgi_acc) {

			$jgi_proper = $jgi_locus;
			return "$jgi_proper";
			
		}
		#else {
		#	$jgi_proper = $jgi_accession;
		#	return "$jgi_proper";
		#}
	}
	$jgi_proper = $jgi_accession;
	return "$jgi_proper";
}

sub modify_accessions {

    # Modify the accession - make it shorter...
    for ( $j = 0 ; $j < $accessionArray_length ; $j++ ) {

        # Assign value from accessionArray to short_accArray
        $short_accArray[$j] = $accessionArray[$j];

        # This should remove all non-alphanumeric chars, better than above!
        $short_accArray[$j] =~ tr/\W*\.*//d;

        # Create sub string of value...
        if ( $acc_head_tail eq "head" ) {
            $short_accArray[$j] =
              substr( $short_accArray[$j], 0, $accession_length );
        }
        else {
            $short_accArray[$j] =
              substr( $short_accArray[$j], -$accession_length );
        }
    }
    return;
}

sub modify_taxa_names {

    # Modify Species Name!!!
    @spName = @species_nameArray;
    for ( $k = 0 ; $k < $spnameArray_length ; $k++ ) {
        $type = $spName[$k];

        # Remove any non-alphanumerics
        $type =~ tr/\W*\.*\:*//d;

        #$type =~ s/^\s+//;
        $spName[$k] = $type;
        switch ($type) {

            # DOE JGI Format
            case m/^jgi/i {
                $length = length $type;

                # Remove jgi tag...
                $spName[$k]       = substr( $spName[$k], 3, $length );
		@species = split( / /, $spName[$k] );
                $sp   = substr( $species[0], 0, $sp_length[0] );
                $name = substr( $species[1], 0, $sp_length[1] );
                $short_spName[$k] = "$sp$name";
                #$short_spName[$k] = substr( $spName[$k], 0, $sp_length[0] );
            }

         # GenBank / EMBL Data Lib / DDBJ / Swiss-Prot / NCBI Reference Sequence
            case m/^gb/i {
                $length = length $type;

                # Remove gb tag...
                $spName[$k] = substr( $spName[$k], 2, $length );

                @species = split( / /, $spName[$k] );
                $sp   = substr( $species[0], 0, $sp_length[0] );
                $name = substr( $species[1], 0, $sp_length[1] );
                $short_spName[$k] = "$sp$name";
            }

            # Short gi from NCBI
            case m/^gis/i {
                $length = length $type;

                # Remove gis tag...
                $spName[$k] = substr( $spName[$k], 3, $length );

                @species = split( / /, $spName[$k] );
                $sp   = substr( $species[0], 0, $sp_length[0] );
                $name = substr( $species[1], 0, $sp_length[1] );
                $short_spName[$k] = "$sp$name";
            }

            # gi identifiers from NCBI
            case m/^gil/i {

# There is an assumed preference for text in [] followed by () then whatever is left
# if $type contains square braces []
                if ( $type =~ m/(.*?)(\[.*?\])/is ) {
                    $species = $2;

# if within the [] there are round braces (), ignore e.g. [Oryz sativa (japonica cultivar)]
                    if ( $species =~ m/(.*?)(\(.*\))/is ) {
                        $species = $1;

                        # Trim all braces and spaces
                        $species =~ tr/\(\)\[\]/ /d;
                        $species =~ s/^\s+//;
                        $species =~ s/\s+$//;

                        $spName[$k] = $species;
                        @species = split( / /, $species );
                        $sp   = substr( $species[0], 0, $sp_length[0] );
                        $name = substr( $species[1], 0, $sp_length[1] );
                        $short_spName[$k] = "$sp$name";
                    }
                    else {
                        $species =~ tr/\[\]/ /d;
                        $species =~ s/^\s+//;
                        $species =~ s/\s+$//;

                        $spName[$k] = $species;
                        @species = split( / /, $species );
                        $sp   = substr( $species[0], 0, $sp_length[0] );
                        $name = substr( $species[1], 0, $sp_length[1] );
                        $short_spName[$k] = "$sp$name";
                    }
                }

                # Else if it contains round braces
                elsif ( $type =~ m/(\(.*\))/is ) {
                    $species = $1;

                    # Remove braces and spaces
                    $species =~ tr/\(\)/ /d;
                    $species =~ s/^\s+//;
                    $species =~ s/\s+$//;

                    $spName[$k] = $species;
                    @species = split( / /, $species );
                    $sp   = substr( $species[0], 0, $sp_length[0] );
                    $name = substr( $species[1], 0, $sp_length[1] );
                    $short_spName[$k] = "$sp$name";
                }
                else {
                    $length = length $type;

                    # Remove gil tag...
                    $spName[$k] = substr( $spName[$k], 3, $length );

                    @species = split( / /, $spName[$k] );
                    $sp   = substr( $species[0], 0, $sp_length[0] );
                    $name = substr( $species[1], 0, $sp_length[1] );
                    $short_spName[$k] = "$sp$name";
                }
            }
        }
    }
    return;
}

sub check_for_duplicate_IDs {

    # Create + populate refgenID array
    @refgenIDArray = ();
    for ( $c = 0 ; $c < $accessionArray_length ; $c++ ) {
        $refgenIDArray[$c] = "$short_accArray[$c]$separator$short_spName[$c]";
    }
    $refgenIDArray_length = @refgenIDArray;

    # Gives 18278 possible changes which should be more
    # than enough for any large phylogenies
    @count   = ( "A" .. "ZZZ" );
    $counter = 0;

    @dupeMatchRef = ();
    for ( $d = 0 ; $d < $refgenIDArray_length ; $d++ ) {
        $testID = $refgenIDArray[$d];
        for ( $f = $d + 1 ; $f < $refgenIDArray_length ; $f++ ) {
            $matchID = $refgenIDArray[$f];
            if ( $testID eq $matchID ) {

                $matchID_length    = length $matchID;
                $matchID           = substr( $matchID, 0, $matchID_length - 1 );
                $matchID           = "$matchID$count[$counter]";
                $refgenIDArray[$f] = $matchID;

                push( @dupeMatchRef, "$matchID\t\t\t$accessionArray[$d]" );
                $counter++;
            }
        }
    }
    $dupeMatchRef_length = @dupeMatchRef;
    return;
}

sub remove_duplicates {
    $accessionArray_length = @accessionArray;
    undef(@duplicates);
    for ( $i = 0 ; $i < $accessionArray_length ; $i++ ) {
        $test_line = $accessionArray[$i];
        for ( $j = $i + 1 ; $j < $accessionArray_length ; $j++ ) {
            $match_line = $accessionArray[$j];
            if ( $test_line eq $match_line ) {
                $duplicates[$i] = $accessionArray[$j];
                delete $accessionArray[$j];
                delete $species_nameArray[$j];
                delete $dnaaaArray[$j];
            }
        }
    }

    # Remove white space
    @dnaaaArray        = grep { /\S/ } @dnaaaArray;
    @accessionArray    = grep { /\S/ } @accessionArray;
    @species_nameArray = grep { /\S/ } @species_nameArray;
    @duplicates        = grep { /\S/ } @duplicates;

    # Reset lengths
    $accessionArray_length = @accessionArray;
    $spnameArray_length    = @species_nameArray;
    $dnaaaSize             = @dnaaaArray;
    $duplicates_length     = @duplicates;
    return;
}

sub output_files {
    $outName     = "$file\_refgen$ext";
    $outNEXUS    = "$file\_refgen.nex";
    $speciesName = "$file\_species.csv";
    $tree_graphName = "$file\_tree_graph.txt";
    &refgen_out;
    &nexus_out;
    &species_out;
    &tree_graph_name_table;
    &messages_out;
    $fileName = $outName;
    &changeEOL;
    $fileName = $outNEXUS;
    &changeEOL;
    return;
}

sub refgen_out {

    # REFGEN File
    open( $output, '>', "$upload_dir\/results\/rg\/$id\/$outName\_preEOL" )
      or die("Error: Unable to save file $outName\_preEOL: $!\n");
    for ( $m = 0 ; $m < $refgenIDArray_length ; $m++ ) {
        print $output ">$refgenIDArray[$m]\n";
        @pieces = $dnaaaArray[$m] =~ m/(.{1,70})/gs;
        foreach (@pieces) {
            print $output "$_\n";
        }
        print $output "\n";
    }
    print $output "\n\n";

    @non = grep { /\S/ } @non;
    $non_length = @non;
    for ( $n = 0 ; $n < $non_length ; $n++ ) {
        print $output "$non_gb_jgi[$n]\n";
        @pieces = $non[$n] =~ m/(.{1,70})/gs;
        foreach (@pieces) {
            print $output $_ . "\n";
        }
        print $output "\n";
    }
    close($output);
    return;
}

sub nexus_out {
    ### NEXUS File
    $maxLength = 0;
    $datatype  = "DNA";
    ## Get the largest sequence length for NEXUS nchar value and check for PROTEIN or DNA
    ## Idealy this should only be done if all sequence lengths are the same - i.e. when I can be bothered to
    ## add a checkbox or some form of checking system for now it just happens.
    for ( $nn = 0 ; $nn < $non_length ; $nn++ ) {
        $seq_length = length( $non[$nn] );
        if ( $seq_length ge $maxLength ) {
            $maxLength = $seq_length;
        }
    }
    for ( $mm = 0 ; $mm < $refgenIDArray_length ; $mm++ ) {
        $code = $dnaaaArray[$mm];
        if ( $code =~ m/[^BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvxxyz]/i ) {
            $datatype = "PROTEIN";
        }
        $seq_length = length( $refgenIDArray_length[$mm] );
        if ( $seq_length ge $maxLength ) {
            $maxLength = $seq_length;
        }

    }

    open( $nexus_output, '>',
        "$upload_dir\/results\/rg\/$id\/$outNEXUS\_preEOL" )
      or die("Error: Unable to open file $outNEXUS\_preEOL: $!\n");

    $ntax = $non_length + $refgenIDArray_length;

    print $nexus_output "#NEXUS\n";
    print $nexus_output "[File created by REFGEN Version $VERSION on "
      . localtime() . "]\n";
    print $nexus_output "BEGIN DATA;\n";
    print $nexus_output "  DIMENSIONS NTAX=$ntax NCHAR=$maxLength;\n";
    print $nexus_output "  FORMAT DATATYPE=$datatype\n  ;\n";
    print $nexus_output "MATRIX\n";
    $count = 1;
    for ( $mmm = 0 ; $mmm < $refgenIDArray_length ; $mmm++ ) {
        print $nexus_output "\[$count\] $refgenIDArray[$mmm]\n";
        @pieces = $dnaaaArray[$mmm] =~ m/(.{1,70})/gs;
        foreach (@pieces) {
            print $nexus_output "$_\n";
        }
        print $nexus_output "\n";
        $count++;
    }
    print $nexus_output "\n";
    for ( $nnn = 0 ; $nnn < $non_length ; $nnn++ ) {
        print $nexus_output "\[$count\] $non_gb_jgi[$nnn]\n";
        @pieces = $non[$nnn] =~ m/(.{1,70})/gs;
        foreach (@pieces) {
            print $nexus_output $_ . "\n";
        }
        print $nexus_output "\n";
        $count++;
    }
    print $nexus_output ";\nEND;\n";
    close($nexus_output);
    ###
    return;
}

sub tree_graph_name_table {

    # Tree Graph name table...
    open( $treegraph_out, '>', "$upload_dir\/results\/rg\/$id\/$tree_graphName" )
      or die("Error 6a: Unable to save file $tree_graphName: $!\n");
    for ( $spec = 0 ; $spec < $accessionArray_length ; $spec++ ) {
        print $treegraph_out
          "$refgenIDArray[$spec]\t$spName[$spec]\n";
    }
    close($treegraph_out);
    return;
}

sub species_out {

    #Species List
    open( $species_out, '>', "$upload_dir\/results\/rg\/$id\/$speciesName" )
      or die("Error 6: Unable to save file $speciesName: $!\n");
    for ( $spec = 0 ; $spec < $accessionArray_length ; $spec++ ) {
        print $species_out
          "$accessionArray[$spec],$refgenIDArray[$spec],$spName[$spec],\n";
    }
    close($species_out);
    return;
}

sub messages_out {

    # Duplicates list...
    open( $messages_out, '>', "$upload_dir\/results\/rg\/$id\/messages.txt" )
      or die("Error 7: Unable to save file $fileName: $!\n");

    print $messages_out "REMOVED SEQUENCES \($duplicates_length\)\n";
    foreach my $dupe (@duplicates) {
        if ( $dupe ne "" ) {
            print $messages_out "$dupe\n";
        }
    }
    print $messages_out "\nREFGEN ID DUPLICATES \($dupeMatchRef_length\)\n";
    print $messages_out "New REFID\t\tOriginal Accession\n";
    foreach my $match (@dupeMatchRef) {
        print $messages_out "$match\n";
    }
    print $messages_out "\nNON GenBank or DOE JGI FORMATTED SEQUENCES\n";
    foreach (@non_gb_jgi) {
        print $messages_out $_ . "\n";
    }
    print $messages_out
"\nYou will need to alter these sequences yourself to give them a\nREFGEN ID. All sequences and FASTA deflines will be appeneded to\nthe end of the REFGEN file. Don't forget to include the code and \nrelevant information in the .csv file (using a text editor) in the\nformat 'Accession,ID,Species Name' where ID is the ID you created in\nthe REFGEN file. Don't forget all FASTA deflines start with '>'";

    close($messages_out);

    # Test if non array has values for message in final output...
    if (@non) {
        $non_message =
"<font color=\"#8E0A00\">Warning:</font> Your FASTA file contains sequences that cannot be processed due to not being in a standard NCBI or DOE JGI format.<br />Please check the messages file to resolve this issue.";
    }
    else {
        $non_message = "";
    }
    return;
}

# Reverse EOL line changes
sub changeEOL {

    # Open the saved formatted file
    open( $pre_eol, '<', "$upload_dir\/results\/rg\/$id\/$fileName\_preEOL" )
      or die("Error 8: Unable to open file $fileName: $!\n");

    # Open file for writing
    open( $post_eol, '>', "$upload_dir\/results\/rg\/$id\/$fileName" )
      or die("Error 9: Unable to open file $fileName: $!\n");

    $os = $query->param("opsys");
    switch ($os) {

        case "win" {
            while (<$pre_eol>) {
                my ($line) = $_;

                #convert to win
                $line =~ s/\n/\015\012/gs;
                print $post_eol "$line";
            }
        }
        case "nix" {

            # Nothing to do but dump contents to new file
            while (<$pre_eol>) {
                my ($line) = $_;
                print $post_eol "$line";
            }
        }
        case "mac" {
            while (<$pre_eol>) {
                my ($line) = $_;

                # Convert to mac EOL
                $line =~ s/\n/\015/gs;
                print $post_eol "$line";
            }
        }
        else {

            # Leave it as it is...
            # Nothing to do but dump contents to new file
            while (<$pre_eol>) {
                my ($line) = $_;
                print $post_eol "$line";
            }
        }
    }

    #Tidy up
    close($pre_eol);
    close($post_eol);

    # Remove preEOL temporary file
    unlink("$upload_dir\/results\/rg\/$id\/$fileName\_preEOL");
    return;
}

sub usage {

    $user_IP      = $ENV{REMOTE_ADDR};
    $user_browser = $ENV{HTTP_USER_AGENT};
    $time         = localtime time;
    open( $usage_out, '>>', "$upload_dir\/results\/usage\/refgen.txt" )
      or die("Error 10: Unable to write to file $upload_dir\/results\/usage\/refgen.txt: $!\n");
    print $usage_out "$time\t$user_IP\t$id\t$user_browser\n";
    return;
}

sub html_out {

    open( $html_out, '>', "$upload_dir\/results\/rg\/$id\/$id.html" )
      or die("Error 11: Unable to save file $fileName: $!\n");

    print $html_out <<"ENDOFTEXT";
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">

<html xmlns="http://www.w3.org/1999/xhtml">
<head>
  <title>REFGEN Results - Centre for Eukaryotic Evolutionary
  Microbiology</title>
  <link href="$domain/css/ceem.css" rel=
  "stylesheet" type="text/css" />
  <link href="$domain/css/form.css" rel=
  "stylesheet" type="text/css" />
</head>

<body>
  <div align="center" class="container">
    <form>
      <table align="center" id="table1">
        <tbody>
          <tr>
            <td><span class="pi_name"><img src=
            "$domain/images/LogoTOOLS.png" /></span></td>
          </tr>

          <tr>
            <td><span class="pi_sub_title">REFGEN Results for
            $file$ext</span></td>
          </tr>

          <tr>
            <td class="help_text">Please right click and
            choose "Save as..." to save your files.</td>
          </tr>

          <tr>
            <td>
              <fieldset>
                <legend>ID Code</legend>

                <p class="help_text">Your unique identifier code is
                <font color="#449E00" size="4">$id</font>.</p>

                <p class="help_text">Please bookmark <a href=
                "$domain/results/rg/$id/$id.html">
                <font color="#449E00" size="3">this</font></a> URL
                if you wish to review your results at another
                time.<br />
                Your results will be stored and will be available
                for around <em><font color="#8E0A00">one
                week</font></em>.</p>
              </fieldset>
            </td>
          </tr>

          <tr>
            <td>
              <fieldset>
                <legend>Original Files</legend>

                <p class="help_text"><a href=
                "$domain/results/rg/$id/$file$ext">
                Your Original Sequence File</a></p>
              </fieldset>
            </td>
          </tr>
          <tr>
            <td>
              <fieldset>
                <legend>Messages</legend>
                <p class="help_text">&gt;&gt; $duplicates_length
                duplicate sequences removed.<br />
                &gt;&gt; $dupeMatchRef_length REFGEN IDs
                altered.</p>

                <p class="help_text">$non_message</p>

                <p class="help_text"><a href=
                "$domain/results/rg/$id/messages.txt">
                Messages File</a> - A detailed version of the
                above.</p>
              </fieldset>
            </td>
          </tr>
          <tr>
            <td>
              <table width="100%">
                <tr>
                  <td width="50%">
                    <fieldset>
                      <legend>Output Files</legend>

                      <p class="help_text"><a href=
                      "$domain/results/rg/$id/$file\_refgen$ext">
                      REFGEN Formatted File</a></p>

                      <p class="help_text"><a href=
                      "$domain/results/rg/$id/$file\_species.csv">
                      Species List</a> for use with <a href=
                      "$domain/treenamer.html">
                      <font color=
                      "#449E00">TREENAMER</font></a></p>
                    </fieldset>
                  </td>

                  <td width="50%">
                    <fieldset>
                      <legend>Experimental</legend>

                      <p class="help_text"><a href=
                      "$domain/results/rg/$id/$file\_refgen.nex">
                      REFGEN NEXUS File</a></p>

                      <p class="help_text">
                      <a href=
                      "$domain/results/rg/$id/$file\_tree_graph.txt">
                      Tree Graph 2 - Name Table</a>
                      </p>
                    </fieldset>
                  </td>
                </tr>
              </table>
            </td>
          </tr>

          <tr>
            <td>
              <fieldset>
                <legend>Information</legend>

                <p class="help_text">Thank you for using REFGEN
                $VERSION.</p>

                <p class="help_text">Other Tools: <a href=
                "$domain/treenamer.html"><font color="#449E00">
                TREENAMER</font></a> - A tool to change the unique
                identifiers made by REFGEN to a more 'human'
                readable format (accession &amp; species name) in
                your tree files.</p>

                <p class="help_text">Please come back again and
                don't forget to suggest these tools to your
                friends!</p>
                
                                <p class="help_text">Citation: <a href="http://www.la-press.com/refgen-and-treenamer-automated-sequence-data-handling-for-phylogenetic-a1451">Leonard, G., et al. (2009). REFGEN and TREENAMER: Automated Sequence Data Handling for Phylogenetic Analysis in the Genomic Era. <em>Evolutionary Bioinformatics Online</em>, 5:1-4.</a> 
								<br />PDF: <a href="pdf/Leonard_(2009).pdf">Click here to download...</a> 
								</p> 
              </fieldset>
            </td>
          </tr>

    <tr> 
      <td class="pi_menu"><a href="$domain/index.html">CEEM Home</a> &#8226; <a href=
      "$domain/ceem_staff_tom.html">Tom Home</a> &#8226; <a href="$domain/ceem_staff_tom_people.html">People</a> 
      &#8226; <a href="$domain/ceem_staff_tom_research.html">Research</a> &#8226; <a href=
      "$domain/ceem_staff_tom_tools.html">Tools</a> &#8226; <a href=
      "$domain/ceem_staff_tom_publications.html">Publications</a> &#8226; <a href=
      "javascript:history.go(-1)">Go Back</a></td> 
    </tr> 

          <tr>
            <td>
              <div align="center">
                <span class="copyright">Copyright CEEM MMVII</span>
              </div>
            </td>
          </tr>
        </tbody>
      </table>
    </form>
  </div>
</body>
</html>
ENDOFTEXT
    return;
}

sub redirect_head {
    return <<"END_OF_TEXT";
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
  <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
  <title>Centre for Eukaryotic Evolutionary Microbiology</title>
<META HTTP-EQUIV="Refresh" CONTENT="2; URL=$domain/results/rg/$id/$id.html">
  <link href="$domain/css/ceem.css" rel="stylesheet" type="text/css" />
</head>
END_OF_TEXT

}

sub process_form {
    print &PrintHeader;
    print &redirect_head;

    print <<"ENDOFTEXT";
<body>
<table align="center" border="0" cellpadding="0" cellspacing="0" width="600">
  <tbody>
    <tr>
      <td><span class="logo_title"><img src="$domain/images/LogoTOOLS.png" /></span></td>
    </tr>
  <tr>
    <td><span class="pi_sub_title"><strong>Generating Results</strong></span></td>
  </tr>
  <tr>
    <td height="450" class="text"><p class="pi_sub_title"><br />REFGEN Results for $file$ext</span></p>
		<p class="profile_text">Your unique identifier code is $id...</p>
        <p class="profile_text"><img src="$domain/images/spinning_wheel.gif" /><img src="$domain/images/working.png" /><img src="$domain/images/spinning_wheel.gif" /></p><p class="profile_text">$VERSION</p></td>
  </tr>
  <tr>
    <td>
      <div class="pi_menu"><a href="$domain/index.html" class="pi_menu">CEEM Home</a>
&bull; <a href="$domain/CEEM_STAFF_PAGESTAR.html" class="pi_menu">Tom
Home</a> &bull; <a href="$domain/CEEM_STAFF_PAGESTARpeople.html" class="pi_menu">People</a><a> &bull; </a><a href="$domain/CEEM_STAFF_PAGESTARres.html" class="pi_menu">
Research</a> &bull; <a href="$domain/CEEM_STAFF_PAGESTARtools.html" class="pi_menu">
Tools </a> &bull; <a href="$domain/CEEM_STAFF_PAGESTARpubs.html" class="pi_menu"> Publications</a> &bull; <a href="javascript:history.go(-1)" class="pi_menu">Go
Back</a></div>
      </td>
  </tr>
  <tr>
    <td><div align="center"><span class="copyright">Copyright CEEM MMVII - MMX</span></div></td>
  </tr>
</table>
ENDOFTEXT
    print &HtmlBot;
    &html_out;
    return;
}

sub tidy_up {

    # Clear up old file results
    $time_now = time();
    $d        = "$upload_dir\/results\/rg";
    opendir( DIR, $d ) or die("Cannot open directory $dir: $!\n");
    @dirs = grep { !/^\./ && -d "$d/$_" } readdir(DIR);
    closedir DIR;
    $dir_size = @dirs;
    for ( $i = 0 ; $i < $dir_size ; $i++ ) {
        $diff = $time_now - $dirs[$i];

        # Delete after a week (in seconds)
        if ( $diff >= 604800 ) {
            rmtree( "$d\/$dirs[$i]", 0, 0 );
        }
    }
    return;
}
