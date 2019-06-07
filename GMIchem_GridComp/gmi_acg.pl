#!/usr/bin/perl
# --------------------------------------------------------------------------
# ACG for GMICHEM, Usage notes below.
#
# Revision history:
# Tue 12 Apr 2011   Nielsen    Origination.
# Tue 27 May 2011   Nielsen    Added reaction rate exports.
# Tue  2 Jun 2011   Nielsen    Added chemical rate of change exports.
# Wed 17 Aug 2011   Nielsen    Modifications for Joules' grid comp classes.
# Fri 25 May 2012   Nielsen    Export specs for AGCM_ and PHYS_GridComp
# Thu 12 Sep 2013   Nielsen    Now using setkin_lchem.h instead of setkin_chem_mech.txt
# --------------------------------------------------------------------------
use File::Basename;
use Getopt::Std;         # command line options

getopts('d:Hm:p:Rr:v');
usage() if ( $opt_H );

$Iam = "gmi_acg";
$i = -9999;
$m = -1;
@process = ( "SCAV_" , "DD_" , "WD_" );
@getPointerFile = ( "Deposition" , "Reactions" );

# Echo our progress
# -----------------
print "\n";
if ( $opt_R ) {
 print "$Iam".": Creating additional registry files ...\n";
} else {
 print "$Iam".": Prepare files for filling exports ... \n";
}
print "\n";

# Satisfy required file names
# ---------------------------
$registryFile   = "GMI_GridComp/GMI_Registry.rc"                           unless ( $registryFile = $opt_r );
$mechFileName   = "GMI_GridComp/GmiChemistry/stratTropMech/setkin_lchem.h" unless ( $mechFileName = $opt_m );
$deposParamFile = "GMI_GridComp/GmiChemistry/stratTropMech/setkin_depos.h" unless ( $deposParamFile = $opt_d );
$mechParamFile  = "GMI_GridComp/GmiChemistry/stratTropMech/setkin_par.h"   unless ( $mechParamFile = $opt_p );

# Find how many species and how many thermal 
# and photolysis reactions are in the mechanism.
# ----------------------------------------------
open FILE, $mechParamFile or die "Cannot find $mechParamFile\n";
while ( <FILE> ) {
 chomp($_);
 @tokens = split();
 chop($tokens[3]);
 if ( $tokens[1] eq "(NSP"   && $tokens[2] eq "=" ) { $numSpecies = $tokens[3]; }
 if ( $tokens[1] eq "(NUM_K" && $tokens[2] eq "=" ) { $num_k      = $tokens[3]; }
 if ( $tokens[1] eq "(NUM_J" && $tokens[2] eq "=" ) { $num_j      = $tokens[3]; }
 $i++;
}
close(FILE);
if ( $opt_v ) {
 print "Query of $mechParamFile shows the GMIchem mechanism has\n";
 print "$numSpecies species \n";
 print "$num_k thermal reactions \n";
 print "$num_j photolytic reactions \n";
 print "\n";
}

# Obtain the species names from the chemical kinetics text file
# -------------------------------------------------------------
open FILE, $mechFileName or die "Cannot find $mechFileName\n";
print "Downloading species names from $mechFileName\n";

# Read file line-by-line
# ----------------------
while ( <FILE> ) {

# Remove newline and grab tokens
# ------------------------------
 chomp($_);
 @tokens = split();

# List of species names begins after "Species indicies:"
# ------------------------------------------------------
 if ( $tokens[1] eq "All" && $tokens[2] eq "species" ) { $i = -1; }

# The next numSpecies lines contain the species names
# ---------------------------------------------------
 if ( $i > 0 && $i <= $numSpecies ) {

# Save the species formula in uppercase
# -------------------------------------
  $m++;
  $l = length( $tokens[2] ) - 4;
  $species[$m] = substr( $tokens[2], 2, $l );
  $species[$m] = uc $species[$m];

# Make species names conformal with GEOS-5
# ----------------------------------------
  if ( $species[$m] eq "O3"     ) { $species[$m] = "OX" ;      }
  if ( $species[$m] eq "CFCL3"  ) { $species[$m] = "CFC11" ;   }
  if ( $species[$m] eq "CF2CL2" ) { $species[$m] = "CFC12" ;   }

 }
 $i++;
}

close(FILE);

# Find the maximum length of the species names
# --------------------------------------------
$maxChars = 0;
for $i ( 0..$#species ) {
 $l = length( $species[$i] );
 if ( $l > $maxChars ) { $maxChars = $l; }
}

# Create blank-fill string for each species.  This 
# will allow us to make the generated text look pretty.
# -----------------------------------------------------
for $i ( 0..$#species ) {
 $l = length( $species[$i] );
 $blankFill[$i] = "";
 if ( $l < $maxChars ) {
  for $n ( 1..$maxChars-$l ) { $blankFill[$i] = $blankFill[$i]." "; }
 }
}

# Print the list of species if in verbose mode
# --------------------------------------------
if ( $opt_v ) {
 $n = $numSpecies/10;
 for $row ( 1..$n + 1 ) {
  for $col ( 1..10 ) {
   $i = 10 * ( $row - 1) + $col - 1;
   if ( $i < $numSpecies ) { print "$blankFill[$i]$species[$i] "; }
  }
  print "\n";
 }
 print "\n";
}

# Grab the order of the species in the internal state
# ---------------------------------------------------
open FILE, $registryFile or die "Cannot find $registryFile\n";
$i = -9999;
$m = -1;

# Read file line-by-line
# ----------------------
while ( <FILE> ) {

# Remove newline and grab tokens
# ------------------------------
 chomp($_);
 @tokens = split();

# List of internal state begins after "<InternalSpec" + 5 lines
# -------------------------------------------------------------
 if ( $tokens[0] eq "<InternalSpec" ) { $i = 0; }
 if ( $i > 5 && $tokens[0] eq "#" ) { $i = -9999; }

 if ( $i > 4 ) {
  $m++;

# Species name is token 0
# -----------------------
  $internal_spec_name[$m] = $tokens[0];

# Count the number of tokens in this line
# ---------------------------------------
  $count = 0;
  foreach (@tokens) { $count++ };

# Which tokens are vertical lines "|"
# -----------------------------------
  $n = -1;
  $count = -1;
  foreach (@tokens) {
   $count ++;
   if ( $tokens[$count] eq "|" ) {
    $n ++;
    $vlpos[$n] = $count;
   }
  }

# Unit lies between 1st and 2nd vertical line, i.e., vlpos[0] and vlpos[1]
# ------------------------------------------------------------------------
  $internal_spec_unit[$m] = $tokens[2];
  $string = "";
  for $n ( $vlpos[0]+1..$vlpos[1]-1 ) { $string = "$string"." "."$tokens[$n]"; }
  $internal_spec_unit[$m] = $string;

# Friends are listed after the 11th vertical line, i.e., vlpos[10]+1
# ------------------------------------------------------------------
  $friends[$m] = $tokens[$vlpos[10]+1];
 }

 $i++;
 
}
$num_internal_state = $m + 1;
if ( $opt_v ) {
 print "Query of $registryFile shows the internal state has\n";
 print "$num_internal_state members\n\n";
}

close(FILE);

# Obtain the list of thermal reactions from the setkin_lchem include file
# -----------------------------------------------------------------------
open FILE, $mechFileName or die "Cannot find $mechFileName\n";
print "Downloading the thermal reaction list from $mechFileName\n";
$m = -1;
$sets = int( 1 + $num_k / 10 );
$s = -1;
$i = -2;

# Read file line-by-line
# ----------------------
while ( <FILE> ) {

# Remove newline and grab tokens
# ------------------------------
 chomp($_);
 @tokens = split();

# List of thermal reactions begins after "Thermal reaction labels"
# ----------------------------------------------------------------
 if ( $tokens[1] eq "Thermal" && $tokens[2] eq "reaction" ) { $s = 1; }

# There are $sets-1 groups of 10 reactions and 1 group with between 1 and 9 reactions
# ----------------------------------------------------------------------------------- 
 if ( $s > 0 && $s <= $sets) {
  if ( $s < $sets ) { $n = 10; } else { $n = $num_k % 10; }

# After two garbage lines, the next $n lines contain reactions
# ------------------------------------------------------------
  if ( $i > 0 && $i <= $n) {

# Strip leading and trailing characters and save reaction formula
# ---------------------------------------------------------------
   $m++;
   if ( $i == $n ) { $l = length( $_ ) - 11; } else { $l = length( $_ ) - 12; }
   $kReaction[$m] = substr ( $_ ,8 , $l);
   if ( $opt_v ) { print "  "."$kReaction[$m]\n"; }
  }

  if ( $i == $n ) { $i = -1; $s++; } else { $i++; }

 }
}

close(FILE);
if ( $opt_v ) { print " \n"; }

# Obtain the list of photolytic reactions from the setkin_lchem include file
# --------------------------------------------------------------------------
open FILE, $mechFileName or die "Cannot find $mechFileName\n";
print "Downloading the photolytic reaction list from $mechFileName\n";
$m = -1;
$sets = int( 1 + $num_j / 10 );
$s = -1;
$i = -2;

# Read file line-by-line
# ----------------------
while ( <FILE> ) {

# Remove newline and grab tokens
# ------------------------------
 chomp($_);
 @tokens = split();

# List of photolytic reactions begins after "Photolytic reaction labels"
# ----------------------------------------------------------------------
 if ( $tokens[1] eq "Photolytic" && $tokens[2] eq "reaction" ) { $s = 1; }

# There are $sets-1 groups of 10 reactions and 1 group with between 1 and 9 reactions
# ----------------------------------------------------------------------------------- 
 if ( $s > 0 && $s <= $sets) {
  if ( $s < $sets ) { $n = 10; } else { $n = $num_j % 10; }

# After two garbage lines, the next $n lines contain reactions
# ------------------------------------------------------------
  if ( $i > 0 && $i <= $n) {

# Strip leading and trailing characters and save reaction formula
# ---------------------------------------------------------------
   $m++;
   if ( $i == $n ) { $l = length( $_ ) - 11; } else { $l = length( $_ ) - 12; }
   $jReaction[$m] = substr ( $_ ,8 , $l);
   if ( $opt_v ) { print "  "."$jReaction[$m]\n"; }
  }

  if ( $i == $n ) { $i = -1; $s++; } else { $i++; }

 }
}

close(FILE);
if ( $opt_v ) { print " \n"; }

# 0. Derive Registries
# --------------------
    if ( $opt_R ) {
        register_exports();
        exit(0);
    }

# 1. Write the files for filling exports for scavenging,
#    dry deposition, and wet deposition, in that order.
# ------------------------------------------------------
for $i ( 0..2 ) {

    if ( $i == 0 ) { $hstar_string = "hstar_wet"; }
    if ( $i == 1 ) { $hstar_string = "hstar_dry"; }
    if ( $i == 2 ) { $hstar_string = "hstar_wet"; }

    print "Setting hstar to $hstar_string \n";

    find_depos();

    write_dep_h_files();
}

# 2. Write the files for filling exports for rates constants
#    and rates from the thermal and photolytic reactions.
# ----------------------------------------------------------
for $i ( 0..3 ) { write_q_h_files(); }

# 3. Change "EXPORT" to "expChem" in Get Pointers files
# -----------------------------------------------------
for $i ( 0..1 ) { revert(); }

print "\n";
exit(0);

sub revert {

# Mark GMI_GridComp when necessary
# --------------------------------
 $string = "GMI_GridComp/"."$getPointerFile[$i]";

# Open the two files
# ------------------
 $oldFile = "$string"."_GetPointer___.h";
 open OF, $oldFile or die "Cannot find $oldFile\n";

 $newFile = "$string"."_GetPointer2___.h";
 open NF, ">", $newFile or die "Cannot open $newFile\n";;
 print "Modifying $oldFile\n";

# Examine each line of the file
# -----------------------------
 while ($line = <OF>) {

# Revert state name to expChem or impChem
# ---------------------------------------
  $ec = index( $line, "EXPORT" );
  if ( $ec > 0 ) { $line =~ s/EXPORT/expChem/; }
  $ec = index( $line, "call M" );
  if ( $ec > 0 ) { $line =~ s/call M/CALL M/; }

  print NF $line;

 }

# Done. Close files.
# ------------------
close(OF);
close(NF);
}

sub write_dep_h_files {
 
# Choose file name and open
# -------------------------
 $fileName = "GMI_GridComp/$process[$i]"."FillExports___.h";
 open FILE, ">", $fileName or die "Cannot open $fileName\n";
 print "Building $fileName\n";

# Cycle through the species
# -------------------------
 for $m ( 0..$numSpecies-1 ) {
  $n = $m + 1;
  
# A spacer for tidiness if less than 1000 species
# -----------------------------------------------
  $prefix = "";
  if ( $n <  100 ) { $prefix = " "  ; }
  if ( $n <   10 ) { $prefix = "  " ; }

# Scan the list of desired exports and continue only if species(m) found in list
# ------------------------------------------------------------------------------
  $match = 0;
  for $d ( 0..$numDepos-1 ) {
   $deposName = $species[$deposSpeciesNumber[$d]];
   if ( $deposName eq $species[$m] ) { $match = 1; }
  }

# Concatenate species name to process abbreviation, right justified
# -----------------------------------------------------------------
  if ( $match == 1 ) {
   $string = "$blankFill[$m]"."$process[$i]"."$species[$m]";

# Complete line to be added to file
# ---------------------------------
   if ( $i == 0 ) {
    $line = "  IF(ASSOCIATED("."$string".")) "."$string"."(i1:i2,j1:j2,k) = var4dDBL(i1:i2,j1:j2,kReverse,"."$prefix$n".")\n";
   } else {
    $line = " IF(ASSOCIATED("."$string".")) "."$string"."(i1:i2,j1:j2) = var3d(i1:i2,j1:j2,"."$prefix$n".")\n";
   }

# Append content to file
# ----------------------
   print FILE $line;
  }
 }

# Done. Close file.
# -----------------
 close(FILE);
}

sub write_q_h_files {

# The four streams
# ----------------
 $stream[0] = "QK";
 $stream[1] = "QQK";
 $stream[2] = "QJ";
 $stream[3] = "QQJ";

# Choose file name and open
# -------------------------
 $fileName = "GMI_GridComp/$stream[$i]"."_FillExports___.h";
 open FILE, ">", $fileName or die "Cannot open $fileName\n";
 print "Building $fileName\n";

# Terminal index is $num_k or $num_j
# ----------------------------------
 $lastM = $num_k;
 if ( $i > 1 ) { $lastM = $num_j; }

# Cycle through the equations
# ---------------------------
 for $m ( 0..$lastM - 1 ) {
  $n = $m + 1;
  
# Spacers for tidiness
# --------------------
  $prefix = "";
  if ( $n <  100 ) { $prefix = "0"  ; }
  if ( $n <   10 ) { $prefix = "00" ; }
  $spaces = "";
  if ( $n <  100 ) { $spaces = " "  ; }
  if ( $n <   10 ) { $spaces = "  " ; }

# Pointer name
# ------------
  $string = "$stream[$i]"."$prefix"."$n";

# Complete line to be added to file
# ---------------------------------
  $line = "    IF(ASSOCIATED("."$string".")) "."$string"."(i1:i2,j1:j2,1:km) = gmi"."$stream[$i]"."("."$spaces"."$n".")%pArray3D(i1:i2,j1:j2,km:1:-1)\n";

# Append content to file
# ----------------------
   print FILE $line
 }

# Done. Close file.
# -----------------
 close(FILE);
}

sub find_depos {

open FILE, $deposParamFile or die "Cannot find $deposParamFile\n";
# print "Opened $deposParamFile\n";

# These constants are unique to setkin_depos.h and 
# must be modified if the format of the table is changed.
# -------------------------------------------------------
$colsInTable = 5;
$firstChar = 8;
$valueLength = 10;
$filler = 3;

# The number of lines to parse
# ----------------------------
use integer;
$numLines = 1 + ( $numSpecies-1 ) / $colsInTable;

# Number of columns in the last line
# ----------------------------------
$lc = $numSpecies - $colsInTable*($numLines-1) -1 ;
no integer;

# Initialize
# ----------
$proc = 0;
$n = -1000;
$numDepos = 0;
$sn = -1;
$zeroString = " 0.000D+00";

# Read file line-by-line
# ----------------------
while ( $line = <FILE> ) {
 
# Find "data hstar" and process the next $numLines lines
# ------------------------------------------------------
 $string = "data $hstar_string";
 $rc = index( $line, $string );
 if ( $rc > 0 ) { $n = 0; $proc = 1; } 

# While the processing flag is on
# -------------------------------
 if ( $n > 0 && $proc == 1 ) {

# For safety, be explicit about columns of data to examine
# --------------------------------------------------------
  if ( $n == $numLines ) { $lastColumn = $lc; } else { $lastColumn = $colsInTable - 1; }

# Column by column
# ----------------
  for $col ( 0..$lastColumn ) {

# Current species number
# ----------------------
   $sn++;

# Extract value of hstar
# ----------------------
   $firstCharNum = $firstChar + ( $valueLength + $filler ) * $col;
   $value = substr $line, $firstCharNum, $valueLength;

# Store the species number is a non-zero value is found
# -----------------------------------------------------
   if ( $value != $zeroString ) { 
    $numDepos++;
    $deposSpeciesNumber[$numDepos-1] = $sn;
   }

# Next column on this line
# ------------------------
  }

# After the designated number of lines
# has been processed, turn the flag off.
# --------------------------------------
  if ( $n >= $numLines ) { $proc = 0; }
 }

# Increment line number
# ---------------------
 $n++;
}

close(FILE);

}
sub register_exports {
# -------------------------------------------------
# Create individual export specification files for:
# 1. Tendencies
# 2. Scavenging and dry and wet deposition
# 3. Reaction rate constants and reaction rates
# -------------------------------------------------

# Start out with the rates of change registries
# ---------------------------------------------
 $stream[0] = "GMI";
 $stream[1] = "PHYS";
 $stream[2] = "AGCM";
 $stream[3] = "IM";
 $stream[4] = "IT";
 $component[0] = "chemistry";
 $component[1] = "physics";
 $component[2] = "dynamics";
 $component[3] = "moist";
 $component[4] = "turbulence";

 for $n ( 0..0 ) {
  $fileName = "$stream[$n]"."_Tendency_Registry___.rc";
  open NF, ">", $fileName or die "Cannot open $fileName\n";
  print "Building "."$fileName\n";
  export_spec_header();

  $dim = "xyz";
  $vloc = "C";
  for $m ( 0..$num_internal_state-1 ) {

   $longName  = "$component[$n]"."_rate_of_change_"."$internal_spec_name[$m]";
   $suffix    = "       ";
   $unit      = "$internal_spec_unit[$m]"." s-1 ";

   if ( $n < 3 ) {
    $shortName = "$internal_spec_name[$m]"."_"."$stream[$n]"."TEND";
      } else {
    $shortName = "$internal_spec_name[$m]"."$stream[$n]";
   }

   $keep = 0;
   if ( $friends[$m] eq "D,T,C" ) { $keep = 1; }

   $line = "  "."$shortName"."$suffix"." | "."$unit"." | "."$dim"." | "."$vloc"." |    |   |   |     | "."$longName\n";
   if ( $keep == 1 ) {
    $mLast = $m;
    print NF $line;
   }
  }

  $line = "</ExportSpec>\n";
  print NF $line;
  close(NF);
 }

# Scavenging and dry and wet deposition
# -------------------------------------
 $fileName = "Deposition_Registry___.rc";
 open NF, ">", $fileName or die "Cannot open $fileName\n";
 print "Building "."$fileName\n";
 export_spec_header();

# Each process, scav, dd, wd
# --------------------------
 for $i ( 0..2 ) {

    if ( $i == 0 ) { $hstar_string = "hstar_wet"; }
    if ( $i == 1 ) { $hstar_string = "hstar_dry"; }
    if ( $i == 2 ) { $hstar_string = "hstar_wet"; }
#   print "Setting hstar to $hstar_string \n";

    find_depos();

# Cycle through the species
# -------------------------
  for $m ( 0..$numSpecies-1 ) {
   $n = $m + 1;

# Scan the list of desired exports and continue only if species(m) found in list
# ------------------------------------------------------------------------------
   $match = 0;
   for $d ( 0..$numDepos-1 ) {
    $deposName = $species[$deposSpeciesNumber[$d]];
    if ( $deposName eq $species[$m] ) { $match = 1; }
   }

# Concatenate species name to process abbreviation, right justified
# -----------------------------------------------------------------
   if ( $match == 1 ) {
    $shortName = "$process[$i]"."$species[$m]"."$blankFill[$m]";
    if ( $i == 0 ) { 
     $dim = "xyz";
     $vloc = "C";
     $suffix = "    ";
     $unit = "kg m-3 s-1    ";
    } else {
     $dim = "xy ";
     $vloc = " ";
     $suffix = "      ";
     $unit = "kg m-2 s-1    ";
    }
    if ( $i == 0 ) { $longName = "scavenging_of_"."$species[$m]"; }
    if ( $i == 1 ) { $longName = "dry_deposition_of_"."$species[$m]"; }
    if ( $i == 2 ) { $longName = "wet_deposition_of_"."$species[$m]"; }

# Complete line to be added to file
# ---------------------------------
    $line = "  "."$shortName"."$suffix"." | "."$unit"." | "."$dim"." | "."$vloc"." |    |   |   |     | "."$longName\n";

# Append content to file
# ----------------------
    print NF $line;
   }
  }
 }

 $line = "</ExportSpec>\n";
 print NF $line;
 close(NF);

# Reaction rate constants and reaction rates
# ------------------------------------------
 $fileName = "Reactions_Registry___.rc";
 open NF, ">", $fileName or die "Cannot open $fileName\n";
 print "Building "."$fileName\n";
 export_spec_header();

# Each process, scav, dd, wd
# --------------------------
 $dim = "xyz";
 $vloc = "C";

# Cycle through the thermal equations
# -----------------------------------
 for $m ( 0..$num_k-1 ) {
  $n = $m + 1;

# A spacer for tidiness
# ---------------------
  $prefix = "";
  if ( $n <  100 ) { $prefix = "0"  ; }
  if ( $n <   10 ) { $prefix = "00" ; }
  
# Add two lines for each equation
# -------------------------------
  $shortName = "QK$prefix$n";
  $longName  = "rate_constant: $kReaction[$m]";
  $unit      = "2-3body_varies";
  $suffix    = "          ";
  $line = "  "."$shortName"."$suffix"." | "."$unit"." | "."$dim"." | "."$vloc"." |    |   |   |     | "."$longName\n";
  print NF $line;

  $shortName = "QQK$prefix$n";
  $longName  = "reaction_rate: $kReaction[$m]";
  $unit      = "mole m-3 s-1  ";
  $suffix    = "         ";
  $line = "  "."$shortName"."$suffix"." | "."$unit"." | "."$dim"." | "."$vloc"." |    |   |   |     | "."$longName\n";
  print NF $line;
 }

# Cycle through the photolysis equations
# --------------------------------------
 for $m ( 0..$num_j-1 ) {
  $n = $m + 1;

# A spacer for tidiness
# ---------------------
  $prefix = "";
  if ( $n <  100 ) { $prefix = "0"  ; }
  if ( $n <   10 ) { $prefix = "00" ; }
  
# Add two lines for each equation
# -------------------------------
  $shortName = "QJ$prefix$n";
  $longName  = "rate_constant: $jReaction[$m]";
  $unit      = "s-1           ";
  $suffix    = "          ";
  $line = "  "."$shortName"."$suffix"." | "."$unit"." | "."$dim"." | "."$vloc"." |    |   |   |     | "."$longName\n";
  print NF $line;

  $shortName = "QQJ$prefix$n";
  $longName  = "reaction_rate: $jReaction[$m]";
  $unit      = "mole m-3 s-1  ";
  $suffix    = "         ";
  $line = "  "."$shortName"."$suffix"." | "."$unit"." | "."$dim"." | "."$vloc"." |    |   |   |     | "."$longName\n";
  print NF $line;
 }

 $line = "</ExportSpec>\n";
 print NF $line;
 close(NF);

}

sub export_spec_header {
# -----------------------------
# Write export spec file header
# -----------------------------
 $line = "  COMP_NAME: GMICHEM\n";
 print NF $line;
 $line = "  MAPL_REGISTRY_VERSION: 1.00\n";
 print NF $line;
 $line = "<ExportSpec name=\"GMICHEM\", cols=\"short_name,units,dims,vlocation,stat,refresh_interval,averaging_interval,num_subtiles,long_name\">\n";
 print NF $line;

}
sub usage {
   print <<"EOF";

NAME
     gmi_acg.pl - Generates code fragments for GMIchem_GridComp and GMI_GridComp.
          
SYNOPSIS

     gmi_acg.pl [OPTIONS]
          
DESCRIPTION

     Generates the following #include files for GMIchem and GMI grid components:

     Deposition_DeclarePointer___.h
     Deposition_ExportSpec___.h
     Deposition_GetPointer___.h
     Deposition_Registry___.rc
     SCAV_FillExports___.h
     DD_FillExports___.h
     WD_FillExports___.h
     QK_FillExports___.h
     QQK_FillExports___.h
     QJ_FillExports___.h
     QQJ_FillExports___.h
     Reactions_DeclarePointer___.h
     Reactions_ExportSpec___.h
     Reactions_GetPointer___.h
     Reactions_Registry___.rc
     Tendency_DeclarePointer___.h
     Tendency_ExportSpec___.h
     Tendency_FillExports___.h
     Tendency_GetPointer___.h
     Tendency_InitialState___.h
     Tendency_Registry___.rc

     These files are ignored by the CVS ESMA repository at NASA/GSFC and are not to 
     be checked in.

     gmi_acg.pl works in tandem with, but independently from, mapl_acg.pl to account 
     for properties unique to GMICHEM.  A multiple-step process is required to create
     the above files:
     
     Command
     ------------------------------------------------------
     gmi_acg.pl -R
     mapl_acg.pl -N GMICHEM Deposition_Registry___.rc
     mapl_acg.pl -N GMICHEM Reactions_Registry___.rc
     mapl_acg.pl -N GMICHEM Tendency_Registry___.rc
     mapl_acg.pl -N GMICHEM GMI_GridComp/GMI_Registry___.rc
     gmi_acg.pl -v

     After each of the first two uses of mapl_acg.pl, move the files.  For example:
      /bin/mv -f GMICHEM_DeclarePointer___.h GMI_GridComp/Deposition_DeclarePointer___.h
      /bin/mv -f GMICHEM_GetPointer___.h GMI_GridComp/Deposition_GetPointer___.h
      /bin/mv -f GMICHEM_ExportSpec___.h Deposition_ExportSpec___.h
     
     After the third use of mapl_acg.pl, rename the files:
      /bin/mv -f GMICHEM_DeclarePointer___.h Tendency_DeclarePointer___.h
      /bin/mv -f GMICHEM_GetPointer___.h Tendency_GetPointer___.h
      /bin/mv -f GMICHEM_ExportSpec___.h Tendency_ExportSpec___.h

     After the fourth use of mapl_acg.pl, these files can be removed:
      /bin/rm GMICHEM_GetPointer___.h  GMICHEM_DeclarePointer___.h

     After the second use of gmi_acg.pl, restore names for GetPointer files:
      /bin/mv -f GMI_GridComp/Deposition_GetPointer2___.h GMI_GridComp/Deposition_GetPointer___.h
      /bin/mv -f GMI_GridComp/Reactions_GetPointer2___.h GMI_GridComp/Reactions_GetPointer___.h
      /bin/mv -f GMICHEM_ExportSpec___.h Tendency_ExportSpec___.h
     
     In addition to GMI_Registry.rc, gmi_acg.pl parses the following text files:
     
     File contents          Default file name     Where to find
     ---------------------  --------------------  --------------------------
     Chemical mechanism     setkin_chem_mech.txt  GmiChemistry/stratTropMech
     Deposition parameters  setkin_depos.h        GmiChemistry/stratTropMech
     Mechanism parameters   setkin_par.h          GmiChemistry/stratTropMech

OPTIONS

     -d   Deposition parameters file name. Default: See above.
     -H   Get this help page.
     -m   Mechanism text file name. Default: See above.
     -p   Mechanism parameters file name.  Default: See above.
     -R   Create registry files.
     -r   Registry.  Default: GMI_GridComp/GMI_Registry.rc
     -v   Verbose mode.

EOF
exit(1)
}
