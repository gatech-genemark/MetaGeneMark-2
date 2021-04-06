#!/usr/bin/perl
# -----------------------
# Alex Lomsadze
# GeorgiaTech 2020
#
# Predict genes in metagenome
# Discriminate between sequences of genetic code 11 and 4
# -----------------------

use warnings;
use strict;

use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use Cwd qw(abs_path cwd);
use File::Temp qw(tempfile tempdir);

# -----------------------
# algorithm thresholds

# All contigs are assigned to genetic code 11 by default.
# After that genetic code 4 is tested.
# Two tests are used for genetic code change.
#   by average gene length
# and
#   by logodd score

# by average gene length
# apply average gene length criteria to contigs longer than this minimum
my $MIN_CONTIG_LENGTH_FOR_GCODE_BY_AVERAGE = 50000;
# it is genetic code 11 if average length is above this minimum
my $MIN_GCODE_11_AVERAGE = 700;
# for shorter average length use this separation line
my $GCODE_CONST_A = 1.33;
my $GCODE_CONST_B = 20;

# by logodd score
# apply logood criteria to contigs longer than this minimum
my $MIN_CONTIG_LENGTH_FOR_GCODE_BY_LOGODD = 5000;
# minimum difference in logood value to call genetic code 4
my $MIN_LOGODD_DELTA_FOR_4 = 10;
# minimum logodd score to call not 11
my $MIN_LOGODD_FOR_NOT_11 = 20;

# -----------------------
my $VERSION = "1.00";

# required
my $seq_file = '';
my $out_file = '';

# optional
my $format = "gtf";
my $nt_file = '';
my $aa_file = '';
my $clean = '';
my $tmpf = '';
my $verbose = 0;
my $debug = 0;
my $pf_summary = "summary.csv";

# local
my $work_dir = cwd;
my $command = '';

# -----------------------

Usage() if @ARGV < 1;
ParseCMD();
CheckBeforeRun();

$tmpf = SetTmpFolder($tmpf);

# temporary files genetic code 11
my $out_11 = $tmpf ."/mgm_11.out";
my $nt_11  = $tmpf ."/mgm_11.nt";
my $aa_11  = $tmpf ."/mgm_11.aa";

$command = "$Bin/gmhmmp2  -s $seq_file -f $format  -M $Bin/mgm2_11.mod  --gid_per_contig";
$command .= " --out $out_11";
$command .= " --NT  $nt_11" if $nt_file;
$command .= " --AA  $aa_11" if $aa_file;

print "# predicing with genetic code 11 ...\n" if $verbose;
RunSystem( $command );

# temporary files genetic code 4
my $out_4 = $tmpf ."/mgm_4.out";
my $nt_4  = $tmpf ."/mgm_4.nt";
my $aa_4  = $tmpf ."/mgm_4.aa";

$command = "$Bin/gmhmmp2  -s $seq_file  -f $format -M $Bin/mgm2_4.mod  --gid_per_contig";
$command .= " --out $out_4";
$command .= " --NT  $nt_4" if $nt_file;
$command .= " --AA  $aa_4" if $aa_file;

print "# predicing with genetic code 4 ...\n" if $verbose;
RunSystem( $command );

# compare predictions and select the best
my %data_11 = ();
my %data_4 = ();
my %data = ();

LoadMetaData( $out_11, \%data_11 );
LoadMetaData( $out_4, \%data_4 );
SelectBestGeneticCode( \%data_11, \%data_4, \%data );
MergeMGM_Predictions( \%data, $out_11, $out_4, $out_file );

# `echo "p4,p11,pf" > $pf_summary`;
# 
# 
# my $p4;
# my $p11;
# for ($p4 = -100; $p4 < 100; $p4 += 5) {
#     for ($p11 = -100; $p11 < 100; $p11 += 5) {
#         #%data_11 = ();
#         #%data_4 = ();
#         %data = ();
#         #LoadMetaData( $out_11, \%data_11 );
#         #LoadMetaData( $out_4, \%data_4 );
#         $MIN_LOGODD_DELTA_FOR_4 = $p4;
#         # minimum logodd score to call not 11
#         $MIN_LOGODD_FOR_NOT_11 = $p11;
# 
#         my $newout = "$out_file" . "_" . $p4 . "_" . $p11;
#         SelectBestGeneticCode( \%data_11, \%data_4, \%data );
#         MergeMGM_Predictions( \%data, $out_11, $out_4, $newout);
# 
#         `echo "$p4,$p11,$newout" >> $pf_summary`;
#     }
# }


# MergeMGM_FASTAs( \%data, $nt_11, $nt_4, $nt_file ) if $nt_file;
# MergeMGM_FASTAs( \%data, $aa_11, $aa_4, $aa_file ) if $aa_file;

CleanTMP() if $clean;

print "# Done\n" if $verbose;

# ---------------------
sub CleanTMP
{
	unlink $out_11 if ( -e $out_11);
	unlink $out_4  if ( -e $out_4);
	unlink $nt_11  if ( -e $nt_11);
	unlink $aa_11  if ( -e $aa_11);
	unlink $nt_4   if ( -e $nt_4);
	unlink $aa_4   if ( -e $aa_4);

	rmdir $tmpf;
}
# ----------------------
sub MergeMGM_FASTAs
{
	;
}
# ----------------------
sub MergeMGM_Predictions
{
	my $ref = shift;
	my $fname_11 = shift;
	my $fname_4 = shift;
	my $fname = shift;

	print "# creating output ...\n" if $verbose;

	open( my $OUT, ">", $fname ) or die "error on open file $fname\n";

	open( my $IN_11, $fname_11 ) or die "error on open file $fname_11\n";
	while(<$IN_11>)
	{
		if ( /^#/ )
		{
			if ( /^##gff-version/ )
			{
				print $OUT $_;
			}
			elsif ( /^# GeneMark.hmm/ )
			{
				print $OUT $_;
			}
			elsif ( /^# File with sequence:/ )
			{
				print $OUT $_;
			}
			elsif ( /^# File with MetaGeneMark parameters:/ )
			{
				; # add both 11 and 4 ?
			}
			elsif ( /^# translation table:/ )
			{
				; # thid is nou a feature of each contig
			}
			elsif ( /# output date start:/ )
			{
				print $OUT $_;
			}
			elsif ( /^##sequence-region\s+(\S+)\s*/ )
			{
				my $id = $1;
				my $line = $_;

				if ( exists $ref->{$id} )
				{
					if ( $ref->{$id} == 11 )
					{
						print $OUT "\n";
						print $OUT $_;
						print $OUT "# seqid: $id  genetic code: 11\n";
						print $OUT "\n";
					}
					elsif ( $ref->{$id} == 0 )
					{
						print $OUT "\n";
						print $OUT $_;
						print $OUT "# seqid: $id  genetic code: 0\n";
						print $OUT "\n";
					}
				}
				else
					{ die "error, contig not found in genetic code classifier: $id\n"; }
			}
			elsif ( /^#\s+(\S+)\s+total_logodd/ )
			{
				my $id = $1;

				if ( exists $ref->{$id} )
				{
					if ( $ref->{$id} == 11 )
					{
						print $OUT "\n";
						print $OUT $_;
					}
					elsif ( $ref->{$id} == 0 )
					{
						print $OUT $_;
					}
				}
			}
			else
				{ die "error, unexpected line found: $_"; }
		}
		else
		{
			if ( /^(\S+)/ )
			{
				my $id = $1;

				print $OUT $_ if ($ref->{$id} == 11);
				print $OUT $_ if ($ref->{$id} == 0);
			}		
		}
	}
	close $IN_11;

	open( my $IN_4, $fname_4 ) or die "error on open file $fname_4\n";
	while(<$IN_4>)
	{
		if ( /^#/ )
		{
			if ( /^##sequence-region\s+(\S+)\s*/ )
			{
				my $id = $1;
				my $line = $_;

				if ( exists $ref->{$id} )
				{
					if ( $ref->{$id} == 4 )
					{
						print $OUT "\n";
						print $OUT $_;
						print $OUT "# seqid: $id  genetic code: 4\n";
						print $OUT "\n";
					}
				}
				else
					{ die "error, contig not found in genetic code classifier: $id\n"; }
			}
			elsif ( /^#\s+(\S+)\s+total_logodd/ )
			{
				my $id = $1;

				if ( exists $ref->{$id} )
				{
                                        if ( $ref->{$id} == 4 )
					{
						print $OUT "\n";
						print $OUT $_;
					}
				}
			}
		}
                else
                {
			if ( /^(\S+)/ )
			{
				my $id = $1;
				print $OUT $_ if ($ref->{$id} == 4);
			}
                }
	}
	close $IN_4;

	close $OUT;
}
# ----------------------
sub BestByLogodd
{
	my $score_11 = shift;
	my $score_4 = shift;
	my $contig_length = shift;

	if ( $score_11 == 0 and $score_4 == 0 )
	{
		return 0;
	}

	if ( $contig_length < $MIN_CONTIG_LENGTH_FOR_GCODE_BY_LOGODD )
	{
		return 11;
	}	

	my $best_code = 11;

	if (( $score_4 - $score_11 > $MIN_LOGODD_DELTA_FOR_4 )and( $score_4 > $MIN_LOGODD_FOR_NOT_11 ))
	{
		$best_code = 4;
	}

	return $best_code;
}
# ----------------------
sub BestByAverage
{
	my $score_11 = shift;
	my $score_4 = shift;
	my $contig_length = shift;

	if ( $score_11 == 0 and $score_4 == 0 )
	{
		return 0;
	}

	if( $contig_length < $MIN_CONTIG_LENGTH_FOR_GCODE_BY_AVERAGE )
	{
		return 11;
	}

	if ( $score_11 > $MIN_GCODE_11_AVERAGE )
	{
		return 11;
	}

	my $best_code = 11;

	if ( $score_4 > $GCODE_CONST_A * $score_11 + $GCODE_CONST_B )
	{
		$best_code = 4;
	}

	return $best_code;
}
# ----------------------
sub SelectBestGeneticCode
{
	my $ref_11 = shift;
	my $ref_4 = shift;
	my $ref = shift;

	print "# determining genetic code ...\n" if $verbose;

	foreach my $key (keys %{$ref_11})
	{
		if ( exists $ref_4->{$key} )
		{
			my $gcode_by_average = BestByAverage( $ref_11->{$key}{"average"}, $ref_4->{$key}{"average"}, $ref_11->{$key}{"size"} );
			my $gcode_by_logodd  = BestByLogodd ( $ref_11->{$key}{"logodd"},  $ref_4->{$key}{"logodd"},  $ref_11->{$key}{"size"} );

			if ( $gcode_by_average ne $gcode_by_logodd )
			{
				print "# warning, difference in calls: by_average $gcode_by_average ne by_logodd $gcode_by_logodd\n" if $debug;
			}
			
			$ref->{$key} = $gcode_by_logodd;
		}
		else
			{ die "error, FASTA ID is not found in the second hash: $key\n"; }
	}

	if ($verbose)
	{
		my %h = ();

		$h{'0'} = 0;
		$h{'4'} = 0;
		$h{'11'} = 0;

		foreach my $key (keys %{$ref})
		{
			$h{ $ref->{$key} } += 1;
		}
		

		foreach my $key (keys %h)
		{
			 print "# gcode $key assigned to $h{$key} contigs\n";
		}
	}
}
# ----------------------
sub LoadMetaData
{
	my $fname = shift;
	my $ref = shift;

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while(<$IN>)
	{
		if ( /##sequence-region\s+(\S+)\s+1\s+(\d+)\s*/ )
		{
			$ref->{$1}{"size"}    = $2;
			$ref->{$1}{"logodd"}  = 0;
			$ref->{$1}{"average"} = 0;
		}
		elsif ( /^#\s+(\S+)\s+total_logodd\s+(\S+)\s+average_length\s+(\S+)\s*/ )
		{
			$ref->{$1}{"logodd"}  = $2;
			$ref->{$1}{"average"} = $3;
		}
	}
	close $IN;

	if ($debug)
	{
		foreach my $key (keys %{$ref})
		{
			die "error, missing contig size for $key" if ( ! exists $ref->{$key}{"size"} );
			die "error, missing logodd for $key"      if ( ! exists $ref->{$key}{"logodd"} );
			die "error, missing average for $key"     if ( ! exists $ref->{$key}{"average"} );
		}

		my $size = scalar (keys %{$ref});
		print "# $size contigs in file $fname\n";
	}
}
# ----------------------
sub RunSystem
{
	my $command = shift;
	print "# $command\n" if $debug;
	system( $command ) and die "error on last system call: $command\n";
}
# ----------------------
sub CreateTmpWorkspace
{
	my $dir = shift;
	my $label = shift;

	my $tmpdir = tempdir( $label ."_XXXXX", DIR => $dir );
	if ( ! -d  $tmpdir ) { die "error, temporary folder not found: $tmpdir\n"; }

	return $tmpdir;
}
# ----------------------
sub SetTmpFolder
{
	my $fname = shift;

	if ( $fname )
	{
		if ( ! -d $fname )
		{
			mkdir $fname or die "error on create temporary folder $fname\n";
		}
        }
	else
	{
		$fname = CreateTmpWorkspace( $work_dir, "mgm" );
        }

	return $fname;
}
# ----------------------
sub CheckBeforeRun
{
	print "# checking before run\n" if $debug;

	if( !$seq_file ) { die "error, required input file name is missing, check option --seq\n"; }
	if( !$out_file ) { die "error, required output file name is missing, check option --out\n"; }

	if (( $format ne "gtf" )and( $format ne "gff3" ))
		{ die "error, unexpected value was specified for option --format: $format\n"; }
}
# ----------------------
sub ParseCMD
{
	my $cmd = $0;
	foreach my $str (@ARGV) { $cmd .= ( "\t". $str ); }

	my $opt_results = GetOptions
	(
		'seq=s'     => \$seq_file,
		'out=s'     => \$out_file,
		'nt=s'      => \$nt_file,
		'aa=s'      => \$aa_file,
		'format=s'  => \$format,
		'clean'     => \$clean,
		'tmpf=s'    => \$tmpf,
		'verbose'   => \$verbose,
		'debug'     => \$debug,
        'p4=i'        => \$MIN_LOGODD_DELTA_FOR_4,
        'p11=i'       => \$MIN_LOGODD_FOR_NOT_11,
        'pf-summary=s' => \$pf_summary
	);

	if( !$opt_results ) { die "error on command line $cmd\n"; }
	if( @ARGV > 0 ) { die "error, unexpected argument found on command line: $cmd\n"; }

	$verbose = 1 if $debug;

	print "# command line parsing ... done\n" if $debug;
}
# -----------------------
sub Usage
{
	print qq(# -------------------
Usage:  $0  --seq [name]  --out [name]

Required options:

  --seq  [name]
     nucleotide sequence of metagenome in FASTA format.
  --out  [name]
     output file with coordinates of predicted protein coding genes.

Output options:

  --nt  [name]
     output file with nucleotide sequences of predicted genes in FASTA format.
  --aa  [name]
     output file with protein sequences of predicted genes in FASTA format.
  --format  [$format]
     format of output file with gene coordinates: gtf or gff3.
  --clean
     delete temporay files

Other parameters:
  --verbose
Developer options:
  --tmpf  [name]
     folder name for temporary files
  --debug

Version  $VERSION
# -------------------
);
         exit 1;
}
# -----------------------

