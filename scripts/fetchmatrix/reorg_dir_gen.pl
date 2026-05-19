#!/usr/bin/perl
#reorg_dir_gen.pl -- 30,000 file directories are an NFS deathtrap, rearrange some directory of PDB-like outputs into a "two" structure
#i.e. 1enh.blah goes in /en/1enh.blah.

use warnings;
use strict;
use File::Copy;

$ARGV[0] or die "usage: reorg_dir_gen.pl <one column list of some files you'd like broken up into manageble chunks>\n";

open (IN, $ARGV[0]) or die "Could not open $ARGV[0] for reading:$!\n";

while (my $ln = <IN>) { 
	$ln =~ /^#/ and next;
	my @F = split(/\s+/, $ln);

	
	if (-f $F[0]) { 
		if ($F[0] =~ /\.\/([\w\/]+)\/\w(\w{2})\w/) { 
			my $dir_junk = $1;
			my $two = $2;
			if (-d "./$dir_junk/$two") { 
				move("$F[0]", "./$dir_junk/$two/");
				print "move $F[0] ./$dir_junk/$two/\n";
			}else{
				if (-d "./$dir_junk") {
					mkdir("./$dir_junk/$two");
					print "mkdir ./$dir_junk/$two\n";
					move("$F[0]", "./$dir_junk/$two/");
					print "move $F[0] ./$dir_junk/$two/\n";
				}
			}
		}else{
			print "$F[0] doesn't look like it has a PDB code in it, skipping\n";
		}

	}else{
		print "$F[0] not found, skipping\n";
	}
}
