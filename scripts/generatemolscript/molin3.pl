#! /usr/bin/perl -w

if( $#ARGV != 0 ){
	print STDERR "ERROR: need to input the pdb file name\n";
	exit;
}

$ARGV[0] =~ /.pdb/;
$pdbid = $`;

system( "molauto -ss_pdb $ARGV[0] > .$pdbid.in" );

open( IN, ".$pdbid.in" ) || die "ERROR: can't open input file .$pdbid.in";
open( OUT, ">$pdbid.in" ) || die "ERROR: can't open output file $pdbid.in";
@coils = ();
@strands = ();
@helices = ();
@others = ();
while( <IN> ){
	chomp;
	if( /^!/ ){
		print OUT "$_\n";
		next;
	}elsif( /^plot/ ){
		print OUT "$_\n\n";

		print OUT "window 42;\n";
		print OUT "!open in photoshop at 300\n";
		print OUT "! window 49;\n";
		print OUT "!open in photoshop at 350\n";
		print OUT "! window 56;\n";
		print OUT "!open in photoshop at 400\n\n";

		print OUT "!Photoshop manipulations\n";
		print OUT "!transfer to the field of 650\n";
		print OUT "!Image Adjust Curves: Pull the line up through the 4th dot in the 3rd square\n";
		print OUT "!N and C for termini: Black Helvetica Regular 15\n";
		print OUT "!a,b,c for panels: Black Helvetica Bold Oblique 25\n";
		print OUT "!A,B,C,a,b,c for helices and strands: Helvetica Bold Oblique 10\n\n";

		print OUT "set segments 50;\n";
		print OUT "!set colourdepthcue 0.5;\n";
		print OUT "set shading        0.5;\n";
		print OUT "set shadingexponent 1.0;\n";
		print OUT "set splinefactor  1.0;\n";
		print OUT "set strandthickness 0.4;\n";
		print OUT "set strandwidth   2.0;\n";
		print OUT "set smoothsteps    5;\n";
		print OUT "set coilradius   0.25;\n";
		print OUT "set linewidth   1.5;\n\n";
		next;
	}elsif( /read mol/ ){
		print OUT "$_\n";
		next;
	}elsif( /transform atom/ ){
		chop;
		print OUT "$_\n";
		print OUT "by rotation x  0\n";
		print OUT "by rotation y 0\n";
		print OUT "by rotation z  0\n";
		print OUT "by translation 0 0 0\n;\n\n"; 
		next;
	}elsif( /set segments 2/ ){
	 	print OUT "! the white parts\n";
		print OUT "set helixwidth 1.4;\n\n";

		print OUT "\n!to emphasize (red)\n";
		print OUT "set planecolour rgb 1.0 0.2 0.2;\n";
		print OUT "set plane2colour rgb 0.5 0.1 0.1;\n\n";

		print OUT "\n!to emphasize (purple)\n";
		print OUT "set planecolour rgb 0.7 0.4 1.0;\n";
		print OUT "set plane2colour rgb 0.4 0.2 0.5;\n\n";
	}

	if( /set planecolour hsb/ ){
		next;
	}elsif( /coil from .* to .*/ ){
		@coils = ( @coils, $_ );
		next;
	}elsif( /turn from .* to .*/ ){
		@coils = ( @coils, $_ );
		next;	
	}elsif( /strand from .* to .*/ ){
		@strands = (@strands, $_ );
		next;
	}elsif( /helix from .* to .*/ ){
		@helices = ( @helices, $_ );
		next;
	}elsif( (not /end_plot/) && (not /^$/) ){
		@others = ( @others, $_ );
		next;
	}
}
		
print OUT "!the strands\n";
print OUT "set planecolour rgb 0.972549 0.972549 1.0;\n";
print OUT "set plane2colour rgb 0.486275 0.486275 0.5;\n";
print OUT "set strandwidth 1.8;\n";
print OUT "set strandthickness 0.4;\n";
foreach $strand (@strands){
	print OUT "$strand\n";
}


print OUT "\n!the helices\n";
print OUT "set planecolour rgb 0.972549 0.972549 1.0;\n";
print OUT "set plane2colour rgb 0.486275 0.486275 0.5;\n";
print OUT "set helixwidth 1.9;\n";
foreach $helix (@helices){
	print OUT "$helix\n";
}

print OUT "\n!the coils\n";
print OUT "set planecolour rgb 0.65 0.85 0.65;\n";
print OUT "set smoothsteps 10;\n";
foreach $coil (@coils){
	print OUT "$coil\n";
}
print OUT "\n";

foreach $other (@others){
	print OUT "$other\n";
}	


print OUT "!colour ss from blue via green to red\n";
print OUT "!by b-factor from 0.3 to 1.8;\n\n";

print OUT "!set colourparts on;\n";
print OUT "!set planecolour rgb 0.2 0.2 0.2;\n";
print OUT "!set atomradius atom * 1.0;\n";
print OUT "!set linewidth 1.0;\n";
print OUT "!set planecolour rgb 1.0 1.0 1.0;\n";
print OUT "!set atomcolour atom * rgb 1.0 0.5 0.0;\n";
print OUT "!  cpk in residue 502;\n";
print OUT "!  cpk in residue 503;\n\n";

print OUT "!!lines(bonds) presentation for cofactors\n";
print OUT "!set planecolour rgb 1.0 0.5 0.0;\n";
print OUT "!set planecolour rgb 1.0 1.0 1.0;\n";
print OUT "!set linecolour rgb 1.0 1.0 1.0;\n";
print OUT "!set stickradius 0.2;\n";
print OUT "!set atomradius atom * 0.3;\n";
print OUT "!!set linecolour rgb 1.0 0.6 0.0;\n";
print OUT "!set linecolour rgb 1.0 0.2 0.2;\n";
print OUT "!set linewidth 4.0;\n";
print OUT "!  bonds in residue 102;\n\n";

print OUT "!!ball-and-stick for side-chains\n";
print OUT "!set linewidth 1.5;\n";
print OUT "!set planecolour rgb 1.0 1.0 1.0;\n";
print OUT "!set linecolour rgb 0.0 0.0 0.0;\n";
print OUT "!set atomradius atom * 1.2;\n";
print OUT "!set atomcolour atom O* rgb 1.0 0.2 0.2;\n";
print OUT "!set atomcolour atom C* rgb 0.3 0.3 0.3;\n";
print OUT "!set atomcolour atom N* rgb 0.3 0.5 1.0;\n";
print OUT "!set atomcolour atom S* rgb 1.0 1.0 0.0;\n";
print OUT "!set atomcolour atom P* rgb 1.0 1.0 1.0;\n\n";

print OUT "!ball-and-stick require\n";
print OUT "!not either atom H*, atom C, atom N, atom P or atom O\n";
print OUT "!and in either\n";
print OUT "!        residue A396,\n";
print OUT "!        residue A401,\n";
print OUT "!        residue A484,\n";
print OUT "!or      residue A518;\n\n";

print OUT "!! connect 2 breaking points\n";
print OUT "!set linecolour rgb 0.54 0.59 0.67;\n";
print OUT "!set linedash 6.0;\n";
print OUT "!set linewidth 7.0;\n";
print OUT "!  line position res-atom A333 CA \n";
print OUT "!	   to position res-atom A370 CA;\n";
print OUT "!  !line position in residue A333 to position in residue A370;\n\n";

print OUT "end_plot\n";


close( IN );
close( OUT );

system( "rm -f .$pdbid.in" );
system( "grep -v title $pdbid.in  > jnk" );
system( "mv jnk $pdbid.in" );
system( "gvim $pdbid.in &" );
system( "molscript -gl -in $pdbid.in" );

