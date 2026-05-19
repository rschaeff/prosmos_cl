#!/usr/bin/perl
use strict;
use warnings;
#this is a update version, transform one matrix file in one dir

my $numargs=@ARGV;

if($numargs<2){

	print "transform genMatrix file\n";
	#print "this is a simple version, only transform one matrix file containg one pdb\n";
	print "Usage:perl transmatrix.pl [matrixfile][resultfile]\n";
	exit;
}

my $matrixfile=$ARGV[0];
my $resultfile=$ARGV[1];

		print "parsing $matrixfile\n";
		unless(open(MATRIX,"$matrixfile")){
			die("Can not open file $matrixfile\n");
      exit;
		}
		unless(open(RESULT,">>$resultfile")){
	 		die("Can not open file $resultfile\n");
      exit;
		}
		while(my $line=<MATRIX>){
			chop($line);
			if($line=~/sheet/){
				print RESULT $line."\n";
				next;
			}
  		elsif($line=~/\*/){
  			my @array=split(/\*/,$line);
  			push(@array,"*");
  			my $arraycount=@array;
  			for(my $i=1;$i<@array;$i++){
  				my $blank="";	
  				for(my $j=1;$j<$i;$j++){
  					$blank=$blank." ";
  				}
  				if($i==($arraycount-1)){
  					print RESULT $blank.$array[$i]."\n";
  				}else{	
  					print RESULT $blank."*".$array[$i]."\n";
  				}	
  			}
  			
  		}elsif($line=~/^$/){
  
  		}else{
  			my $startline=substr($line,0,36);
  			print RESULT $startline."\n";
  			my @array=split(/\s/,$startline);
  			my $count=@array;
  			my $elenumber=$array[$count-1];
  			my $startnumber=36;
  			my $addnumber=68;
  			for(my $i=0;$i<$elenumber;$i++){
  				$startnumber=36+$i*$addnumber;
  				my $eleline=substr($line,$startnumber,$addnumber);
 					print RESULT $eleline."\n";
    		}	
  
  		}	
		}
		close(MATRIX);
		close(RESULT);

print "transform metamatrix to formated metamatrix finished\n";
exit;