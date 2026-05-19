#!/usr/bin/perl
use strict;
use warnings;

my $numargs=@ARGV;
if($numargs<6){

	print "Prepare pdb and log file for ProSMoS pdb hits .\n";
	print "Usage:perl generateinsightIIlog.pl [option][pdbhitslist][ProSMoS result dir][elenumber][pdbdir][resultdir]\n";
	print "option meaning:\n"; 
	print "if you want to look just 1 motif, the first motif in ProSMoS output file will be used, please set option as -s \n";
	print "if you want to look all the motifs, all the motifs in ProSMoS output file will be used, please set option as -m \n";
	print "Example:perl generateinsightIIlog.pl -s pdbhitslist ./ProSMoSoutput/ 5 ./pdbstructure/ ./insightlogoutput/ \n";
	exit;
}
my $option=$ARGV[0];
my $pdbhitslist=$ARGV[1];
my $prosmosoutputdir=$ARGV[2];
my $elenumber=$ARGV[3];
my $pdbdir=$ARGV[4];
my $resultdir=$ARGV[5];

if($option ne "-s" and $option ne "-m"){
	print "please set option as -s or -m\n";
	exit;
}	
unless(open(LIST,"$pdbhitslist")){
	 die("Can not open file $pdbhitslist\n");
   exit;
}

while(my $line=<LIST>){
	chop($line);
	if($line=~/pdb(\w+)\.txt/){
			my $targetfile=$prosmosoutputdir.$line;
			print "generate insightII log file for $targetfile\n";		
	    my $pdbid=$1;
	    my $pdbfile=$pdbdir.$pdbid.".pdb";
	    if(-e $pdbfile){
	    	my ($pdbmotifref,$allchainhashref)=readmotif($targetfile);
	    	my @pdbmotif=@$pdbmotifref;
	    	my %allchainhash=%$allchainhashref;
	    	my $atomchainhashref=readpdb($pdbid,\%allchainhash);
	    	my %atomchainhash=%$atomchainhashref;
	    	generatelog(\@pdbmotif,$pdbid,\%atomchainhash);
	    }else{
	    	print "pdb file $pdbfile does not exist in $pdbdir,please check it\n";
	    	
	    }
	    
	    
	    		
	}else{
		print "the format in your list should be pdb+pdbid+.txt,for example pdb1ubq.txt\n";
		print "please check your $pdbhitslist\n";
		exit;
	}		
}	
sub readpdb{
	my $pdbid=shift(@_);
	my $pdbfile=$pdbdir.$pdbid.".pdb";
	my $allchainhashref=shift(@_);
	my %allchainhash=%$allchainhashref;
	my %atomchainhash;
	foreach my $pdbchain(keys %allchainhash){
		unless(open(PDB,"$pdbfile")){
	 		die("Can not open file $pdbfile\n");
    	exit;
		}
		my $pdboutfile;
		if($pdbchain eq "-"){
			$pdboutfile=$resultdir.$pdbid.".pdb";
		}else{
			$pdboutfile=$resultdir.$pdbid."_".$pdbchain.".pdb";	
		}		
		
		unless(open(TMP,">>$pdboutfile")){
	 		die("Can not open file $pdboutfile\n");
    	exit;
		}
	  my $atomindex=0;
	  my %atomhash;
	  my %atomposvsatomindex;
	  my %atomindexvsatompos;
		while(my $line=<PDB>){
			chop($line);			
			if($line=~/^ATOM/ or $line=~/^HETATM/){
  
					my $nowchain=substr($line,21,1);
					if($nowchain eq " "){
						$nowchain="-";
					}	
		      if($nowchain ne $pdbchain){
						next;
					}	
				  my $indexnum=substr($line,22,5);
					if($indexnum=~/(-?\w+)/){
						$indexnum=$1;
					}	
    		  print TMP $line."\n";
    		  
	        my $type=substr($line,12,4);
				  if($type eq " CA "){
				  	my $indexinfo=$pdbchain.":".$indexnum;
				  
				  	if(defined $atomhash{$indexinfo}){
				  		next;
				  	}	
			    
				  	$atomindex++;
				  	$atomposvsatomindex{$indexnum}=$atomindex;
				  	$atomindexvsatompos{$atomindex}=$indexnum; 
					}	  
				 
	  	}	
		}	
		close(PDB);
		print TMP "END\n";
		close(TMP);
		my @tmpatominfoarray;
		push(@tmpatominfoarray,\%atomposvsatomindex);
		push(@tmpatominfoarray,\%atomindexvsatompos);
		$atomchainhash{$pdbchain}=\@tmpatominfoarray;
	}
	return(\%atomchainhash);	
}	
sub readmotif{
	my $hitfile=shift(@_);
	unless(open(HITFILE,"$hitfile")){
	 		die("Can not open file $hitfile\n");
      exit;
	}
	my @matrixarray;
	while(my $matrixline=<HITFILE>){
		chop($matrixline);
		if(($matrixline=~/pdb/) or ($matrixline=~/sub-matrix/) or ($matrixline=~/MOTIF:/) or ($matrixline=~/^END$/)){
				
		}else{
			push(@matrixarray,$matrixline);
		}		
	}	
	my @pdbmotif;

	my %allchainhash;
		for(my $putflag=0;$putflag<@matrixarray;$putflag=$putflag+$elenumber){
			my @pdbmatrix;
			my $alllength=0;
			for(my $i=$putflag;$i<$putflag+$elenumber;$i++){
				my @tempinfo=();
			
				my $startpos=substr($matrixarray[$i],36,5);
		    my $endpos=substr($matrixarray[$i],43,5);
				my $chaininfo=substr($matrixarray[$i],49,1);
				if($chaininfo eq " "){
					$chaininfo="-";
				}	
				if(defined $allchainhash{$chaininfo}){
					
				}else{
					$allchainhash{$chaininfo}=1;
				}		
				
				if($startpos=~/(-?\w+)/){
					$startpos=$1;
				}	
				if($endpos=~/(-?\w+)/){
					$endpos=$1;
				}	
				$startpos=$chaininfo."_".$startpos;
				$endpos=$chaininfo."_".$endpos;
			  if((defined $startpos) and (defined $endpos)){
				
				}else{
						print "error in $hitfile\n";
						exit;
				}		

		  
		 		push(@tempinfo,$startpos);
				push(@tempinfo,$endpos);
        push(@tempinfo,$chaininfo);
				push(@pdbmatrix,\@tempinfo);
			}
			push(@pdbmotif,\@pdbmatrix);
			if($option eq "-s"){
				last;
			}	
	  }
		close(HITFILE);
		return(\@pdbmotif,\%allchainhash);
}	

sub generatelog{
	
	my $pdbmotifref=shift(@_);
	my $pdbid=shift(@_);
	my $atomchainhashref=shift(@_);
	my %atomchainhash=%$atomchainhashref;

	my @pdbmotif=@$pdbmotifref;
	my @simplemotif;
	
  for(my $i=0;$i<@pdbmotif;$i++){
  	my $count=@simplemotif;
  	if($count>=1){
  		my $flag=0;
  		for(my $j=0;$j<@simplemotif;$j++){
  			if(($simplemotif[$j]->[0]->[0] eq $pdbmotif[$i]->[0]->[0]) and ($simplemotif[$j]->[0]->[1] eq $pdbmotif[$i]->[0]->[1]) and  ($simplemotif[$j]->[0]->[2] eq $pdbmotif[$i]->[0]->[2]) and ($simplemotif[$j]->[$elenumber-1]->[0] eq $pdbmotif[$i]->[$elenumber-1]->[0] and $simplemotif[$j]->[$elenumber-1]->[1] eq $pdbmotif[$i]->[$elenumber-1]->[1] and $simplemotif[$j]->[$elenumber-1]->[2] eq $pdbmotif[$i]->[$elenumber-1]->[2])){
  				$flag=1;
  				last;
  			}
  		}
  		if($flag==0){
  			push(@simplemotif,$pdbmotif[$i]);
  		}
  	}else{
  		push(@simplemotif,$pdbmotif[$i]);
  	}
  }
  my %logchainhash;
  for(my $ch=0;$ch<@simplemotif;$ch++){
  	my $chain=$simplemotif[$ch]->[0]->[2];
     
  	if(defined $logchainhash{$chain}){
  		my $prearrayref=$logchainhash{$chain};
  		my @prearray=@$prearrayref;
  		push(@prearray,$simplemotif[$ch]);
  		$logchainhash{$chain}=\@prearray;
  		
  	}else{
  		my @tmparray;
  		push(@tmparray,$simplemotif[$ch]);
  		$logchainhash{$chain}=\@tmparray;
  	}		
  	
  }	
  foreach my $thischain (keys %logchainhash){
    my $thisatomchaininfoarrayref=$atomchainhash{$thischain};
    my @thisatomchaininfoarray=@$thisatomchaininfoarrayref;
    my $atomposvsatomindexref=$thisatomchaininfoarray[0];
    my $atomindexvsatomposref=$thisatomchaininfoarray[1];
    my %atomposvsatomindex=%$atomposvsatomindexref;
    my %atomindexvsatompos=%$atomindexvsatomposref;
    
  	my $pname;
  	my $logname;
  	my $nowchain=$thischain;
  	my $pdbidandchain;
  	if($nowchain ne "-"){
  		$pname="P_".$pdbid."_".$nowchain;
  		$pdbidandchain=$pdbid."_".$nowchain;
  		$logname=$pdbid."_".$nowchain;
  	}else{
  		$logname=$pdbid;
  		$pdbidandchain=$pdbid;
  		$pname="P_".$pdbid;
  	}	

  	my $logfile=$resultdir.$logname.".log";	
  	unless(open(INLOG,">>$logfile")){
	 		die("Can not open file $logfile\n");
      exit;
		}
		print INLOG "Get Molecule PDB User $pdbidandchain ".$pname." Heteroatom -Reference_object\n";
		print INLOG "Display Molecule Only Atoms Trace ".$pname."*\n";
  	print INLOG "Center -World $pname Center_of_mass\n";
  	print INLOG "Clip -Auto 250 0 -Perspective INSIGHTII\n";
  	print INLOG "Color Molecule Atoms $pname Specified Specification 192,192,192\n";
  	my $thismotifarrayref=$logchainhash{$thischain};
  	my @thismotifarray=@$thismotifarrayref;
  	for(my $tm=0;$tm<@thismotifarray;$tm++){
  		my $logmotifref=$thismotifarray[$tm];
  		my @logmotifarray=@$logmotifref;
  		for(my $ele=0;$ele<$elenumber;$ele++){
  			my $startpos=$logmotifarray[$ele]->[0];
  			my $endpos=$logmotifarray[$ele]->[1];
  			my @stmp=split(/_/,$startpos);
  			$startpos=$stmp[1];
  			my @etmp=split(/_/,$endpos);
  			$endpos=$etmp[1];
  			my $startindex=$atomposvsatomindex{$startpos};
  			my $startmarkindex=$startindex+1;
  			my $startmark=$atomindexvsatompos{$startmarkindex};
  			
  			if($nowchain ne "-"){
  					$startpos=$nowchain.$startpos;
  					$endpos=$nowchain.$endpos;
  					$startmark=$nowchain.$startmark;
  			}
  			if($ele==0){
  				print INLOG "Color Molecule Atoms ".$pname.":".$startpos.":CA Specified Specification green\n";
  			  print INLOG "Color Molecule Atoms ".$pname.":".$startmark."-".$endpos.":CA Specified Specification 255,0,0\n";
  			}	else{
  				print INLOG "Color Molecule Atoms ".$pname.":".$startpos.":CA Specified Specification 0,255,255\n";
					print INLOG "Color Molecule Atoms ".$pname.":".$startmark."-".$endpos.":CA Specified Specification 255,0,0\n";
    		}		
  			
  		}	
  		
  	}
  	print INLOG "Move ".$pname." Center_of_mass Relative 0,0,0 Screen\n";		
  	 close(INLOG);
  }	
}



print "job finished, thanks you\n";
