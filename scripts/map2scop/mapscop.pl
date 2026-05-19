#!/usr/bin/perl
use strict;
use warnings;
use DBI;
my $driver="DBI:mysql:SCOP:localhost"; # or put your MySQL database info for SCOP databasel
my $user="yourusrname";  #put your usrname;
my $password="yourpassword"; #put your password;
my $dbh=DBI->connect ($driver,$user,$password,{RaiseError=>1});

my $numargs = @ARGV;


if($numargs<5){
	
	print "Aim: Map pdbhits generated from ProSMoS program to SCOP on superfamily and fold level;\n";
	print "Usage: [pdbhitslist][pdbhitsdir][superfamilydir][folddir][elementnumber]\n";
	print "For example:\n";
	print "perl mapscop.pl .\\pdbhitlist .\\pdbhits\\ .\\superfam\\ .\\fold\\ 5"."\n";
	exit;
}
my $pdbhitlistfile;
my $pdbhitdir;
my $superfmdir;
my $folddir;
my $elementnumber;
my $scopdomaininfofile;

	$pdbhitlistfile=$ARGV[0];
	$pdbhitdir=$ARGV[1]; 	
  $superfmdir=$ARGV[2];
	$folddir=$ARGV[3];
	$elementnumber=$ARGV[4];
if(-e "$pdbhitlistfile"){
	
}else{
	print "error can not found pdbhitslist\n";
	exit;
}		
if(-e "$pdbhitdir"){
	
}else{
	print "error can not found directory of pdbhits";
	exit;
}	
if(-e "$superfmdir"){
	
}else{
	print "please mkdir for $superfmdir"."\n";
}	
if(-e "$folddir"){
	
}else{
	print "please mkdir for $folddir"."\n";
}		
	
my $notfoundfile=$superfmdir."pdbhitsnorecordinscop";
my $notfoundfile2=$superfmdir."pdbhitsfailmaptoscop";
my $foundfile=$superfmdir."pdbdomainmapped";
unless(open(HIT,"$pdbhitlistfile")){
	 die("Can not open file $pdbhitlistfile\n");
      exit;
}
unless(open(NOTFOUND,">>$notfoundfile")){
	 die("Can not open file $notfoundfile\n");
      exit;
}
unless(open(NOTFOUNDTWO,">>$notfoundfile2")){
	 die("Can not open file $notfoundfile\n");
      exit;
}
unless(open(FOUND,">>$foundfile")){
	 die("Can not open file $foundfile\n");
      exit;
}
my $errorlogfile=$superfmdir."domainmapcheckinfo";
unless(open(ERROR,">>$errorlogfile")){
	 die("Can not open file $errorlogfile\n");
      exit;
}


my @alldomainid;
while(my $line=<HIT>){
	chop($line);
	my $pdbid=substr($line,3,4);
	my $pdbsth=$dbh->prepare("select nodeid,pdbchain from dirdes171 where pdbid='$pdbid'");
	my $pdbcount=0;
	$pdbcount=$pdbsth->execute();
	if($pdbcount==0){
		print NOTFOUND $pdbid."\n";
		print $pdbid."\tnot found in scop\n";  
		next;     
		$pdbsth->finish();                                            
	}else{
		print "mapping $pdbid to scop domain....\n";
		my $pdbfile=$pdbhitdir.$line;
		unless(open(PDBFILE,"$pdbfile")){
	 		die("Can not open file $pdbfile\n");
      exit;
		}
		my @matrixarray;
		while(my $matrixline=<PDBFILE>){
			chop($matrixline);
			if(($matrixline=~/pdb/) or ($matrixline=~/sub-matrix/) or ($matrixline=~/MOTIF:/) or ($matrixline=~/^END$/)){
				
			}else{
				push(@matrixarray,$matrixline);
			}		
		}	
		my @pdbmotif;
		
		for(my $putflag=0;$putflag<@matrixarray;$putflag=$putflag+$elementnumber){
			my @pdbmatrix;
			my $alllength=0;
			for(my $i=$putflag;$i<$putflag+$elementnumber;$i++){
				my @tempinfo=();
			
				my $startpos=substr($matrixarray[$i],36,5);
		    my $endpos=substr($matrixarray[$i],43,5);
				my $chaininfo=substr($matrixarray[$i],49,1);
				if($chaininfo eq " "){
					$chaininfo="-";
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
						print "error".$pdbid."\n";
						exit;
				}		
			  my $elelength=substr($matrixarray[$i],58,4);
			  if($elelength=~/(\d+)/){
			  	$elelength=$1;
			  }	
		    $alllength=$elelength+$alllength;
		  
		 		push(@tempinfo,$startpos);
				push(@tempinfo,$endpos);
				push(@tempinfo,$elelength);
				push(@tempinfo,$chaininfo);
				push(@pdbmatrix,\@tempinfo);
			}
			push(@pdbmatrix,$alllength);
			my $temphitflag="not find";
			push(@pdbmatrix,$temphitflag);
			push(@pdbmotif,\@pdbmatrix);
	  }
		close(PDBFILE);
		my @pdbdomain=();
	
		while(my ($nodeid,$pdbchain)=$pdbsth->fetchrow_array()){
				my @tempdomain=();
			  push(@tempdomain,$nodeid);
				push(@tempdomain,$pdbchain);
				my %thisscopdomainrange=();
				my $domainsth=$dbh->prepare("select domainatom from pdbstyle171 where nodeid=$nodeid");
				my $existflag=$domainsth->execute();		   
				
				if($existflag==1){
					my ($domainatom)=$domainsth->fetchrow_array();
					my %thisscopdomainrange=();
					if($domainatom eq "noinfo"){
						$thisscopdomainrange{"length"}=0;
						push(@tempdomain,\%thisscopdomainrange);
						push(@pdbdomain,\@tempdomain);		
					}else{	
						my @domainatomarray=split(/\|\|/,$domainatom);
						my %domainhash;
						my $atomcount=@domainatomarray;
						$atomcount=$atomcount-1;
						for(my $dindex=1;$dindex<@domainatomarray;$dindex++){
							$domainhash{$domainatomarray[$dindex]}=$dindex;
						}	
						$domainhash{"length"}=$atomcount;
				
						%thisscopdomainrange=%domainhash;
						push(@tempdomain,\%thisscopdomainrange);		
						push(@pdbdomain,\@tempdomain);		  
					}
				}else{
				  $thisscopdomainrange{"length"}=0;
					push(@tempdomain,\%thisscopdomainrange);
					push(@pdbdomain,\@tempdomain);			   
					print ERROR "error: can not find $pdbid $nodeid in pdbstyle171\n";
					
				}
				$domainsth->finish();
				
		}	
		
		$pdbsth->finish();
		my %alldomainrangeck=();
		for(my $tdomaincount=0;$tdomaincount<@pdbdomain;$tdomaincount++){
			my $tnowdomainrangeref=$pdbdomain[$tdomaincount]->[2];
			my %tnowdomainrange=%$tnowdomainrangeref;
			foreach my $nowkey (keys %tnowdomainrange){
				$alldomainrangeck{$nowkey}=1;
			}	
		}	
		for(my $domaincount=0;$domaincount<@pdbdomain;$domaincount++){
			my $domainid=$pdbdomain[$domaincount]->[0];
			my $domainchain=$pdbdomain[$domaincount]->[1]; #pdbchain also a range;
			my $nowdomainrangeref=$pdbdomain[$domaincount]->[2];
			my %nowdomainrange=%$nowdomainrangeref;
		  my $domainlength=$nowdomainrange{"length"};
		  my %thisscopdomainterminalinfo=();
		  if($domainchain eq "-"){
		  	
		  }elsif($domainchain=~/^(\w):$/){
		  	
		  }else{
				my @splitdomain=split(/\,/,$domainchain);
				for(my $splitcount=0;$splitcount<@splitdomain;$splitcount++){
					if($splitdomain[$splitcount]=~/^(-?\d+\w?)-(\d+\w?)$/){
						my $nowstart=$1;
						my $nowend=$2;
						$nowstart="-"."_".$nowstart;
						$nowend="-"."_".$nowend;
		
						my $nowstartindex=$nowdomainrange{$nowstart};
						my $nowendindex=$nowdomainrange{$nowend};
						if(defined $nowstartindex){
						
						}else{
							
							print $domainid."\t".$nowstart;
							exit;
						}		
						if(defined $nowendindex){
							
						}else{
							print $domainid."\t".$nowend;
							exit;
						}		
						$thisscopdomainterminalinfo{$nowstartindex}=1;
						$thisscopdomainterminalinfo{$nowendindex}=1;
					}elsif($splitdomain[$splitcount]=~/^(\w+):(-?\d+\w?)-(\d+\w?)$/){
						my $nowchain=$1;
						my $nowstart=$2;
						my $nowend=$3;
						my $nowupstart=$nowchain."_".$nowstart;
						my $nowupend=$nowchain."_".$nowend;
						
						my $nowstartindex=$nowdomainrange{$nowupstart};
						if(defined $nowstartindex){
							
						
						}else{
							my $lcchain=lc($nowchain);
							my $nowlcstart=$lcchain."_".$nowstart;
					
							$nowstartindex=$nowdomainrange{$nowlcstart};
							if(defined $nowstartindex){
								
							}else{
								print $domainid."\t".$nowlcstart;
								exit;	
							}		
							
						}		
					
						my $nowendindex=$nowdomainrange{$nowupend};
						if(defined $nowendindex){
							
						}else{
							my $lcchain=lc($nowchain);
							my $nowlcend=$lcchain."_".$nowend;
						
							$nowendindex=$nowdomainrange{$nowlcend};
							if(defined $nowendindex){
								
							}else{
								print $domainid."\t".$nowlcend;
								exit;	
							}		
						}	
						
						$thisscopdomainterminalinfo{$nowstartindex}=1;
						$thisscopdomainterminalinfo{$nowendindex}=1;
					}elsif($splitdomain[$splitcount]=~/^(\w):$/){
						my $nowstartindex=1;
						my $nowendindex=keys %nowdomainrange;
						$thisscopdomainterminalinfo{$nowstartindex}=1;
						$thisscopdomainterminalinfo{$nowendindex}=1;
				  }else{
						print "$splitdomain[$splitcount]\t$domainid\n";
						print "unparsed situation, program bugs, please contact me\n";
						exit;
					}			
				}	
			}
		  for(my $motifcount=0;$motifcount<@pdbmotif;$motifcount++){
				my $findlength=0;
				my $motifwholelength=$pdbmotif[$motifcount][$elementnumber];
				my $thismotifhitflag=$pdbmotif[$motifcount][$elementnumber+1];
				if($thismotifhitflag ne "not find"){
					next;
				}	
				for(my $elementcount=0;$elementcount<$elementnumber;$elementcount++){
					my $elementstart=$pdbmotif[$motifcount][$elementcount][0];
					my $elementend=$pdbmotif[$motifcount][$elementcount][1];
					my $elementlength=$pdbmotif[$motifcount][$elementcount][2];
					my $elementchain=$pdbmotif[$motifcount][$elementcount][3];
					if($domainlength==0){
						if($domainchain eq "-"){
							$findlength=$findlength+$elementlength;
						}elsif($domainchain=~/^(\w):$/){
							my $tempchain=$1;
							my $templcchain=lc($tempchain);
							if($tempchain eq $elementchain){
								$findlength=$findlength+$elementlength;
							}elsif($templcchain eq $elementchain){
								$findlength=$findlength+$elementlength;
							}else{
								
							}	
					  }else{
					  	print ERROR "altought domainlength==0 but it has split range (not whole chain) $domainid\n";
		
					  }
					  next;		
					}else{
						if($domainchain eq "-"){
							$findlength=$findlength+$elementlength;
						}elsif($domainchain=~/^(\w):$/){
							my $tempchain=$1;
							my $templcchain=lc($tempchain);
							if($tempchain eq $elementchain){
								$findlength=$findlength+$elementlength;
							}elsif($templcchain eq $elementchain){
								$findlength=$findlength+$elementlength;
							}else{
								
							}		
						}	
						elsif(defined $nowdomainrange{$elementstart} and defined $nowdomainrange{$elementend}){				
								$findlength=$findlength+$elementlength;
						}elsif(defined $nowdomainrange{$elementstart}){
							  my $tmpsize=keys %thisscopdomainterminalinfo;
							  if($tmpsize==0){
							  	print "ERROR can not find scopdomainterminal info for $domainid , please contact me\n";
							  	exit;
							  }
							  my $elementstartindomainindex=$nowdomainrange{$elementstart};
							  my @termarray;
							  foreach my $termkey (keys %thisscopdomainterminalinfo){
									push(@termarray,$termkey);
								}
								my @newtermarray=sort {$a <=> $b} @termarray;
								for(my $ti=0;$ti<@newtermarray-1;$ti++){
									if(($newtermarray[$ti]<=$elementstartindomainindex) and ($elementstartindomainindex<=$newtermarray[$ti+1])){
										$findlength=$findlength+$newtermarray[$ti+1]-$elementstartindomainindex+1;
										last;
									}	
								}	
								if(defined $alldomainrangeck{$elementend}){
									
								}else{
									print ERROR "can not find end\t".$pdbid."\t".$elementend."\n";
								}		
						}elsif(defined $nowdomainrange{$elementend}){
							  my $tmpsize=keys %thisscopdomainterminalinfo;
							  if($tmpsize==0){
							  	print "ERROR can not find scopdomainterminal info for $domainid , please contact me\n";
							  	exit;
							  }
							  my $elementendindomainindex=$nowdomainrange{$elementend};
							  my @termarray;
							  foreach my $termkey (keys %thisscopdomainterminalinfo){
									push(@termarray,$termkey);
								}
								my @newtermarray=sort {$a <=> $b} @termarray;
								for(my $ti=0;$ti<@newtermarray-1;$ti++){
									if(($newtermarray[$ti]<=$elementendindomainindex) and ($elementendindomainindex<=$newtermarray[$ti+1])){
										$findlength=$findlength+$elementendindomainindex-$newtermarray[$ti]+1;
										last;
									}	
								}	
							
									
								if(defined $alldomainrangeck{$elementstart}){
									
								}else{
									print ERROR "can not find start\t".$pdbid."\t".$elementstart;
								}		
								
						}else{
								
								if(defined $alldomainrangeck{$elementstart}){
									
								}else{
									print ERROR "can not find start\t".$pdbid."\t".$elementstart."\n";
								}	
								if(defined $alldomainrangeck{$elementend}){
									
								}else{
									print ERROR "can not find end\t".$pdbid."\t".$elementend."\n";
								}			
								
						}
					}	
				}	
				if($findlength>$motifwholelength/2){
					$pdbmotif[$motifcount][$elementnumber+1]=$domainid;
				}	
			}	
	  }
	  for(my $cmotifcount=0;$cmotifcount<@pdbmotif;$cmotifcount++){
	  	my $nowfindflag=$pdbmotif[$cmotifcount][$elementnumber+1];
	  	if($nowfindflag ne "not find"){
	  		print FOUND $pdbid."\t".$nowfindflag;
	  		for(my $celementcount=0;$celementcount<$elementnumber;$celementcount++){
					my $celementstart=$pdbmotif[$cmotifcount][$celementcount][0];
					my $celementend=$pdbmotif[$cmotifcount][$celementcount][1];
					my $celementlength=$pdbmotif[$cmotifcount][$celementcount][2];
					my $celementchain=$pdbmotif[$cmotifcount][$celementcount][3];
					print FOUND "\t"."chain:$celementchain ".$celementstart."--".$celementend." Length:".$celementlength;
	  		}
	  		print FOUND "\n";
	  		push(@alldomainid,$nowfindflag);
	  	}	else{
	  		print NOTFOUNDTWO $pdbid."\n";
	  	}	
	  }	
	}
}


close(FOUND);
close(NOTFOUND);
close(ERROR);
#%allnodehash=();
my %foldhash;
my %superfamilyhash;
my %domainhash;

print "mapping to superfamily and fold level.....\n";
for(my $domainnumber=0;$domainnumber<@alldomainid;$domainnumber++){
	my $nowdomainid=$alldomainid[$domainnumber];
	my $domainsth=$dbh->prepare("select sid,pdbid,sccsid,description from dircla171 where nodeid=$nowdomainid");
	my $count=$domainsth->execute();
	my $foldid;
	my $superfamilyid;
	my $domainsid;
	my $pdbid;
	my $pdbchain;
  my $domainsccsid;
	if($count==1){
		my ($nowsid,$nowpdbid,$nowsccsid,$nowdescription)=$domainsth->fetchrow_array();
		my @category=split(/,/,$nowdescription);
		$domainsid=$nowsid;
		$pdbid=$nowpdbid;
		$domainsccsid=$nowsccsid;

		for(my $i=0;$i<@category;$i++){
			if($category[$i]=~/cf=(\d+)/){
				$foldid=$1;
				
			}elsif($category[$i]=~/sf=(\d+)/){
				$superfamilyid=$1;
			}else{
				
				
			}				
		}
		$domainsth->finish();
		my $domaindessth=$dbh->prepare("select pdbchain from dirdes171 where nodeid=$nowdomainid");
		my $descount=$domaindessth->execute();
		if($descount==1){
			my ($nowpdbchain)=$domaindessth->fetchrow_array();
			$pdbchain=$nowpdbchain;
			
			$domaindessth->finish();
		}else{
			print "des error\t".$nowdomainid."\n";
			$domaindessth->finish();
			exit;
		}	
			
	}else{
		print "error\t".$nowdomainid;
		$domainsth->finish();
		exit;
	}
	my @domaininfo;

	push(@domaininfo,$domainsccsid);
	push(@domaininfo,$domainsid);
	push(@domaininfo,$pdbid);
	push(@domaininfo,$pdbchain); 		
  $domainhash{$nowdomainid}=\@domaininfo;
  if($superfamilyhash{$superfamilyid}){
  	my $temp=$superfamilyhash{$superfamilyid};
  	my @temparray=@$temp;
  	push(@temparray,$nowdomainid);
  	$superfamilyhash{$superfamilyid}=\@temparray;
  	
  }else{
  	my @temparray=();
  	push(@temparray,$nowdomainid);
  	$superfamilyhash{$superfamilyid}=\@temparray;
  }		
  if($foldhash{$foldid}){
  	my $temp=$foldhash{$foldid};
  	my @temparray=@$temp;
  	my $flag=0;
  	for(my $k=0;$k<@temparray;$k++){
  		if($temparray[$k]==$superfamilyid){
  			$flag=1;
  			last;
  		}	
  	}	
  	if($flag==0){
  		push(@temparray,$superfamilyid);
  		$foldhash{$foldid}=\@temparray;
  	}
  	
  }else{
  	my @temparray=();
  	push(@temparray,$superfamilyid);
  	$foldhash{$foldid}=\@temparray;
  }		
	
}	
my $superfmlistpath;
$superfmlistpath=$superfmdir."superfamilyfilelist";
unless(open(SUPERFAMFILE,">>$superfmlistpath")){
	 die("Can not open file $superfmlistpath\n");
      exit;
}   
my $superfmallpath;
$superfmallpath=$superfmdir."superfamilyallinfo";
unless(open(SUPERFAMALL,">>$superfmallpath")){
	 die("Can not open file $superfmallpath\n");
      exit;
}   

foreach my $superfamilyout (sort keys(%superfamilyhash)){
	
	my $superdes="";
	my $supercom="";
	my $superdessth=$dbh->prepare("select description from dirdes171 where nodeid=$superfamilyout");
	my $descount=$superdessth->execute();
	if($descount==1){
		$superdes=$superdessth->fetchrow_array();
		
		chop($superdes);
	}	
	$superdessth->finish();
	
	my $supercomsth=$dbh->prepare("select description from dircom171 where nodeid=$superfamilyout");
	my $comcount=$supercomsth->execute();
	if($comcount==1){
		$supercom=$supercomsth->fetchrow_array();
		chop($supercom);
	}	
	$supercomsth->finish();
	
	
  my $name=$superfamilyout;
	print SUPERFAMFILE $superfamilyout."\t".$superdes."\t".$supercom."\n";
	my $superfmsinglepath=$superfmdir.$name;
	
	unless(open(SUPERFAM,">>$superfmsinglepath")){
	 die("Can not open file $superfmsinglepath \n");
      exit;
  }    
  print SUPERFAM "superfamily::\t".$superfamilyout."\t".$superdes."\t".$supercom."\n";
  print SUPERFAMALL "superfamily::\t".$superfamilyout."\t".$superdes."\t".$supercom."\n";
  my $alldomainid=$superfamilyhash{$superfamilyout}; 
  my @alldomainarray=@$alldomainid;
  for(my $k=0;$k<@alldomainarray;$k++){
  	if($domainhash{$alldomainarray[$k]}){
  		my $domaininfo=$domainhash{$alldomainarray[$k]};
  		my @specdomainarray=@$domaininfo;
  		print SUPERFAM $alldomainarray[$k];
  		for(my $j=0;$j<@specdomainarray;$j++){
  			print SUPERFAM "\t".$specdomainarray[$j];
  			print SUPERFAMALL "\t".$specdomainarray[$j];
  		}	
  		print SUPERFAM "\n";
  		print SUPERFAMALL "\n";
  	}else{
  		
  		print "error".$alldomainarray[$k]."\n";
  		exit;
  	}		
  }
   close(SUPERFAM);   
}
close(SUPERFAMALL);
close(SUPERFAMFILE);
my $foldlistpath;
$foldlistpath=$folddir."foldlist";
unless(open(FOLDFILE,">>$foldlistpath")){
	 die("Can not open file $foldlistpath\n");
      exit;
}   
foreach my $foldout (sort keys(%foldhash)){
	my $allsupfamid=$foldhash{$foldout};
	my @allsupfamarray=@$allsupfamid;
	my $folddessth=$dbh->prepare("select sccsid,description from dirdes171 where nodeid=$foldout");
	my $folddescount=$folddessth->execute();
	my $foldsccsid="";
	my $folddes="";
	my $foldcom="";
	if($folddescount==1){
		my ($nowsccsid,$nowdes)=$folddessth->fetchrow_array();
		$foldsccsid=$nowsccsid;
		$folddes=$nowdes;
		chop($folddes);
	}else{
		print "error in fold des\t".$foldout."\n";
	}		
	$folddessth->finish();
	my $foldcomsth=$dbh->prepare("select description from dircom171 where nodeid=$foldout");
	my $foldcomcount=$foldcomsth->execute();
	if($foldcomcount==1){
		my ($desnow)=$foldcomsth->fetchrow_array();
		$foldcom=$desnow;
		chop($foldcom);
	}else{
		
	}		
	$foldcomsth->finish();
	print FOLDFILE $foldout."\t".$foldsccsid."\t".$folddes."\t".$foldcom;
	for(my $n=0;$n<@allsupfamarray;$n++){
		print FOLDFILE "\t".$allsupfamarray[$n];
		
	}	
	print FOLDFILE "\n";
	
	
}	
close(FOLDFILE);
print "Congratulations, Mission Complete\n";
exit;
	
