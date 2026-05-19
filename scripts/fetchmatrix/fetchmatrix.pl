#!/usr/bin/perl
use strict;
use warnings;
use Math::Trig;

my $numargs=@ARGV;

if($numargs<2){
	print "Aim: get list of SSE,sheet info,interaction matrix and generate handedness of each triple element\n";
	print "Usage: perl getfinalmx.pl [formatedmatrixfile][resultfile]\n";
	exit;
}
my $mxfile=$ARGV[0];
my $resultfile=$ARGV[1];

unless(open(RESULT,">>$resultfile")){
	 	die("Can not open file $resultfile\n");
    exit;
}
unless(open(MATRIXFILE,"$mxfile")){
	 	die("Can not open file $mxfile\n");
    exit;
}
my $M_PI=3.1415926;
my $pdbid;
my @allsheetarray;
my @allinteractionarray;
my @allelearray;
my %chainhash;
while(my $line=<MATRIXFILE>){
	chop($line);
	
	if($line=~/ssd/){
		my @namearray=split(/\./,$line);
		$pdbid=$namearray[0];
		next;
	}
	if($line=~/^$/){
		next;
	}
	if($line=~/sheet/){
		my @tmparray=split(/\s+/,$line);
		my @tmpsheetarray;
		for(my $i=3;$i<@tmparray;$i++){
			push(@tmpsheetarray,$tmparray[$i]);
		}
		push(@allsheetarray,\@tmpsheetarray);
		next;
	}
	if($line=~/\*/){
		my @tmpstrarray=split(//,$line);
		for(my $j=0;$j<@tmpstrarray;$j++){
			if($tmpstrarray[$j] eq " "){   
				$tmpstrarray[$j]=" ";
			}	
		}	
		push(@allinteractionarray,\@tmpstrarray);
		next;
	}else{
			my $eletype=substr($line,0,1);
			my $elechain=substr($line,1,1);
			if($elechain eq " "){
				$elechain="_";
			}	
			if(defined $chainhash{$elechain}){
				
			}else{
				$chainhash{$elechain}=1;
			}		
			my $elepos=substr($line,2,12);
		  my $elelength=substr($line,15,4);
			if($elelength=~/(\d+)/){
				$elelength=$1;
			}	
			my $fx=substr($line,20,8);
			if($fx=~/(-?\d+\.?\d+)/){
				$fx=$1;
			}
			my $fy=substr($line,28,8);
			if($fy=~/(-?\d+\.?\d+)/){
				$fy=$1;
			}
			my $fz=substr($line,36,8);
			if($fz=~/(-?\d+\.?\d+)/){
				$fz=$1;
			}
			my $lx=substr($line,44,8);
			if($lx=~/(-?\d+\.?\d+)/){
				$lx=$1;
			}
			my $ly=substr($line,52,8);
			if($ly=~/(-?\d+\.?\d+)/){
				$ly=$1;
			}
			my $lz=substr($line,60,8);
			if($lz=~/(-?\d+\.?\d+)/){
				$lz=$1;
			}
			my $fcor=$fx.",".$fy.",".$fz;
			my $lcor=$lx.",".$ly.",".$lz;
			my $info=$elechain.":".$elepos.":".$eletype.":".$fcor.":".$lcor;
			push(@allelearray,$info);
			next;
		}					
}	
close(MATRIXFILE);

my $chainsize=keys %chainhash;

my @indexallelearray;
for(my $ia=0;$ia<@allelearray;$ia++){
	my $allindex=$ia+1;
	my $alleleinfo=$allindex.":".$allelearray[$ia];
	push(@indexallelearray,$alleleinfo);
}	

if($chainsize>1){
	print RESULT "there are $chainsize chain in $pdbid\n";
	foreach my $chainkey (sort keys %chainhash){
		print "getting info for $pdbid $chainkey\n";
		print RESULT "the matrix info including SSEs list,sheet and handedness of chain $chainkey:\n";
		print "getting ele list of $pdbid $chainkey\n";
		my $thiselearrayref=getelearrayforthischain(\@allelearray,$chainkey);
		my @thiselearray=@$thiselearrayref;
		print "getting interactionmatrix of $pdbid $chainkey\n";
		my $thisinteractionarrayref=getinteractionarrayforthischain(\@thiselearray,\@allinteractionarray);
		my @thisinteractionarray=@$thisinteractionarrayref;
		print "getting sheetinfo of $pdbid $chainkey\n";
		my $thissheetarrayref=getsheetarrayforthischain(\@allsheetarray,\@thiselearray);
		my @thissheetarray=@$thissheetarrayref;
		print "getting handedness of each triple ele for $pdbid $chainkey\n";
	  my $handednessarrayref=gethandednessarrayforthischain(\@thiselearray);
	  my @handednessarray=@$handednessarrayref;
		for(my $i=0;$i<@thiselearray;$i++){
			my @tmparray=split(/:/,$thiselearray[$i]);
			my $chain=$tmparray[1];
			my $position=$tmparray[2];
			my $type=$tmparray[3];
			if($chain eq "_"){
				$chain=" ";
			}	
			#my $fcor=$tmparray[4];
			#$my $lcor=$tmparray[5];
			#my @ftmparray=split(/\,/,$fcor);
			#my @ltmparray=split(/\,/,$lcor);
			my $outindex=$i+1;
			my $outindexlen=length($outindex);
			if($outindexlen==1){
				$outindex=$outindex."    ";
			}	
			if($outindexlen==2){
				$outindex=$outindex."   ";
			}	
			if($outindexlen==3){
				$outindex=$outindex."  ";
			}
			if($outindexlen==4){
				$outindex=$outindex." ";
			}		
			my $eleinfo=$outindex.$type.$chain." ".$position;
		  print RESULT $eleinfo."\n";	
		}	
		for(my $ti=0;$ti<@thisinteractionarray;$ti++){
			my $tmpinteractionarrayref=$thisinteractionarray[$ti];
			my @tmpinteractionarray=@$tmpinteractionarrayref;
			my $thisinterindex=$ti+1;
			my $thisinterindexlen=length($thisinterindex);
			if($thisinterindexlen==1){
				$thisinterindex=$thisinterindex."   ";
			}	
			if($thisinterindexlen==2){
				$thisinterindex=$thisinterindex."  ";
			}	
			if($thisinterindexlen==3){
				$thisinterindex=$thisinterindex." ";
			}		
			my $thisinteroutinfo=$thisinterindex." ";
			for(my $tj=0;$tj<@tmpinteractionarray;$tj++){
				$thisinteroutinfo=$thisinteroutinfo.$tmpinteractionarray[$tj];
			}
			print RESULT $thisinteroutinfo."\n";
		}	
		
		for(my $j=0;$j<@thissheetarray;$j++){
			my $tmpsheetarrayref=$thissheetarray[$j];
			my @tmpsheetarray=@$tmpsheetarrayref;
			#my $tmpcount=@tmpsheetarray;
			my $sheetindex=$j+1;
			my $sheetindexlen=length($sheetindex);
			if($sheetindexlen==1){
				$sheetindex=$sheetindex."   ";
			}	
			if($sheetindexlen==2){
				$sheetindex=$sheetindex."  ";
			}	
			if($sheetindexlen==3){
				$sheetindex=$sheetindex." ";
			}	
			print RESULT "sheet".$sheetindex;
			
			for(my $jt=0;$jt<@tmpsheetarray;$jt++){
				print RESULT " ".$tmpsheetarray[$jt];
			}	
			print RESULT "\n";
		}	
		for(my $k=0;$k<@handednessarray;$k++){
			print RESULT "handedness\t".$handednessarray[$k]."\n";
		}	
		print RESULT "\n\n";	
	}	
	
	print "getting info of $pdbid allchain\n";
	print "getting handedness info of $pdbid allchain\n";
	my $allhandednessarrayref=gethandednessarrayforthischain(\@indexallelearray);
	my @allhandednessarray=@$allhandednessarrayref;
	print RESULT "the matrix info including SSEs list,sheet and handedness of $pdbid allchain (in case you want to defined matrix involed with different chain\n"; 
	for(my $ai=0;$ai<@allelearray;$ai++){
			my @tmparray=split(/:/,$allelearray[$ai]);
	
			my $chain=$tmparray[0];
			if($chain eq "_"){
				$chain=" ";
			}	
			my $position=$tmparray[1];
			my $type=$tmparray[2];
			my $outindex=$ai+1;
			my $outindexlen=length($outindex);
			if($outindexlen==1){
				$outindex=$outindex."    ";
			}	
			if($outindexlen==2){
				$outindex=$outindex."   ";
			}	
			if($outindexlen==3){
				$outindex=$outindex."  ";
			}
			if($outindexlen==4){
				$outindex=$outindex." ";
			}		
			my $eleinfo=$outindex.$type.$chain." ".$position;
		  print RESULT $eleinfo."\n";	
	}
	for(my $an=0;$an<@allinteractionarray;$an++){
		my $tmpinterarrayref=$allinteractionarray[$an];
		my @tmpinterarray=@$tmpinterarrayref;
		my $outinterindex=$an+1;
		my $outinterindexlen=length($outinterindex);
		if($outinterindexlen==1){
			$outinterindex=$outinterindex."   ";
		}	
		if($outinterindexlen==2){
			$outinterindex=$outinterindex."  ";
		}	
		if($outinterindexlen==3){
			$outinterindex=$outinterindex." ";
		}	
		my $outinterinfo=$outinterindex." ";
		for(my $ani=0;$ani<@tmpinterarray;$ani++){
			$outinterinfo=$outinterinfo.$tmpinterarray[$ani];
		}	
		print RESULT $outinterinfo."\n";
	}	
	for(my $aj=0;$aj<@allsheetarray;$aj++){
			my $tmpsheetarrayref=$allsheetarray[$aj];
			my @tmpsheetarray=@$tmpsheetarrayref;
			my $sheetindex=$aj+1;
			my $sheetindexlen=length($sheetindex);
			if($sheetindexlen==1){
				$sheetindex=$sheetindex."   ";
			}	
			if($sheetindexlen==2){
				$sheetindex=$sheetindex."  ";
			}	
			if($sheetindexlen==3){
				$sheetindex=$sheetindex." ";
			}	
			print RESULT "sheet".$sheetindex;
			for(my $ajt=0;$ajt<@tmpsheetarray;$ajt++){
				print RESULT " ".$tmpsheetarray[$ajt];
			}	
			print RESULT "\n";
	}
	for (my $ak=0;$ak<@allhandednessarray;$ak++){
		print RESULT "handedness\t".$allhandednessarray[$ak]."\n";
		
	}	
		
}
elsif($chainsize==1){
	print "getting info of $pdbid (only one chain)\n";
	print "getting handedness info of $pdbid (only one chain)\n";
	my $allhandednessarrayref=gethandednessarrayforthischain(\@indexallelearray);
	my @allhandednessarray=@$allhandednessarrayref;
	print RESULT "the matrix info including SSEs list,sheet and handedness of $pdbid (only one chain)\n"; 
	for(my $ai=0;$ai<@allelearray;$ai++){
			my @tmparray=split(/:/,$allelearray[$ai]);
	
			my $chain=$tmparray[0];
			if($chain eq "_"){
				$chain=" ";
			}	
			my $position=$tmparray[1];
			my $type=$tmparray[2];
			my $outindex=$ai+1;
			my $outindexlen=length($outindex);
			if($outindexlen==1){
				$outindex=$outindex."    ";
			}	
			if($outindexlen==2){
				$outindex=$outindex."   ";
			}	
			if($outindexlen==3){
				$outindex=$outindex."  ";
			}
			if($outindexlen==4){
				$outindex=$outindex." ";
			}		
			my $eleinfo=$outindex.$type.$chain." ".$position;
		  print RESULT $eleinfo."\n";	
	}
	for(my $an=0;$an<@allinteractionarray;$an++){
		my $tmpinterarrayref=$allinteractionarray[$an];
		my @tmpinterarray=@$tmpinterarrayref;
		my $outinterindex=$an+1;
		my $outinterindexlen=length($outinterindex);
		if($outinterindexlen==1){
			$outinterindex=$outinterindex."   ";
		}	
		if($outinterindexlen==2){
			$outinterindex=$outinterindex."  ";
		}	
		if($outinterindexlen==3){
			$outinterindex=$outinterindex." ";
		}	
		my $outinterinfo=$outinterindex." ";
		for(my $ani=0;$ani<@tmpinterarray;$ani++){
			$outinterinfo=$outinterinfo.$tmpinterarray[$ani];
		}	
		print RESULT $outinterinfo."\n";
	}	
	for(my $aj=0;$aj<@allsheetarray;$aj++){
			my $tmpsheetarrayref=$allsheetarray[$aj];
			my @tmpsheetarray=@$tmpsheetarrayref;
			my $sheetindex=$aj+1;
			my $sheetindexlen=length($sheetindex);
			if($sheetindexlen==1){
				$sheetindex=$sheetindex."   ";
			}	
			if($sheetindexlen==2){
				$sheetindex=$sheetindex."  ";
			}	
			if($sheetindexlen==3){
				$sheetindex=$sheetindex." ";
			}	
			print RESULT "sheet".$sheetindex;
			for(my $ajt=0;$ajt<@tmpsheetarray;$ajt++){
				print RESULT " ".$tmpsheetarray[$ajt];
			}	
			print RESULT "\n";
	}
	for (my $ak=0;$ak<@allhandednessarray;$ak++){
		print RESULT "handedness\t".$allhandednessarray[$ak]."\n";
		
	}	
}else{
	print "error for parsing matrixfile\n";
	exit;
}		
print "all job finished\n";
close(RESULT);

sub getelearrayforthischain{
	my $elearrayref=shift(@_);
	my @elearray=@$elearrayref;
	my $thischain=shift(@_);
	
	my @thiselearray;
	for(my $i=0;$i<@elearray;$i++){
		my @tmparray=split(/:/,$elearray[$i]);
		my $tmpchain=$tmparray[0];
		if($tmpchain eq $thischain){
			my $tmpindex=$i+1;
			my $info=$tmpindex.":".$elearray[$i];
			push(@thiselearray,$info);
		}	
	}
	return(\@thiselearray);
}	
sub getsheetarrayforthischain{
	my $sheetarrayref=shift(@_);
	my @sheetarray=@$sheetarrayref;
	my $thiselearrayref=shift(@_);
	my @thiselearray=@$thiselearrayref;
  
  my @thissheetarray;
  my %thiselehash;
  for(my $i=0;$i<@thiselearray;$i++){
  	my @tmparray=split(/:/,$thiselearray[$i]);
  	my $thisindex=$tmparray[0];
  	$thiselehash{$thisindex}=$i+1;
  }	
	
	for(my $j=0;$j<@sheetarray;$j++){
		my $tmpsheetref=$sheetarray[$j];
		my @tmpsheet=@$tmpsheetref;
		my @outsheet;
		for(my $k=0;$k<@tmpsheet;$k++){
			if(defined $thiselehash{$tmpsheet[$k]}){
				push(@outsheet,$thiselehash{$tmpsheet[$k]});
			}	
		}	
		my $outsheetcount=@outsheet;
		if($outsheetcount==0){
			
		}else{
			push(@thissheetarray,\@outsheet);
		}		
	}	
	return(\@thissheetarray);
}	

sub gethandednessarrayforthischain{
	my $elearray=shift(@_);
	my @elearray=@$elearray;
	my @handednessarray;
	
	for(my $i=0;$i<@elearray;$i++){
		for(my $k=0;$k<@elearray;$k++){
			  if($k==$i){
			  	next;
			  }	
				for(my $j=0;$j<@elearray;$j++){
					if(($j==$k) or ($j==$i)){
						next;
					}	
					my $segA=$elearray[$i];
					my ($segAfx,$segAfy,$segAfz,$segAlx,$segAly,$segAlz)=getcor($segA);
					
				  my $segB=$elearray[$k];
				  my ($segBfx,$segBfy,$segBfz,$segBlx,$segBly,$segBlz)=getcor($segB);
				  
					my $segC=$elearray[$j];
					my ($segCfx,$segCfy,$segCfz,$segClx,$segCly,$segClz)=getcor($segC);
					
					
					my $N1x; my $N1y; my$N1z;
    			my $N2x;
    			my $N2y;
    			my $N2z;
    			my $N3x;
    			my $N3y;
    			my $N3z;
    			my $Nx;
    			my $Ny;
    			my $Nz;
    			my $Mx;
    			my $My;
    			my $Mz;
    			my $ang;
    			my $m1x;
    			my $m1y;
    			my $m1z;
    			my $m2x;
    			my $m2y;
    			my $m2z;
    			my $m3x;
    			my $m3y;
    			my $m3z;
    			$N1x = $segAlx-$segAfx;
    			$N1y = $segAly-$segAfy;
    			$N1z = $segAlz-$segAfz;

          $ang= angle($segAfx,$segAfy,$segAfz,$segAlx,$segAly,$segAlz,$segBfx,$segBfy,$segBfz,$segBlx,$segBly,$segBlz);     
          if ($ang<90.001 ){
                $N2x = $segBlx-$segBfx;
                $N2y = $segBly-$segBfy;
                $N2z = $segBlz-$segBfz;

        	}else{
                $N2x = $segBfx-$segBlx;
                $N2y = $segBfy-$segBly;
                $N2z = $segBfz-$segBlz;
        	}
        	$ang = angle($segAfx,$segAfy,$segAfz,$segAlx,$segAly,$segAlz,$segCfx,$segCfy,$segCfz,$segClx,$segCly,$segClz);
        	if ($ang< 90.001){
                $N3x = $segClx-$segCfx;
                $N3y = $segCly-$segCfy;
                $N3z = $segClz-$segCfz;
        	}else{
                $N3x = $segCfx-$segClx;
                $N3y = $segCfy-$segCly;
                $N3z = $segCfz-$segClz;
        	}

        	$Nx = ( $N1x + $N2x + $N3x ) / 3;
        	$Ny = ( $N1y + $N2y + $N3y ) / 3;
        	$Nz = ( $N1z + $N2z + $N3z ) / 3;

			    $m1x = ($segAlx+$segAfx)/2;
    			$m1y = ($segAly+$segAfy)/2;
    			$m1z = ($segAlz+$segAfz)/2;
    			$m2x = ($segBlx+$segBfx)/2;
    			$m2y = ($segBly+$segBfy)/2;
    			$m2z = ($segBlz+$segBfz)/2;
    			$m3x = ($segClx+$segCfx)/2;
    			$m3y = ($segCly+$segCfy)/2;
    			$m3z = ($segClz+$segCfz)/2;
    			$Mx = ($m2y-$m1y)*($m3z-$m1z)-($m2z-$m1z)*($m3y-$m1y);
    			$My = ($m2z-$m1z)*($m3x-$m1x)-($m2x-$m1x)*($m3z-$m1z);
    			$Mz = ($m2x-$m1x)*($m3y-$m1y)-($m2y-$m1y)*($m3x-$m1x);
    			$ang = angle(0,0,0,$Nx,$Ny,$Nz,0,0,0,$Mx,$My,$Mz);
    		   
    		   my $indexi=$i+1;
    		   my $indexj=$j+1;
    		   my $indexk=$k+1; 
    			if($ang<90){
       			my $handedness=$indexi." ".$indexk." ".$indexj." "."R";
       			push(@handednessarray,$handedness);
    			}  
    			elsif($ang>90.001){
    				my $handedness=$indexi." ".$indexk." ".$indexj." "."L";
       			push(@handednessarray,$handedness);
       			
    			}   
    			else{
    	 			my $handedness=$indexi." ".$indexk." ".$indexj." "."N";
       			push(@handednessarray,$handedness);
    			}   
			}	
		}	
	}	
	return(\@handednessarray);
}	

sub getcor{
	my $eleinfo=shift(@_);
	my @tmparray=split(/:/,$eleinfo);
	my $fcor=$tmparray[4];
	my $lcor=$tmparray[5];
	
	my @tmpfarray=split(/\,/,$fcor);
	my $fcorx=$tmpfarray[0];
	my $fcory=$tmpfarray[1];
	my $fcorz=$tmpfarray[2];
	
	my @tmplarray=split(/\,/,$lcor);
	my $lcorx=$tmplarray[0];
	my $lcory=$tmplarray[1];
	my $lcorz=$tmplarray[2];
	
	return($fcorx,$fcory,$fcorz,$lcorx,$lcory,$lcorz);
}	
sub angle{

	my $mix=shift(@_);
	my $miy=shift(@_);
	my $miz=shift(@_);
	my $ix=shift(@_);
	my $iy=shift(@_);
	my $iz=shift(@_);
	my $mjx=shift(@_);
	my $mjy=shift(@_);
	my $mjz=shift(@_);
	my $jx=shift(@_);
	my $jy=shift(@_);
	my $jz=shift(@_);
	my $dot;
	my $ma;
	my $mb;
  my $cos_theta;
  my $theta;
  $dot = ( $ix - $mix ) * ( $jx - $mjx ) + ( $iy - $miy ) * ( $jy - $mjy ) 
                + ( $iz - $miz ) * ( $jz - $mjz );

        $ma = measure_distance ( $mix, $miy, $miz, $ix, $iy, $iz );
        $mb = measure_distance ( $mjx, $mjy, $mjz, $jx, $jy, $jz );

        $cos_theta = $dot / ( $ma * $mb ); 
        if ( ( $cos_theta > 1.01 )|| ( $cos_theta < -1.01 ) )
        {       
          
        
        }
        if ( $cos_theta > 1 ){
        	$cos_theta = 1;
        }	
        if ( $cos_theta < -1 ){
                $cos_theta = -1;
        }        
        $theta = abs ( acos ( $cos_theta ) / $M_PI * 180 );
        return $theta;
}	
sub measure_distance{
	my $ix=shift(@_);
	my $iy=shift(@_);
	my $iz=shift(@_);
	my $jx=shift(@_);
	my $jy=shift(@_);
	my $jz=shift(@_);
	my $dis= sqrt ((($jx-$ix)*($jx-$ix))+(($jy-$iy)*($jy-$iy))+(($jz-$iz)*($jz-$iz)));
	return $dis;
}

sub getinteractionarrayforthischain{
	my $elearrayref=shift(@_);
	my @elearray=@$elearrayref;
	my $interactionarrayref=shift(@_);
	my @interactionarray=@$interactionarrayref;
	
	my %elelist;
	for(my $j=0;$j<@elearray;$j++){
		my $noweleinfo=$elearray[$j];
		my @tmparray=split(/:/,$noweleinfo);
		my $indexpos=$tmparray[0]-1;
		$elelist{$indexpos}=1;
	}	
	
	my @thisinteractionarray;
	for(my $i=0;$i<@interactionarray;$i++){
		if(defined $elelist{$i}){
			my $tmpinteractionarrayref=$interactionarray[$i];
			my @tmpinteractionarray=@$tmpinteractionarrayref;
			my @tmpoutinteractionarray;
			for(my $j=0;$j<@tmpinteractionarray;$j++){
				if(defined $elelist{$j}){
					push(@tmpoutinteractionarray,$tmpinteractionarray[$j]);
				}
			}	
		  push(@thisinteractionarray,\@tmpoutinteractionarray);
		}
	}	
	return(\@thisinteractionarray);
	
}	
