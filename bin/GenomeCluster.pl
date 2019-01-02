#!/usr/bin/env perl

use warnings;

###########################################################################
###########################################################################

#              ********************************
#              ** GenomeCluster(version 1.0) **
#              **         stand-alone        **
#              ********************************
#
#       Computational Genomics and Bioinformatics Group
#          University of Granada, Dept. of Genetics
#              Web: http://bioinfo2.ugr.es
#         CGI: http://bioinfo2.ugr.es/GenomeCluster


#  For questions, feedback, etc please contact to: Francisco Dios (frankdiox@gmail.com)
#                                                  Gillermo Barturen (bartg01@gmail.com)
#                                                  Michael Hackenberg (mlhack@gmail.com)
#                                                  José L. Oliver (oliver@ugr.es)

# To see the options of the algorithm, please launch GenomeCluster without any command line arguments

# If you use GenomeCluster, please cite

# Dios F., G. Barturen, R. Lebrón, A. Rueda, M. Hackenberg and J.L. oliver. 2014.
# DNA clustering and genome complexity
# Computational Biology and Chemistry (in preparation)

############################################################################
############################################################################

use strict;


my ($seqst, $ID, $seqlen, $seqlenbruto, $elemnr, $prob, $output, $assembly_dir,
  $d, $plimit, $getd, @dd, @getd_hash ,@dist_n,$chrom_intersec,$genome_intersec,
  $bed_file, $nBedFile, %n_coords, %chrom_bed, %all_dist_n, %n_chrom_index, $elemmn, $cMethod);

#################################
### Parameters ##################
#################################
my $maxN = 0; ## maximal number of Ns
#################################

my @dist_all; # holds distances for the genome
my $elemnr_genome; # number of genomic elements in genome
my $length_genome; # number of dinucleotides in genome
my $elmean_genome; # genomic element mean length in genome
my $gc_name = "_GenomeCluster.txt";

&GetDefault();

####


print "      ***                Getting    Coordinates              ***\n";
GetBED();


####

	$elemnr_genome = 0;
	$length_genome = 0;
  $elmean_genome = 0;

	print "\n      --- Single chromosome:\n";

	for my $chrom (keys %chrom_bed){ ### MAIN LOOP
		print "\n      ".$chrom.":\n";

    print "      ***             Calculating   Seq   Features           ***\n";
    my ($num_dinuc,$last_chunk);
    $elemnr = $#{ $chrom_bed{$chrom}{'cStart'} } + 1; #numero de entradas almacenadas para este chrom
    my $sum1 = 0;
    for (my $i = 0; $i < $elemnr; $i++){
      $sum1 += $chrom_bed{$chrom}{'cEnd'}[$i] - $chrom_bed{$chrom}{'cStart'}[$i];
    }
    my $elmean = int(($sum1/$elemnr) + 0.5); #longitud media de los elementos genomicos
    if($assembly_dir){ #suponemos que nos dan la secuencia
      my $file = $assembly_dir.'/'.$chrom.'.fa';
      print "      ***                  Reading   Sequence                ***\n";
      ($seqst, $ID) = &GetFA($file); #seqst es toda la secuencia sin \n, ID es lo que hay detras de ">"

		  ($num_dinuc, $seqlen) = &GetSeqFeat($seqst);
      $seqlenbruto = length($seqst);
      $last_chunk = $seqlenbruto - $chrom_bed{$chrom}{'cEnd'}[$elemnr - 1];
		}else{ #ESTIMACIONES
      my $sum2 = $chrom_bed{$chrom}{'cStart'}[0];
      for (my $i = 1; $i < $elemnr; $i++){
        my $aux = $chrom_bed{$chrom}{'cStart'}[$i] - $chrom_bed{$chrom}{'cEnd'}[$i-1];
        $sum2 += (($aux > 0) ? $aux : 0);
      }
      $last_chunk = (int(($sum2/$elemnr) + 0.5)); # se hace un round()
      $seqlen = $chrom_bed{$chrom}{'cEnd'}[$elemnr - 1] + $last_chunk;
      $num_dinuc = $seqlen - 1;
      $seqlenbruto = $seqlen;

    }
    my $nm = 0;
    if ($nBedFile){
      $nm = matchN($chrom, $chrom_bed{$chrom}{'cEnd'}[$elemnr - 1], $seqlenbruto);
    }
    push @{ $all_dist_n{$chrom} }, $last_chunk if $nm <= $maxN;

		if($genome_intersec){ #genome
			$elemnr_genome += $elemnr;
			$length_genome += $num_dinuc;
      $elmean_genome += $elmean; # Se divide mas adelante
		}

		## Prob GE ##
    $prob = GetProb($num_dinuc, $elemnr, $elmean);


		#** Begin: single chrom INTERSEC
		if($chrom_intersec){

			print "      ***          Calculating  Chrom  Intersection          ***\n";
			my ($max, $min) = getMinMaxDistance($all_dist_n{$chrom}, $prob);
			$d = $max;
			$getd = "Chromosome Intersection";
			### get protoislands
			print "      ***        Detecting clusters  (Chrom Intersec)        ***\n";
			my @protoislas = &GetProtoIslas($chrom_bed{$chrom}{'cStart'},$chrom_bed{$chrom}{'cEnd'});
			print "      ***       Calculating P-values  (Chrom Intersec)       ***\n";
			@protoislas = &CalcPvalNB(\@protoislas,$prob,$elmean);
			## Get Features like the obs/esp, clustering etc....
			#@protoislas = &GetCGI_features(\@protoislas);

			## Writing output
      $elemmn = $elmean;
			&OUT_f(\@protoislas, $chrom, $output.$chrom."_chromIntersec".$gc_name,$output."all"."_chromIntersec".$gc_name);
		}

		#** End: single chrom INTERSEC




		#** Begin: single chrom PERCENTILE
		if(@getd_hash > 0){
			@dd = sort {$a <=> $b} @{ $all_dist_n{$chrom} };

			foreach(@getd_hash){
				$getd = $_;
				print "      ***               Calculating Percentile $getd            ***\n";
				$d = &GetPerc(\@dd,$_);

				### get protoislands
				print "      ***               Detecting clusters (p$getd)             ***\n";
				my @protoislas = &GetProtoIslas($chrom_bed{$chrom}{'cStart'},$chrom_bed{$chrom}{'cEnd'});
				print "      ***              Calculating P-values  (p$getd)           ***\n";
				@protoislas = &CalcPvalNB(\@protoislas,$prob,$elmean);

				## Get Features like the obs/esp, clustering etc....
				#@protoislas = &GetCGI_features(\@protoislas);

				## Writing output
        $elemmn = $elmean;
        &OUT_f(\@protoislas, $chrom, $output.$chrom."_p".$getd.$gc_name,$output."all"."_p".$getd.$gc_name);


			}
		}

		#** End: single chrom PERCENTILE

		push @dist_all,@{ $all_dist_n{$chrom} } if($genome_intersec); #genome
    undef (@{ $all_dist_n{$chrom} });

	}


	#** Begin: genome
	if($genome_intersec){
		print "\n\n      --- Genome:\n\n";

    $elmean_genome /= (scalar keys %chrom_bed);
    $prob = GetProb($length_genome, $elemnr_genome, $elmean_genome);

    $elemmn = $elmean_genome;
    $elemnr = $elemnr_genome;

    print "      ***          Calculating  Genome Intersection          ***\n\n";
		my ($max, $min) = getMinMaxDistance(\@dist_all, $prob);

    $d = $max;
    $getd = "Genome Intersection";
		foreach my $chrom (keys %chrom_bed){
			print "      ".$chrom."\n";

			### get protoislands
			print "      ***        Detecting clusters (Genome Intersec)        ***\n";
			my @protoislas = &GetProtoIslas($chrom_bed{$chrom}{'cStart'},$chrom_bed{$chrom}{'cEnd'});
			print "      ***       Calculating P-values (Genome Intersec)       ***\n\n";
			@protoislas = &CalcPvalNB(\@protoislas,$prob,$elmean_genome);
			## Get Features like the obs/esp, clustering etc....
			#@protoislas = &GetCGI_features(\@protoislas);

			## Writing output

      &OUT_f(\@protoislas, $chrom, $output.$chrom."_genomeIntersec".$gc_name,$output."all"."_genomeIntersec".$gc_name);
		}
	}

#** End: genome



####################################################################
#######   SUBFUNCTIONS   ###########################################
###################################################################

sub GetDefault{
  print "\n";
  print "---------------------------------------------------------------------------\n";
  print "---------------------------------------------------------------------------\n";
  print "---------                                                         ---------\n";
  print "---------     Computational Genomics and Bioinformatics Group     ---------\n";
  print "---------        University of Granada, Dept. of Genetics         ---------\n";
  print "---------                                                         ---------\n";
  print "---------               Web: http://bioinfo2.ugr.es               ---------\n";
  print "---------        CGI: http://bioinfo2.ugr.es/GenomeCluster        ---------\n";
  print "---------                                                         ---------\n";
  print "---------              GenomeCluster (1.0) 11/30/13               ---------\n";
  print "---------                                                         ---------\n";
  print "---------------------------------------------------------------------------\n";
  print "---------------------------------------------------------------------------\n";
  print "\n";


	if($#ARGV < 3){
	    print "Example of use:\n\n";
      print "perl GenomeCluster.pl <cMethod> <BED> <d> <P-value> [<assembly> [<N_BED> [<maxN>]]]\n";

      print "\ncMethod:   Clustering Method. Type \"element\", \"start\", \"middle\" or \"end\"\n";
      print "             in order to select the method to find clusters: taking into\n";
      print "             account the whole elements or the start, middle or end\n";
      print "             coordinates of each element, respectively.\n";
      print "\nBED:       File input in BED format. Input elements do not need to be sorted\n";
      print "             nor merged. This program will sort, merge and select the input\n";
      print "             depending on the arguments.\n";
	    print "\nd:         The threshold distance on basis of a given percentile.\n";
	    print "             For example: d=25 calculates the percentile 25 of the genomic\n";
	    print "             distance distribution and takes this value as the threshold\n";
	    print "             distance.\n";
	    print "             The recommended value is 50 (median distance).\n";
		  print "             You can add multiple comma-separated percentile values, \"ci\"\n";
		  print "             (chromosome intersection) or \"gi\" (genome intersection).\n";
	   	print "             Example: gi,25,60,ci,50\n";
	    print "\nP-value:   The maximal P-value under which an element group is considered\n";
	    print "             as a cluster.\n";
	    print "             The recommended limit is 1E-5.\n\n";
      print "\nassembly:  Directory containing sequence files in FASTA format. If none is\n";
      print "             provided the program will estimate several parameters.\n";
      print "\nN_BED:     File containing coordinates of N blocks in BED format. If none\n";
      print "             is provided it will assume that sequences do not contain any N.\n";
      print "\nmaxN:      Maximum number of Ns allowed. Default is 0.\n";


	    die "\n";
	}

  #DEFAULT
  $assembly_dir = 0;
  $maxN = 0;
  $nBedFile = 0;


  my $i= 0;

  if ($ARGV[$i] eq "element" or $ARGV[$i] eq "start" or $ARGV[$i] eq "middle" or $ARGV[$i] eq "end"){
    $cMethod = $ARGV[$i];
  }else{
    die "Could not recognize \"$ARGV[$i]\"\n";
  }


  $i++;
  if(-e $ARGV[$i]){
    $bed_file = $ARGV[$i];
  }  else{
    die "Cannot find the input file: $ARGV[$i]\n";
  }

  $i++;
    #ARGV1 - Percentile/ci/gi
	$getd = $ARGV[$i];
	foreach(split(',',$ARGV[$i])){
		if(/\D/){
			if(lc eq 'ci'){
				$chrom_intersec = 1;
			}elsif(lc eq 'gi'){
				$genome_intersec = 1;
			}
		}elsif($_ < 0 or $_ > 100){
			die "The Percentile must be between 0 and 100\n";
		}else{
			push @getd_hash,$_;
		}
	}

  $i++;
  	#ARGV2 - pLimit
    $plimit = $ARGV[$i];
	die "The maximal P-value you have choosen is higher than 1!\nPlease revise the order of the input parameters\n" if($plimit > 1);

  if($#ARGV > $i){
    $i++;
    #ARGV3 - assembly/directory
    if(-d $ARGV[$i]){
      $assembly_dir = $ARGV[$i];

    }else{
      die "Cannot find the input directory: $ARGV[$i]\n";
    }
    if($#ARGV > $i){
      $i++;
      if(-e $ARGV[$i]){
        $nBedFile = $ARGV[$i];
        getNbed($nBedFile);
      }else{
        die "Cannot find the inputo file: $ARGV[$i]";
      }

      if($#ARGV > $i){
        $i++;
        $maxN = $ARGV[$i];
        die "Incorrect value for maxN" if ($maxN < 0);

      }
    }
  }

  $gc_name = "_".$cMethod.$gc_name;

  #Check output names
  $output = $bed_file."_GCresult/";
  if (!(-d $output)){
    mkdir $output, 0755 or die "Cannot create output directory\n";
  }else{ #These files need to be empty (It will concatenate them, not trucante)
    if($chrom_intersec){
      my $filecheck = $output.'all'.'_chromIntersec'.$gc_name;
      unlink $filecheck if(-e $filecheck);
    }
    if($genome_intersec){
      my $filecheck = $output.'all'.'_genomeIntersec'.$gc_name;
      unlink $filecheck if(-e $filecheck);
    }
    foreach my $f (@getd_hash){
      my $filecheck = $output.'all'."_p".$f.$gc_name;
      unlink $filecheck if(-e $filecheck);
    }
  }

  #Data input
  my $command;

  if ($cMethod eq "start"){
    $command = "cut -f1,2,3 ".$bed_file." | sort -u -k1,1V -k2,2n "."| cut -f1,2";
  }elsif ($cMethod eq "end"){
    $command = "cut -f1,2,3 ".$bed_file." | sort -u -k1,1V -k3,3n "."| cut -f1,3";
  }elsif ($cMethod eq "element"){
    $command = "cut -f1,2,3 ".$bed_file." | sort -k1,1V -k2,2n "."| mergeBed -i";
  }elsif ($cMethod eq "middle"){
    $command = "cut -f1,2,3 ".$bed_file." | sort -u -k1,1V -k2,2n "."| awk '{printf(\"%s\t%i\\n\",\$1,\$2+((\$3 - \$2)/2));}' "."| sort -k1,1V -k2,2n ";
  }

  my $bed_file_aux = (split('/',$bed_file))[-1];
  $command .= " > ".$output.$bed_file_aux."_tmp";
  system($command);
  die "Command failed: $!"." - You may need to install 'cut', 'sort', 'awk', 'bedtools' or change Write permissions" if ( $? == -1 );

  $bed_file = $output.$bed_file_aux."_tmp";

}

# Read Sequence
sub GetFA{

  my $seqst_temp = "";
  open (I,$_[0]) or die "Can't open $_[0]";
  my $z = <I>;
  my $tes = substr($z,0,1);
  if($tes ne ">"){
    die "Sequence seems not to be in fasta format !!!!";
  }
  my @z = split(/\s+/," $z");
  $z = $z[1];
  $z=~ s/>//;
  $z =~ s/[\n\t\f\r\s]//g;
  my $ID_temp = $z;
  while($z = <I>){
    $z =~ s/[\n\t\f\r_0-9\s]//g;
    $seqst_temp .= $z;
  }
  close(I);

  return ($seqst_temp,$ID_temp);
}

sub GetBED{

  my $last_coord = 0;
  open (B,$bed_file) or die "Can't open $bed_file";

  while ( my $line = <B> ) {
     chomp($line);
     my $rec = {};
     my @recsplit = split("\t", $line);
     my $chrom = $recsplit[0];
     my $refPoint;
     my $refPoint2;

     $last_coord = 0 if(!exists($chrom_bed{$chrom})); #Reinicializar para cuando cambie de cromosoma

     if($cMethod eq "element"){
       $refPoint = $recsplit[1];
       $refPoint2 = $recsplit[2];

     }elsif($cMethod eq "start"){
       $refPoint = $recsplit[1];
       $refPoint2 = $refPoint + 1;

     }elsif($cMethod eq "middle"){
       #$refPoint= (@{ $chrom_bed{$chrom}{'cStart'} }) + int(((@{ $chrom_bed{$chrom}{'cEnd'} }) - (@{ $chrom_bed{$chrom}{'cStart'} }))/2);
       $refPoint = $recsplit[1];
       $refPoint2 = $refPoint + 1;

     }elsif($cMethod eq "end"){
       $refPoint = $recsplit[1] - 1;
       $refPoint2 = $recsplit[1];
     }

     push @{ $chrom_bed{$chrom}{'cStart'} }, $refPoint;
     push @{ $chrom_bed{$chrom}{'cEnd'} }, $refPoint2;

     #push @{ $chrom_bed{$chrom}{'id'} }, $recsplit[3];
     #push @{ $chrom_bed{$chrom}{'score'} }, $recsplit[4];
     #push @{ $chrom_bed{$chrom}{'strand'} }, $recsplit[5];

     if($refPoint > $last_coord){
       my $nm = 0;
       if ($nBedFile){
        $nm = matchN($chrom, $last_coord, $refPoint);
       }

       push @{ $all_dist_n{$chrom} }, ($refPoint - $last_coord) if($nm <= $maxN);
     }else{
       push @{ $all_dist_n{$chrom} }, 0;
     }
     $last_coord = $refPoint2;

  }

  close(B);
  unlink $bed_file;
}


sub OUT_f {
  my $c=0;

  open (OO,">".$_[2]) or die "Cannot not open ".$_[2];
  open (OOgb,">>".$_[3]) or die "Cannot not open ".$_[3];
  #print OO "chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tlength\tcount\tPValue\tlogPValue\n";
  print OO "Chrom\tFrom\tTo\tLength\tCount\tPValue\tlogPValue\n";

  while($_[0]->[$c]){
	my $log_pvalue = ($_[0]->[$c]->[4] == 0 ? 0 : (log($_[0]->[$c]->[4])/log(10)));
  my $score = 0;
  my $strand = '.';
	#my $patden = (($_[0]->[$c]->[3])*$elemmn/$_[0]->[$c]->[2]);
	printf OO "%s\t%i\t%i\t%i\t%i\t%.2e\t%.2f\n", $_[1], $_[0]->[$c]->[0], $_[0]->[$c]->[1], $_[0]->[$c]->[2], $_[0]->[$c]->[3], $_[0]->[$c]->[4],$log_pvalue;
  #printf OO "%s\t%i\t%i\t%s\t%i\t%s\t%i\t%i\t%.2e\t%.2f\n", $_[1], $_[0]->[$c]->[0], $_[0]->[$c]->[1], $_[1]."_".$_[0]->[$c]->[0],$score,$strand,$_[0]->[$c]->[2], $_[0]->[$c]->[3], $_[0]->[$c]->[4],$log_pvalue;
  #rintf OOgb "%s\t%i\t%i\t%s\t%i\t%s\t%i\t%i\t%.2e\t%.2f\n", $_[1], $_[0]->[$c]->[0], $_[0]->[$c]->[1], $_[1]."_".$_[0]->[$c]->[0],$score,$strand,$_[0]->[$c]->[2], $_[0]->[$c]->[3], $_[0]->[$c]->[4],$log_pvalue;
  printf OOgb "%s\t%i\t%i\t%s\t%i\t%s\n", $_[1], $_[0]->[$c]->[0], $_[0]->[$c]->[1], $_[1]."_".$_[0]->[$c]->[0],$score,$strand;

    $c++;
  }
  close(OO);
  close(OOgb);

  open(O,">".$_[2]."-log.txt") or die "can't open $_[2]-log.txt";
  print O "Basic statistics of the input sequence: $_[1]"."_".$getd."\n";
  printf O "Length: %d\n",$seqlen;
  printf O "Length including Ns: %d\n",$seqlenbruto;
  #my $fg = $seqst =~ s/g/g/ig;
  #my $fc = $seqst =~ s/c/c/ig;
  #my $fa = $seqst =~ s/a/a/ig;
  #my $ft = $seqst =~ s/t/t/ig;
  #printf O "GC content: %0.3f\n",100*($fg+$fc)/$seqlen;
  printf O "Number of elements in sequence: %d\n",$elemnr;
  printf O "Element mean length in sequence: %d\n",$elemmn;
  printf O "Probability to find a genomic element: %.4f\n\n",$prob;
  printf O "Calculated distance threshold: %i\n",$d;
  print O "Parameters used:\n";

  printf O "p-value threshold: $plimit\n";
  if($getd){
    print O "Distance threshold method: ";
    $_ = $getd;
    print O "percentile " if(/\d/);
    print O "$getd\n";
  }

  close(O);
}


## Get CpG cluster
sub GetProtoIslas{

  my @cStart = @{$_[0]};
  my @cEnd = @{$_[1]};
  my @t;
  my ($start, $end, $elementnr);
  my $des = "no";
  for(my $i = 0; $i <= $#cStart - 1; $i++){

    my $dist = $cStart[$i+1]  - $cEnd[$i]; ## Puede ser negativo si solapan sin problema

    if($dist <= $d){
      if($des eq "no"){
	      $start = $cStart[$i];
        $elementnr = 1;
      }
      $end = $cEnd[$i+1];
      $elementnr++;
      $des = "yes";
    }
    elsif($dist > $d and $des eq "yes"){
      $des = "no";
      my @f = ($start, $end, $elementnr);
      push @t,\@f;
    }
  }
  if($des eq "yes"){
    my @f = ($start, $end, $elementnr);
    push @t,\@f;
  }
  return @t;
}

sub GetCGI_features{ # Deprecated

  my @temp;
  my $c=0;
  while(defined($_[0]->[$c])){
    if($_[0]->[$c]->[4] < $plimit){
      my $len = $_[0]->[$c]->[1] - $_[0]->[$c]->[0] +1;
      (my $oe, my $cpgseq, my $gccont)= &CalcObsEsp($_[0]->[$c]->[0] -1,$len);
      my $coord1 = &GetCoord($cpgseq);
      (my $clust, my $meandist) = &GetClust($coord1);

      my $pval = $_[0]->[$c]->[4];
      $_[0]->[$c]->[4] = $oe;
      $_[0]->[$c]->[5] = $meandist;
      $_[0]->[$c]->[6] = $clust;
      $_[0]->[$c]->[7] = $gccont;
      $_[0]->[$c]->[8] = $pval;
      my @t = @{$_[0]->[$c]};
      push @temp,\@t;
    }
    $c++;
  }
  return @temp;
}

sub CalcObsEsp{

  my $cpgseq = substr ($seqst,$_[0],$_[1]);
  my $fc = $cpgseq =~ s/c/c/ig;
  my $fg = $cpgseq =~ s/g/g/ig;
  my $CGICpG = $cpgseq =~ s/CG/CG/ig;
  my $e = $fc*$fg;
  my $oet = $CGICpG*length($cpgseq)/$e;
  return ($oet,$cpgseq,$fc+$fg);
}

sub GetSeqFeat{

  my $n = $_[0];
  #my $elemnr = $n =~ s/CG/Cg/ig;
  my $ln = length($n);
  my $NN = $n =~ s/N/N/ig;
  my $seqlen = $ln-$NN;
  my $num_dinuc = 0;


  my $elem;
  my $limit = $ln - 1;
  for(my $i = 0; $i < $limit; $i++){
  	$elem = substr($n,$i,2);
  	if($elem !~ m/N/i){
  		$num_dinuc++;
  	}
  }
  return ($num_dinuc,$seqlen);
}

sub GetProb{
  my ($length, $elemnr, $elmean) = @_;
  return (($elemnr * $elmean)/$length);
}


sub GetPerc{

  my @t = @{$_[0]};
  my $totnr = @t;
  for(my $i = 0; $i <= $#t; $i++){
    if(100*$i/$totnr >= $_[1]){
      return $t[$i];
    }
  }
  return $t[$#t];
}

sub getNbed{
  open(N, '<'.$_[0]) or die "Cannot open file input: $_[0]\n";
    while ( my $line = <N> ) {
      my $rec = {};
      my @recsplit = split("\t", $line);
      my $chrom = $recsplit[0];

      if(!exists($n_coords{$chrom})){
        $n_coords{$chrom} = [];
      }

      push (@{$n_coords{$chrom}}, [$recsplit[1], $recsplit[2]]);
  }

  close(N);

}

sub matchN{
  my ($chrom, $cStart, $cEnd) = @_;

  foreach my $elem (@{$n_coords{$chrom}}){

    if(($cStart >= $elem->[0] && $cStart <= $elem->[1]) && ($cEnd >= $elem->[0] && $cEnd <= $elem->[1])){
      return ($cEnd - $cStart);
    }elsif($cStart >= $elem->[0] && $cStart <= $elem->[1]){
      return ($elem->[1] - $cStart);
    }elsif($cEnd >= $elem->[0] && $cEnd <= $elem->[1]){
      return ($cEnd - $elem->[0]);
    }elsif($cStart <= $elem->[0] && $cEnd >= $elem->[1]){
      return ($elem->[1] - $elem->[0]);
    }

  }
  return 0; #No hay N
}


##########################################################################
############## SUBFUNCTIONS for P-value calculations  ####################
##########################################################################

######################################################################
## *** Calculates the negative binomial distribution

sub CalcPvalNB{

  my ($pval, @temp);
  my @islas = @{$_[0]};
  my $prob_ge = $_[1];
  my $elmean = $_[2];

  for(my $i = 0; $i <= $#islas; $i++){
    my $l = $islas[$i]->[1] - $islas[$i]->[0];
          #my $str = substr ($seqst,$islas[$i]->[0]-1,$l);
          #my $cpg = $str =~ s/cg/cg/ig;
    my $gel = $islas[$i]->[2];
          #my $pval = &GetNB({{{$l-(2*$cpg}}}),$cpg-1,$_[1]);
    my $nf = int(($l-($elmean*$gel))/$elmean) + 1; #ceiling
    $nf = (($nf > 0) ? $nf : 0);
    my $pval = &GetNB($nf,$gel-1,$prob_ge);
    $pval = sprintf("%.5e",$pval);
    next if $pval > $plimit; # Uso de PLIMIT
    my @t = ($islas[$i]->[0],$islas[$i]->[1],$l,$gel,$pval);
    push @temp, \@t;
  }

  return @temp;
}
sub GetNB{
  my $pval = 0;
  for(my $j = 0; $j <= $_[0]; $j++){
    my $ptemp = &FactorialNB($j,$_[1]) + $_[1]*log($_[2]) + $j*log(1.0-$_[2]);
    $ptemp = exp($ptemp);
    $pval += $ptemp;
    }
  return $pval;
}


sub FactorialNB{
  my $stop = $_[0]+$_[1]-1;
  my $l1 = 0;
  my $l2 = 0;
  for(my $i = $_[0]+1;$i <= $stop; $i++){
    $l1 +=log($i);
  }
  for(my $i = 1;$i <= $_[1]-1; $i++){
    $l2 +=log($i);
  }
  return $l1-$l2;
}
##########################################################################
############## SUBFUNCTIONS for Clustering Calculation  ####################
##########################################################################

sub GetClust{

  my @t = @{$_[0]};
  my $dist = &GetDist($_[0]);
  my $mean = &Normalize($dist);
  return (1, $mean);
}

sub Normalize{

  my @d = @{$_[0]};
  my $tot;
  for(my $i = 0; $i <= $#d; $i++){
    $tot+=$d[$i];
  }
  my $mean = $tot/@d;
  return $mean;
}

sub GetDist{

  my @dist;
  my @d = @{$_[0]};
  for(my $i = 0; $i < $#d; $i++){
    my @f = split (/\-/,$d[$i]);
    my @s = split (/\-/,$d[$i+1]);
    my $dist = $s[0]-$f[1];
    push @dist,$dist;
  }
  return \@dist;
}

sub GetCoord{

  my @coord;
  my @c = split (//,$_[0]);
  for(my $i = 0; $i < $#c; $i++){
    if($c[$i] =~ /[cC]/ and $c[$i+1] =~ /[gG]/){
      my $str = $i.'-'.eval($i+1);
      push @coord,$str;
    }
  }
  return \@coord;
}


sub getMinMaxDistance{
	my @distances = @{$_[0]};
	my $prob = $_[1];

	my @distCount = ();
	my $maxDist = 0;
	my $nrDistances = 0;
	my $stop = @distances;

	for(my $i = 0; $i < $stop; $i++){
		$nrDistances++;
		my $dist = $distances[$i];
		if(defined $distCount[$dist]){
			$distCount[$dist]++;
		}else{
			$distCount[$dist] = 1;
		}

		if($dist > $maxDist){
			$maxDist = $dist;
		}
	}

	my $obsCum = 0;
	my $teoCum = 0;
	my @dif = ();
	for(my $i = 1; $i <= $maxDist; $i++){
		if(defined $distCount[$i]){
			$obsCum += $distCount[$i]/$nrDistances;
		}
		$teoCum += $prob * (1 - $prob)**($i-1);
		$dif[$i] = ($obsCum - $teoCum);
	}

	my $max = -1;
	my $min = 1;
	my $maxD = 0;
	my $minD = 0;
	for(my $i = 1; $i <= $maxDist; $i++){
		if(defined $dif[$i]){
			my $difT = $dif[$i];
			if($difT > $max){
				$max = $difT;
				$maxD = $i;
			}
			elsif($difT < $min){
				$min = $difT;
				$minD = $i;
			}
		}
	}
	return ($maxD, $minD);
}
