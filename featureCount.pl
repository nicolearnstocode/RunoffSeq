#featureCount.pl -m GRCh38_latest_rna.fna -p GRCh38_latest_protein.faa -l1 GR101-pro.txt,2,test -id idlist.txt -o 0
#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $usage = <<USAGE;

perl featureCount.pl <parameters>
    -m 	<str>	the mRNA file
    -p 	<str>	the peptide file
    -l1 <str>	the id list1(name.txt,1,treated) eg. 1-geneid,2-genename,3-uniprotid
    -id <str>   the id file
    -o  <num>   1 or 0 for output geneSname.list

USAGE

my %enGeneID2pid;
my %GeneName2pid;
my %Uniprot2pid;

my %enGeneID2mid;
my %GeneName2mid;
my %Uniprot2mid;

my $mrnaFile;
my $proteinFile;
my $list1;
my $idfile;
my $output;

GetOptions ("m=s" => \$mrnaFile,
            "p=s" => \$proteinFile,
            "l1=s" => \$list1,
            "id=s" => \$idfile,
            "o=i" => \$output
            );

if ($mrnaFile eq "" or $proteinFile eq "" or $list1 eq "" or $idfile eq "" or $output eq "") {
    print $usage,"\n";
    exit;
}

loadID($idfile);

my %comlist;
loadList($list1);


my %mRNA;
load_mRNA($mrnaFile);

my %protein;
loadPeptid($proteinFile);

if($output == 1){
	open(FF,">geneSname.list") || die;
}

#featureScan('AAAAAAAAAAAA','mrna',50);
#featureScan('[G|C]{15,}','mrna',50);
#featureScan('GATAAG','mrna',50);
#featureScan('GACAAA','mrna',50);
#featureScan('GGAGGA','mrna',50);
featureScan('KKK','protein',16);
#featureScan('DK','protein',16);
#featureScan('QQQ','protein',16);
#featureScan('GGG','protein',16);
#featureScan('P[P|G|D]','protein',16);
#featureScan('R\wK','protein',16);
#featureScan('RKK','protein',16);
#featureScan('GI|DI|GD','protein',16);
#the first 100nt
#freeEnergy(0,60);

close(FF);

#loading id file for id conversion
sub loadID{
    my $file = shift;
    open(ID,"$file") || die;
    while(my $ln=<ID>){
        chomp $ln;
        my @idlist = split(/\,/,$ln);
        if($idlist[3]){
            if($idlist[0]){
                $enGeneID2pid{$idlist[0]} = $idlist[3];
            }
            if($idlist[1]){
                $GeneName2pid{$idlist[1]} = $idlist[3];
            }
            if($idlist[2]){
                $Uniprot2pid{$idlist[2]} = $idlist[3];
            }
        }
        if($idlist[4]){
            if($idlist[0]){
                $enGeneID2mid{$idlist[0]} = $idlist[4];
            }
            if($idlist[1]){
                $GeneName2mid{$idlist[1]} = $idlist[4];
            }
            if($idlist[2]){
                $Uniprot2mid{$idlist[2]} = $idlist[4];
            }
        }

    }
    close(ID);
}

#loading the id list
sub loadList{
    my $infor = shift;
    my ($file,$ind,$tag) = split(/\,/,$infor);
    open(LIST,"$file") || die;
    while(my $ln=<LIST>){
        chomp $ln;
	    $ln = uc($ln);
        if($ind == 1){
            if(exists $enGeneID2pid{$ln}){
                $comlist{'p'}{$enGeneID2pid{$ln}}=1;
                $comlist{'m'}{$enGeneID2mid{$ln}}=1;
            }
        }
        if($ind == 2){
            if(exists $GeneName2pid{$ln}){
                $comlist{'p'}{$GeneName2pid{$ln}}=1;
                $comlist{'m'}{$GeneName2mid{$ln}}=1;
            }
        }
        if($ind == 3){
            if(exists $Uniprot2pid{$ln}){
                $comlist{'p'}{$Uniprot2pid{$ln}}=1;
                $comlist{'m'}{$Uniprot2mid{$ln}}=1;
            }
        }
    }
    close(LIST);

}

sub featureScan{
    my $pattern = shift;
    my $type = shift;
    my $formost = shift;

    my $tCount = 0;
    my $fCount = 0;
    my $mCount = 0;


    if($type eq 'mrna'){
        foreach my $key(sort keys %{$comlist{'m'}}){
            $tCount++;
            if(exists $mRNA{$key}){
                if($mRNA{$key}=~/$pattern/){
                    $fCount++;
                    if($output==1){
						print FF $pattern,"\t",$key,"\tfull\n";
					}
                }
                if(substr($mRNA{$key},0,$formost)=~/$pattern/){
                    $mCount++;
                    if($output==1){
						print FF $pattern,"\t",$key,"\tfront\n";
					}
                }
            }
        }
    }

    if($type eq 'protein'){
        foreach my $key(sort keys %{$comlist{'p'}}){
            $tCount++;
            if(exists $protein{$key}){
                if($protein{$key}=~/$pattern/){
                    $fCount++;
                    if($output==1){
						print FF $pattern,"\t",$key,"\tfull\n";
					}
                }
                if(substr($protein{$key},0,$formost)=~/$pattern/){
                    $mCount++;
                    if($output==1){
						print FF $pattern,"\t",$key,"\tfront\n";
					}
                }
            }
        }
    }

    print $pattern,"\t",$tCount,"\t",$fCount,'(',sprintf("%.2f%%", $fCount/$tCount*100),")\t";
    print $mCount,"(",sprintf("%.2f%%", $mCount/$tCount*100),")\n";
}

sub load_mRNA{
    my $file = shift;
    my $id;
    open(MR,"$file") || die;
    while(my $ln=<MR>){
        chomp $ln;
        if($ln=~/\>/){
            $id = [split(/[\>,\.]/,$ln)]->[1];
        } else {
            $mRNA{$id}.=$ln;
        }
    }
    close(MR);

}

sub loadPeptid{
    my $file = shift;
    my $id;
    open(MR,"$file") || die;
    while(my $ln=<MR>){
        chomp $ln;
        if($ln=~/\>/){
            $id = [split(/[\>,\.]/,$ln)]->[1];
        } else {
            $protein{$id}.=$ln;
        }
    }
    close(MR);
}



sub freeEnergy{
	my $dart = shift;
	my $formost = shift;
	#open(TEST,">test.fa") || die;
	#print TEST '>',"clash".$index,"\n",$clash{$key}{$sky}{$tky},"\n";
	#close(TEST);
	#system("RNAfold -p -d2 --noLP <test.fa >test.out");
	my $totalFe = 0;
	my $totalMr = 0;
	my $freeE;
	   foreach my $key(sort keys %{$comlist{'m'}}){
		   if(exists $mRNA{$key}){
			   $totalMr++;
			   my $seq = substr($mRNA{$key},$dart,$formost);
			   my @out = qx(echo $seq | RNAfold);
			   if($out[1] =~/\((\-\d+\.\d+)\)/){
			  	 $freeE = $1;
			 	  if($output==1){
					  print FF $key,"\t$freeE\n";
				   }
				   $totalFe+=$freeE;
			}
		   }
	   }
	print $totalMr,"(average",sprintf("%.2f%%", $totalFe/$totalMr),")\n";
	return($totalFe);
}
