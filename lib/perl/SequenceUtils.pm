package VEuPath::SequenceUtils;

use strict;

my $debug = 0;

sub breakSequence {
  my($seq,$lineLength,$beginSpace) = @_;
	my $ll = $lineLength ? $lineLength : 80;
  ##just in case there are returns...
  $seq =~ s/\s//g;
  my $new = "";
  for (my $i = 0;$i<length($seq);$i+=$ll) {
    $new .= $beginSpace . substr($seq,$i,$ll) . "\n";
  }
  return $new;
}

sub makeFastaFormattedSequence{
	my($defline,$sequence,$lineLength) = @_;
	my $ll = $lineLength ? $lineLength : 80;
	return ">" . $defline . "\n" . &breakSequence($sequence,$ll);
}

sub reverseComplementSequence{	##for reverseComplementing sequences
  my($seq) = @_;
  $seq =~ s/\s//g;

  print STDERR "revCompSeq: incoming:\n$seq\n" if $debug == 1;

  my $revcompseq = reverse $seq;
  $revcompseq =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;

  print STDERR "revCompSeq: returns:\n$revcompseq\n" if $debug == 1;
  return $revcompseq;
}

sub compNuc{
  my($nuc) = @_;
  if ($nuc =~ /A/i) {
    return "T";
  } elsif ($nuc =~ /T/i) {
    return "A";
  } elsif ($nuc =~ /C/i) {
    return "G";
  } elsif ($nuc =~ /G/i) {
    return "C";
  }
  return $nuc;									## - and N get returned as themselves
}

##Methods for doing translation ##

my @names; 
my	@AA;
my	@SC;
my	@B1;
my	@B2;
my	@B3;

$names[1] = "Standard";
$AA[1] = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
$SC[1] = "---M---------------M---------------M----------------------------";
$B1[1] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B2[1] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B3[1] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

$names[2] = "Vertebrate Mitochondrial";
$AA[2] = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
$SC[2] = "--------------------------------MMMM---------------M------------";
$B1[2] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B2[2] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B3[2] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
 
$names[3] = "Yeast Mitochondrial";
$AA[3] = "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
$SC[3] = "----------------------------------MM----------------------------";
$B1[3] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B2[3] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B3[3] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
  
$names[4] = "Mold, Protozoan, Coelenterate Mitochondrial and Mycoplasma/Spiroplasma";
$AA[4] = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
$SC[4] = "--MM---------------M------------MMMM---------------M------------";
$B1[4] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B2[4] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B3[4] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

$names[5] = "Invertebrate Mitochondrial";
$AA[5] = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG";
$SC[5] = "---M----------------------------MMMM---------------M------------";
$B1[5] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B2[5] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B3[5] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

$names[6] = "Ciliate, Dasycladacean and Hexamita Nuclear";
$AA[6] = "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
$SC[6] = "-----------------------------------M----------------------------";
$B1[6] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B2[6] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B3[6] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

$names[9] = "Echinoderm and Flatworm Mitochondrial";
$AA[9] = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
$SC[9] = "-----------------------------------M---------------M------------";
$B1[9] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B2[9] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B3[9] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

$names[10] = "Euplotid Nuclear";
$AA[10] = "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
$SC[10] = "-----------------------------------M----------------------------";
$B1[10] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B2[10] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B3[10] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

$names[11] =  "Bacterial, Archaeal and Plant Plastid";
$AA[11] = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
$SC[11] = "---M---------------M------------MMMM---------------M------------";
$B1[11] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B1[11] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B1[11] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
  
$names[12] = "Alternative Yeast Nuclear";
$AA[12] = "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
$SC[12] = "-------------------M---------------M----------------------------";
$B1[12] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B2[12] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B3[12] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

$names[13] =  "Ascidian Mitochondrial";
$AA[13] = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG";
$SC[13] = "---M------------------------------MM---------------M------------";
$B1[13] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B2[13] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B3[13] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

$names[14] =  "Alternative Flatworm Mitochondrial";
$AA[14] = "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
$SC[14] = "-----------------------------------M----------------------------";
$B1[14] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B2[14] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B3[14] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

$names[16] =  "Chlorophycean Mitochondrial";
$AA[16] = "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
$SC[16] = "-----------------------------------M----------------------------";
$B1[16] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B2[16] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B3[16] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";


$names[21] =  "Trematode Mitochondrial";
$AA[21] = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
$SC[21] = "-----------------------------------M---------------M------------";
$B1[21] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B2[21] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B3[21] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";


$names[22] =  "Scenedesmus obliquus Mitochondrial";
$AA[22] = "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
$SC[22] = "-----------------------------------M----------------------------";
$B1[22] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B2[22] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B3[22] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

$names[23] =  "Thraustochytrium Mitochondrial";
$AA[23] = "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
$SC[23] = "--------------------------------M--M---------------M------------";
$B1[23] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B2[23] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B3[23] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

$names[24] =  "Pterobranchia Mitochondrial";
$AA[23] = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG";
$SC[24] = "---M---------------M---------------M---------------M------------";
$B1[24] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B2[24] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B3[24] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

$names[25] =  "Candidate Division SR1 and Gracilibacteria";
$AA[25] = "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
$SC[25] = "---M-------------------------------M---------------M------------";
$B1[25] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
$B2[25] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
$B3[25] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

my %organismCode = ("Standard" => 1,
										"Vertebrate Mitochondrial" => 2,
										"Yeast Mitochondrial" => 3,
										"Mold, Protozoan, Coelenterate Mitochondrial and Mycoplasma/Spiroplasma" => 4,
										"Invertebrate Mitochondrial" => 5,
										"Ciliate, Dasycladacean and Hexamita Nuclear" => 6,
										"Echinoderm and Flatworm Mitochondrial" => 9,
										"Euplotid Nuclear" => 10,
										"Bacterial, Archaeal and Plant Plastid" => 11,
										"Alternative Yeast Nuclear" => 12,
										"Ascidian Mitochondrial" => 13,
										"Alternative Flatworm Mitochondrial" => 14,
										"Chlorophycean Mitochondria" => 16, 
										"Trematode Mitochondrial" => 21, 
                                                                                "Scenedesmus obliquus Mitochondrial" => 22, 
                                                                                "Thraustochytrium Mitochondrial" => 23, 
                                                                                "Pterobranchia Mitochondrial" => 24, 
                                                                                "Candidate Division SR1 and Gracilibacteria" => 25 );

sub generateCodonHash {
	my($index) = @_;
	my $i = $index ? $index : 1;
	my @aa = split('',$AA[$i]);
	my @b1 = split('',$B1[$i]);
	my @b2 = split('',$B2[$i]);
	my @b3 = split('',$B3[$i]);
	my %cod;
	for (my $a = 0;$a<scalar(@aa);$a++){
		my $codon = $b1[$a] . $b2[$a] . $b3[$a];
		$cod{$codon} = $aa[$a];
	}
	return %cod;
}

sub translateSequence {
	my($seq,$codonTable) = @_;
	my $ct = $codonTable ? $codonTable : 1;  ##default is the standard..
	my %cod = &generateCodonHash($ct);
#	print STDERR "Codons: (".join(', ',keys%cod).")\n";
	my $trans;
	for(my $i = 0;$i<length($seq)-2;$i+=3){
		my $aa = $cod{substr($seq,$i,3)};
		$trans .= $aa ? $aa : "X";
	}
	return $trans;
}

1;
