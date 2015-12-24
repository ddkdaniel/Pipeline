#!/usr/bin/perl -w
#@version 0.1.0
#Branched from bigtable.pl v.0.1.5
#Deal with bug of duplicated mutation_ids
#Re-engineer to use Mutation_ID as hash key to link all three annotation sources
#@version 0.1.1
#    removed ensGene,dbsnp annotations, and other deprecated sift and pp2
#    Annovar empty annotation not using NA anymore
#@version 0.1.2
#    added back dbsnp annotation
#@version 0.1.3
#    remove exon and p. prefix in snpeff transcript field
#
#@version 0.1.4
#    not doing dbsnp annotation
#    update to annovar annotations
#    known issue: DEL is represented as "-" in ALT, incompitable to VCF format
#@version 0.2.0
#    use annovar/2015-03-22 TSV output format
#    use snpEFF 4.1c format
#
use strict;

use Getopt::Std;
use vars qw($opt_o $opt_e $opt_a $opt_c $opt_N $opt_V);
getopts('e:a:o:V!N!c:');

my $VERSION="0.2.0";

my $HELP=qq~$0 [options]
	-e snpEff vcf [required]
	-a annovar table_annovar csv [required]
	-c chr
	-N No table header in output
	-o output tsv [required]
	-V print version
~;

if($opt_V){
	print "Version: $VERSION\n";
	exit(0);
}
die "$HELP\n" unless($opt_e && $opt_a && $opt_o);

my $eff_vcf=$opt_e;
my $annovar=$opt_a;

our %SNPEFF_IMPACT_RANK=('HIGH'=>4,'MODERATE'=>3,'LOW'=>2,'MODIFIER'=>1); 
my @pos=();
my %data=();

my @eff_col=qw(CHR POS REF ALT MutationID snpEff_Impact snpEff_Gene snpEff_Transcript);
my @annovar_col=qw(Func.refGene Gene.refGene ExonicFunc.refGene AAChange.refGene PopFreqMax 1000G_ALL 1000G_AFR 1000G_AMR 1000G_EAS 1000G_EUR 1000G_SAS ExAC_ALL ExAC_AFR ExAC_AMR ExAC_EAS ExAC_FIN ExAC_NFE ExAC_OTH ExAC_SAS ESP6500siv2_ALL ESP6500siv2_AA ESP6500siv2_EA CG46 COSMIC70 NCI60 snp138NonFlagged ljb SIFT_score SIFT_pred pp2_HDIV pp2_HDIV_pred pp2_HAVR pp2_HAVR_pred LRT LRT_pred MutationTaster MutationTaster_pred MutationAssessor MutationAssessor_pred FATHMM FATHMM_pred RadialSVM RadialSVM_pred LR LR_pred VEST3_score CADD_raw CADD_phred GERP_RS PhyloP46way_placental phyloP100way_vertebrate SiPhy29way_logOdds);

open OUT, ">$opt_o" or die "Can't write to $opt_o:$!\n";
&read_eff($eff_vcf);
&read_annovar($annovar);

my @col=(@eff_col, @annovar_col);

my $table = dataToTable(data=>\%data, col=>\@col, key=>\@pos, sep=>"\t");
print OUT $table;


close OUT;

sub read_eff{
	my $in=shift;
	open EFF, "$in" or die "Can't open $in:$!\n";
	while(<EFF>){
		next if(/^##/);
		last if(/^#CHR/);
	}

	while(my $line=<EFF>){
		chomp $line;
		my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info)=split "\t", $line;
		next if($opt_c && $chr ne $opt_c);
	
		my $mutation_id=0;
		my $ann="";
		if($info=~/MutationID=(\d+);ANN=(.*)/){
			($mutation_id, $ann)=($1, $2);
		} else {
			die "Invalid SNPEff format at $chr:$pos - $line\n";
		}
		push @pos, "$chr:$pos:$mutation_id";

		&parseSnpEff($chr, $pos, $ref, $alt, $mutation_id, $ann);
	}
	close EFF;
}


sub read_annovar {
	my $in=shift;
	open ANNOVAR, "$in" or die "Can't open $in:$!\n";
	chomp(my $header=<ANNOVAR>);
	die "I=Incompatible ANNOVAR header:\n$header\n" unless($header eq "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tPopFreqMax\t1000G_ALL\t1000G_AFR\t1000G_AMR\t1000G_EAS\t1000G_EUR\t1000G_SAS\tExAC_ALL\tExAC_AFR\tExAC_AMR\tExAC_EAS\tExAC_FIN\tExAC_NFE\tExAC_OTH\tExAC_SAS\tESP6500siv2_ALL\tESP6500siv2_AA\tESP6500siv2_EA\tCG46\tcosmic70\tnci60\tsnp138NonFlagged\tSIFT_score\tSIFT_pred\tPolyphen2_HDIV_score\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_score\tPolyphen2_HVAR_pred\tLRT_score\tLRT_pred\tMutationTaster_score\tMutationTaster_pred\tMutationAssessor_score\tMutationAssessor_pred\tFATHMM_score\tFATHMM_pred\tRadialSVM_score\tRadialSVM_pred\tLR_score\tLR_pred\tVEST3_score\tCADD_raw\tCADD_phred\tGERP++_RS\tphyloP46way_placental\tphyloP100way_vertebrate\tSiPhy_29way_logOdds\tOtherinfo");

	while(<ANNOVAR>){
		next if($opt_c && $_!~/^$opt_c,/);
		chomp;
		$_.=",";
		&parseAnnovar($_);
	}
	close ANNOVAR;
}


sub parseSnpEff{
	my ($chr, $pos, $ref, $alt, $mutation_id, $ann)=@_;
	my $cur_impact_score=0;
	my $cur_impact="";
	my $cur_gene="";
	my @eff_info=();

##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
#ANN=T|missense_variant|MODERATE|VPS13D|VPS13D|transcript|NM_015378.2|Coding|48/70|c.9799C>T|p.Arg3267Trp|9940/16328|9799/13167|3267/4388||,T|missense_variant|MODERATE|VPS13D|VPS13D|transcript|NM_018156.2|Coding|47/69|c.9724C>T|p.Arg3242Trp|9865/16253|9724/13092|3242/4363||
	my @anns = split ",", $ann;
	foreach my $eff(@anns){
		my($allele,$annotation,$impact, $gene_name, $gene_id, $feature, $feat_id, $transcript_type, $rank_str, $hgvs_coding, $hgvs_prot, $cdna_str, $cds_str, $aa_str, $distance, $err)=split /\|/, $eff;
		my $transcript_id= $feature eq 'transcript' ? $feat_id : "";
		my ($exon_rank, $num_exons)=$rank_str ? split "/", $rank_str : ('','');
		my($cdna_pos, $cdna_len)=$cdna_str ? split "/", $cdna_str : ('','');
		my($cds_pos, $cds_len) = $cds_str ? split "/", $cds_str : ('','');
		my($aa_pos, $aa_len) = $aa_str ? split "/", $aa_str : ('','');		
		if($SNPEFF_IMPACT_RANK{$impact}>$cur_impact_score){
			$cur_impact_score= $SNPEFF_IMPACT_RANK{$impact};
			$cur_impact = $impact;
			$cur_gene = $gene_name;
		}
				
		push @eff_info, [$annotation, $impact, $gene_name, $feature, $transcript_id, $transcript_type, $exon_rank, $num_exons, $hgvs_coding, $hgvs_prot, $cdna_pos, $cdna_len, $cds_pos, $cds_len, $aa_pos, $aa_len, $distance];  
	}
	my $snpeff_info="";
	foreach (@eff_info){ 
		$snpeff_info .= join(":", @$_) . ",";
	}
	$data{"$chr:$pos:$mutation_id"}={MutationID=>$mutation_id, CHR=>$chr, POS=>$pos, REF=>$ref, ALT=>$alt, snpEff_Impact=>$cur_impact, snpEff_Gene=>$cur_gene, snpEff_Transcript=>$snpeff_info };
}


sub parseAnnovar{
	my $val = shift;
	my @col = split "\t", $val;
	my($chr, $pos, $end, $ref, $alt, $Func_refGene, $Gene_refGene, $GeneDetail_refGene, $ExonicFunc_refGene, $AAChange_refGene, $PopFreqMax, $_1000G_ALL, $_1000G_AFR, $_1000G_AMR, $_1000G_EAS, $_1000G_EUR, $_1000G_SAS, $ExAC_ALL, $ExAC_AFR, $ExAC_AMR, $ExAC_EAS, $ExAC_FIN, $ExAC_NFE, $ExAC_OTH, $ExAC_SAS, $ESP6500siv2_ALL, $ESP6500siv2_AA, $ESP6500siv2_EA, $CG46, $cosmic70, $nci60, $snp138NonFlagged, $SIFT_score, $SIFT_pred, $Polyphen2_HDIV_score, $Polyphen2_HDIV_pred, $Polyphen2_HVAR_score, $Polyphen2_HVAR_pred, $LRT_score, $LRT_pred, $MutationTaster_score, $MutationTaster_pred, $MutationAssessor_score, $MutationAssessor_pred, $FATHMM_score, $FATHMM_pred, $RadialSVM_score, $RadialSVM_pred, $LR_score, $LR_pred, $VEST3_score, $CADD_raw, $CADD_phred, $GERP_RS, $phyloP46way_placental, $phyloP100way_vertebrate, $SiPhy29way_logOdds, @otherinfo) = @col;

	my $mutation_id=0;
	if($otherinfo[7]=~/MutationID=(\d+)/){
		$mutation_id=$1;
	} else {
		$mutation_id="NA";
		print STDERR "Missing mutation id in annovar table at $chr, $pos\n";
	}

	my $ljb =($SIFT_pred|| $Polyphen2_HVAR_pred|| $LRT_pred||$MutationTaster_pred||$MutationAssessor_pred||$FATHMM_pred||$RadialSVM_pred||$LR_pred|| $VEST3_score|| $CADD_raw|| $CADD_phred|| $GERP_RS|| $phyloP46way_placental||$phyloP100way_vertebrate|| $SiPhy29way_logOdds) ? 1 :0;
	$pos-- if($alt eq "-"); #DEL - VCF POS = ANNOVAR POS - 1
	my %annovar=(
		CHR=>$chr, POS=>$pos, REF=>$ref, ALT=>$alt, 
		'Func.refGene'      =>$Func_refGene, 
		'Gene.refGene'      =>$Gene_refGene, 
		'GeneDetail.refGene'=>$GeneDetail_refGene,
		'ExonicFunc.refGene'=>$ExonicFunc_refGene, 
		'AAChange.refGene'=>$AAChange_refGene, 
		'PopFreqMax'=>      $PopFreqMax,
		'1000G_ALL' =>      $_1000G_ALL,
		'1000G_AFR' =>      $_1000G_AFR,
		'1000G_AMR' =>      $_1000G_AMR,
		'1000G_EAS' =>      $_1000G_AMR,
		'1000G_EUR' =>      $_1000G_EUR,
                '1000G_SAS' =>      $_1000G_SAS,
		'ExAC_ALL' =>       $ExAC_ALL,
                'ExAC_AFR' =>       $ExAC_AFR,
                'ExAC_AMR' =>       $ExAC_AMR,
                'ExAC_EAS' =>       $ExAC_AMR,
                'ExAC_FIN' =>       $ExAC_FIN,
		'ExAC_NFE' =>       $ExAC_NFE,
		'ExAC_OTH' =>       $ExAC_OTH,
                'ExAC_SAS' =>       $ExAC_SAS,
		'ESP6500siv2_ALL' =>$ESP6500siv2_ALL,
		'ESP6500siv2_AA'  =>$ESP6500siv2_AA,
		'ESP6500siv2_EA'  =>$ESP6500siv2_EA,
		'CG46'            =>$CG46,
		'COSMIC70'        =>$cosmic70,
		'NCI60'           =>$nci60,
		'snp138NonFlagged'=>$snp138NonFlagged, 
		'ljb'             =>$ljb,
		'SIFT_score'      =>$SIFT_score, 
		SIFT_pred         =>$SIFT_pred, 
		pp2_HDIV          =>$Polyphen2_HDIV_score, 
		pp2_HDIV_pred     =>$Polyphen2_HDIV_pred, 
		pp2_HVAR          =>$Polyphen2_HVAR_score, 
		pp2_HVAR_pred     =>$Polyphen2_HVAR_pred, 
		LRT               =>$LRT_score,
		LRT_pred          =>$LRT_pred, 
		MutationTaster    =>$MutationTaster_score,
		MutationTaster_pred=>$MutationTaster_pred, 
		MutationAssessor  =>$MutationAssessor_score,
		MutationAssessor_pred=>$MutationAssessor_pred, 
		FATHMM            =>$FATHMM_score,
		FATHMM_pred       =>$FATHMM_pred, 
		RadialSVM         =>$RadialSVM_score,
		RadialSVM_pred    =>$RadialSVM_pred, 
		LR                =>$LR_score, 
		LR_pred           =>$LR_pred,
		VEST3_score       =>$VEST3_score,
		CADD_raw          =>$CADD_raw,
		CADD_phred        =>$CADD_phred, 
		GERP_RS           =>$GERP_RS, 
		phyloP46way_placental=>$phyloP46way_placental,
		phyloP100way_vertebrate=>$phyloP100way_vertebrate,
		SiPhy29way_logOdds=>$SiPhy29way_logOdds,
	);
	
	$data{"$chr:$pos:$mutation_id"} ={%annovar, %{$data{"$chr:$pos:$mutation_id"}}} if(exists $data{"$chr:$pos:$mutation_id"});
}


sub dataToTable {
	my %h = @_;
	my $ar_col = $h{'col'};	
	my $hr_data = $h{'data'};
	my $ar_key = $h{'key'};
	my $sep = $h{'sep'} || "\t";

	return undef unless($ar_col && $hr_data);

	my @cols = @$ar_col;
	my @keys =$ar_key ? @{$ar_key} : sort keys %$hr_data;
	my $table = "";

	unless($opt_N){	
		map {$table .= $_ . $sep} @cols; 
		$table =~s/$sep$/\n/;
	}
	foreach my $key(@keys){
		foreach my $col(@cols){
			my $val = $hr_data->{$key}{$col};
			if(defined $val){
				$val =~s/"//g;
				$val =~s/\s+$//;
				$val =~s/^\s+//;
				$val = "" if($val eq '.');
			} else{
				$val="";
			}
			$table .= $val . $sep;
		}
		$table =~s/$sep$/\n/;
	}
	
	return $table;
}


