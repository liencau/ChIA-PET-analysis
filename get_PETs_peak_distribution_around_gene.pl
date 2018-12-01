#!/usr/bin/perl -w
use strict;

my ($peak,$gene,$genebody_bin_num,$outer_bin_size,$outer_bin_num,$out)=@ARGV;
my %hash_peak;
open P,"$peak" or die;
while(<P>){
	chomp;
	my @ss=split;
	if($ss[0]=~/^[0-9]/){
		my $summit=int(($ss[1]+$ss[2])/2+0.5);
		$hash_peak{"$ss[0]-$summit"}=1;
	}
}
close P;

open G,"$gene" or die;
my %hash_key;
while(<G>){
	chomp;
	my ($chr,$stt,$end,$id,$strand)=split;
	my $len=$end-$stt+1;
	my $bin_size=$len/$genebody_bin_num;
	if($strand eq "+"){
		for(my $i=-1;$i>=-20;$i--){
			for(my $y=$stt+($i+1)*$outer_bin_size;$y>$stt+$i*$outer_bin_size;$y--){
				if(exists  $hash_peak{"$chr-$y"}){
					$hash_key{$i}+=1/$outer_bin_size;
				}else{next;}
			}
		}
		for(my $i=0;$i<$genebody_bin_num;$i++){
			for(my $y=$stt+$i*$bin_size;$y<$stt+($i+1)*$bin_size;$y++){
				if(exists  $hash_peak{"$chr-$y"}){
                                        $hash_key{$i}+=1/$bin_size;
                                }else{next;}
			}
		}
		for(my $i=0;$i<=$outer_bin_num;$i++){
                        for(my $y=$end+$i*$outer_bin_size;$y<$end+($i+1)*$outer_bin_size;$y++){
                                if(exists  $hash_peak{"$chr-$y"}){
                                        $hash_key{$i+$genebody_bin_num}+=1/$outer_bin_size;
                                }else{next;}
                        }
                }	
	}else{
		for(my $i=-1;$i>=-20;$i--){
                        for(my $y=$end-($i+1)*$outer_bin_size;$y>$end-$i*$outer_bin_size;$y--){
                                if(exists  $hash_peak{"$chr-$y"}){
                                        $hash_key{$i}+=1/$outer_bin_size;
                                }else{next;}
                        }
                }
                for(my $i=0;$i<$genebody_bin_num;$i++){
                        for(my $y=$end-$i*$bin_size;$y<$end-($i+1)*$bin_size;$y++){
                                if(exists  $hash_peak{"$chr-$y"}){
                                        $hash_key{$i}+=1/$bin_size;
                                }else{next;}
                        }
                }
                for(my $i=0;$i<=$outer_bin_num;$i++){
                        for(my $y=$stt-$i*$outer_bin_size;$y<$stt-($i+1)*$outer_bin_size;$y++){
                                if(exists  $hash_peak{"$chr-$y"}){
                                        $hash_key{$i+$genebody_bin_num}+=1/$outer_bin_size;
                                }else{next;}
                        }
                }
	}
}
close G;
my @key_sort=sort {$a<=>$b} (keys %hash_key);
open OUT,"+>$out" or die;
foreach(@key_sort){
	print OUT "$_\t$hash_key{$_}\n";
}
close OUT;
