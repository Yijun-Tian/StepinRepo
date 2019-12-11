#This script can be used to detect split circles from NGS data
# $1 is Read1 file, $2 is Read2 file, $3 is sample prefix
bowtie2 --threads 20 --very-sensitive-local -x ensembl_hg38_clean -U $1,$2 -S $3.primaryalign.sam
extractSoftclipped -l 20 $3.primaryalign.sam > $3.softclip.fastq.gz
bowtie2 --threads 15 --very-sensitive-local -x ensembl_hg38_clean -U $3.softclip.fastq.gz -S $3.softclip.sam
rm $3.SplitReads.txt
rm $3.MixtureInput
samtools view $3.softclip.sam | sed 's/|/\t/g' | awk '{if($2==$7 && $3==$6) print $0;}' | awk '$5~/^[2-9][0-9]S[0-9]{1,3}M$/ {print $0,"Left";}'  > $3.LeftSInput
samtools view $3.softclip.sam | sed 's/|/\t/g' | awk '{if($2==$7 && $3==$6) print $0;}' | awk '$5~/^[0-9]{1,3}M[2-9][0-9]S$/ {print $0,"Right";}'  > $3.RightSInput
cat $3.LeftSInput $3.RightSInput | awk '$10~/^[0-9]{1,3}M$/ {print $0;}'> $3.MixtureInput
while IFS= read -r line
do
side=`echo $line | awk '{print \$NF}'`
chr=`echo $line | awk '{print \$2}'`
Strand=`echo $line | awk '{print \$3}'`
FirstPos=`echo $line | awk '{print \$4}'`
FirstCig=`echo $line | awk '{print \$5}'`
SeconPos=`echo $line | awk '{print \$8}'`
SeconCig=`echo $line | awk '{print \$10}'`
if ([ "$side" == "Left" ] && [ $Strand = 0 ] && [ $SeconPos -gt $FirstPos ]) || ([ "$side" == "Right" ] && [ $Strand = 16 ] && [ $SeconPos -gt $FirstPos ])
then
     Breakpoint1=$FirstPos
     Complement=`echo $SeconCig | sed -s 's/M//g'`
     Breakpoint2=$((SeconPos+Complement))
     echo $chr:${Breakpoint1}-${Breakpoint2} >> $3.SplitReads.txt
     echo "We got 1 split read!!!"
elif ([ "$side" == "Left" ] && [ $Strand = 16 ] && [ $FirstPos -gt $SeconPos ]) || ([ "$side" == "Right" ] && [ $Strand = 0 ] && [ $FirstPos -gt $SeconPos ])
then
     Breakpoint1=$SeconPos
     Complement=`echo $FirstCig | sed -e 's/[2-9][0-9]S//' | sed -s 's/M//g'`
     Breakpoint2=$((FirstPos+Complement))
     echo $chr:${Breakpoint1}-${Breakpoint2} >> $3.SplitReads.txt
     echo "We got 1 split read!!!"
else
     echo "We got 1 mystery read..."
fi
done < "$3.MixtureInput"

for i in $( cat $3.SplitReads.txt | sort | uniq )
do
(
 echo $i >> $3.coordinate.txt
 grep $i $3.SplitReads.txt -w -c >> $3.count.txt
) &
done
wait
paste $3.coordinate.txt $3.count.txt > $3.split.profile





























