##
##@param1: file_list name, each row is one sample following: BEDGRAPH file (from biscuit) and sample_name
##@param2: coverage (suggestion: 5)
##@param3: group1 name
##@param4: group2 name
##@output1: metilene_$group1_$group2.sorted.out - default metilene output
##@output2: metilene_$group1_$group2.sorted.q.out - sorted output by q-values


file_list=$1
cov=$2
group1=$3
group2=$4
list1=""
list2=""
while read filename group
do
        echo $filename
        echo perl /common/bermanblab/bin/wgbs/qsubexc.pl "sort -k1,1 -k2,2n $filename | cut -f 1,2,3,4,5 >$filename".sorted"" bedsort$group
        awk '$5 >= 5'  $filename.sorted | cut -f1,2,3,4  > $filename.sorted.cov$cov
        if [ "$group" == $group1 ]; then
                list1+=$filename".sorted.cov$cov,"
        else
                list2+=$filename".sorted.cov$cov,"
        fi
done < $file_list
list1=${list1::-1}
list2=${list2::-1}
perl /common/bermanblab/bin/wgbs/metilene_v0.2-6/metilene_input.pl --in1 $list1 --in2 $list2 --h1 $group1 --h2 $group2  -b /common/bermanblab/bin/bedtools
/common/bermanblab/bin/wgbs/metilene_v0.2-6/metilene -a N -b P metilene_N_P.input | sort -V -k1,1 -k2,2n > metilene_$group1_$group2.sorted.out
perl /common/bermanblab/bin/wgbs/metilene_v0.2-6/metilene_output.pl -q metilene_$group1_$group2.sorted.out -o metilene_$group1_$group2.sorted.q.out
