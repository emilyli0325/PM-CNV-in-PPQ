########################################################################################## 
# PM II/III analysis 

files=$(ls *.bam)

for file in ${files[*]}
do
bedtools genomecov -ibam $file -g chr14.genome -d | awk '$1 == "Pf3D7_14_v3" && $2 > 260000 && $2 < 327000' > $file.PMs.coverage
done


files=$(ls *.coverage)
for file in ${files[*]}
do
awk '$2>=289570 && $2<=298796{tot+=$3;cnt++}END{print tot/cnt}' $file >>1.coverage
awk '$2>=298796 || $2<=289570{tot+=$3;cnt++}END{print tot/cnt}' $file >>2.coverage
done

paste 1.coverage 2.coverage > PMs.coverage 
