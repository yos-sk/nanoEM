echo "Low methylation site (< 0.25)" 
echo "TP"
TP=`awk '{if ($3 < 0.25 && $6 < 0.25 && NR > 1) count+= 1} END {print count}' $1`
echo $TP
echo "FN"
FN=`awk '{if ($3 < 0.25 && $6 >= 0.25 && NR > 1) count +=1} END {print count}' $1` 
echo $FN
echo "FP"
FP=`awk '{if ($3 >= 0.25 && $6 < 0.25 && NR > 1) count += 1} END {print count}'  $1`
echo $FP
echo "TN"
TN=`awk '{if ($3 >= 0.25 && $6 >= 0.25 && NR > 1) count+=1} END {print count}' $1`
echo $TN

Precision=`echo "scale=5; $TP / ($TP + $FP)" | bc`
Recall=`echo "scale=5; $TP / ($TP + $FN)" | bc`
echo "Precision:" $Precision
echo "Recall:" $Recall

echo "Middle methylation site (>=0.25, <= 0.75)" 
echo "TP"
TP=`awk '{if ($3 >= 0.25 && $3 <= 0.75 && $6 >= 0.25 && $6 <= 0.75 && NR > 1) count+= 1} END {print count}' $1`
echo $TP
echo "FN"
FN=`awk '{if ($3 >= 0.25 && $3 <= 0.75 && ($6 < 0.25 || $6 > 0.75) && NR > 1) count +=1} END {print count}' $1` 
echo $FN
echo "FP"
FP=`awk '{if (($3 < 0.25 || $3 > 0.75) && $6 >= 0.25 && $6 <= 0.75 && NR > 1) count += 1} END {print count}'  $1`
echo $FP
echo "TN"
TN=`awk '{if (($3 < 0.25 || $3 > 0.75) && ($6 < 0.25 || $6 > 0.75) && NR > 1) count+=1} END {print count}' $1`
echo $TN

Precision=`echo "scale=5; $TP / ($TP + $FP)" | bc`
Recall=`echo "scale=5; $TP / ($TP + $FN)" | bc`
echo "Precision:" $Precision
echo "Recall:" $Recall

echo "High methylation site (> 0.75)" 
echo "TP"
TP=`awk '{if ($3 > 0.75 && $6 > 0.75 && NR > 1) count+= 1} END {print count}' $1`
echo $TP
echo "FN"
FN=`awk '{if ($3 > 0.75 && $6 <= 0.75 && NR > 1) count +=1} END {print count}' $1` 
echo $FN
echo "FP"
FP=`awk '{if ($3 <= 0.75 && $6 > 0.75 && NR > 1) count += 1} END {print count}'  $1`
echo $FP
echo "TN"
TN=`awk '{if ($3 <= 0.75 && $6 <= 0.75 && NR > 1) count+=1} END {print count}' $1`
echo $TN

Precision=`echo "scale=5; $TP / ($TP + $FP)" | bc`
Recall=`echo "scale=5; $TP / ($TP + $FN)" | bc`
echo "Precision:" $Precision
echo "Recall:" $Recall
