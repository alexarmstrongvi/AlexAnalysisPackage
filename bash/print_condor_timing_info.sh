OUTPUT_FILE="info_output_timing.txt"
SEARCH_STR=""

echo "Information on the largest/slowest samples" | tee $OUTPUT_FILE
echo -e "\n\nRanked by job time\n\n" | tee -a $OUTPUT_FILE
grep -oP "(?<=Analysis time: Real ).*(?=, CPU)" ./*$SEARCH_STR*out | sort -t ":" -k 2,4 -r | cat -n >> $OUTPUT_FILE
echo -e "\n\nRanked by number of selected events\n\n" | tee -a $OUTPUT_FILE
grep -oP -m 1 "(?<=Cut 10: Tau veto                        : )[0-9]*" ./*$SEARCH_STR*out | sort -t ":" -k 2 -rg | cat -n >> $OUTPUT_FILE
echo -e "\n\nRanked by number of processed events\n\n" | tee -a $OUTPUT_FILE
grep -oP "(?<=Number of events processed: ).*" ./*$SEARCH_STR*out | sort -t ":" -k 2 -rg | cat -n >> $OUTPUT_FILE
echo -e "\n\nRanked by event processing speed (Events per millisecond)\n\n" | tee -a $OUTPUT_FILE
grep -oP "(?<=Analysis speed \[kHz\]: ).*" ./*$SEARCH_STR*out | sort -t ":" -k 2 -rg | cat -n >> $OUTPUT_FILE
echo "DONE"
