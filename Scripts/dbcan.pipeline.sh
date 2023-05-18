## ---------------------------
##
##
##
##
## Author: Rafael López-Sánchez
##
##
## Copyright (c) Rafael, López-Sánchez 2023
## Email: rafael.lopez@ibt.unam.mx
##
## ---------------------------
##
## Notes: This script extracts the CAZymes found in our dbCAN2 output with stringent cuttoff and at least two of the dbCAN 2 tools being mandatory to be HMMER one of them. 
##   
##
## ---------------------------
#########THEME##########################################################
#####dbCAN annotation of  CAZymes###########

###Use dbCAN2.

$ run_dbcan.py *.fna -c cluster --cgc_dis 2 --cgc_sig_genes all --db_dir /tres/DB/dbCAN2 prok --use_signalP all --out_dir output_dbcan2
############################# CAZYME COUNTS FOR EACH GENOME####################################################
#.Use cutoff for bacteria in the hmmer file file, use E-value < 1e-18 and coverage > 0.35
#Run: cat hmmer.out| awk '$5<1e-18&&$10>0.35' > hmmer.out.stringent (this allows you to get the same result as what is produced in our dbCAN2 webpage).
 cat hmmer.out | awk '$5<1e-18&&$10>0.35' > hmmer.out.stringent
#Cols in hmmer.out and hmmer.out.stringent:
#1. Family HMM
#2. HMM length
#3. Query ID
#4. Query length
#5. E-value (how similar to the family HMM)
#6. HMM start
#7. HMM end
#8. Query start
#9. Query end
#10. Coverage
# 7. Filter data.
#7.1 Filter only hmms that were validate by another of the dbCAN2 tool with the cutoff_list.txt. Put results into the hmmer.out.stringent.cutoff file
awk '$6 > 1 { print }' overview.txt |cut -f1 > cutoff_list.txt
cat hmmer.out.stringent | grep -f cutoff_list.txt > hmmer.out.stringent.cutoff
#7.2
cut -f1,3 hmmer.out.stringent.cutoff  | sort -n > proteins.stringent.txt
# 7.3
 sed -i 's/\.hmm//' proteins.stringent.txt
# 7.4
 cut -f1 proteins.stringent.txt > proteins.only.stringent.txt
# 7.5
 uniq -c proteins.only.stringent.txt > proteins.stringent.count
# 7.6 
 cat *.only.stringent.txt > all.only.stringent.txt
# 7.7
 sort all.only.stringent.txt |  uniq -c  > all.uniq.stringent.count
# 7.8
sort all.only.stringent.txt |  uniq   > all.stringent.uniq
#7.9
 for s in $(cat all.stringent.uniq); do (grep -c -w $s proteins.only.stringent.txt >> proteins.stringent.full_count);done
# 7.10
paste all.stringent.uniq proteins.stringent.full_count > cazy_counts.stringent.txt

rm hmmer.out.stringent cutoff_list.txt hmmer.out.stringent.cutoff proteins.stringent.txt proteins.only.stringent.txt proteins.stringent.count all.only.stringent.txt all.uniq.stringent.count all.stringent.uniq proteins.stringent.full_count