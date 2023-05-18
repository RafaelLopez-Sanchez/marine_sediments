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
## Notes: This script extracts the CAZymes with a signal Peptide found in our dbCAN2 output with stringent cuttoff and at least two of the dbCAN 2 tools being mandatory to be HMMER one of them. 
##   
##
## ---------------------------
#########THEME##########################################################
#####dbCAN annotation of extracellular CAZymes###########

###Use dbCAN2.

$ run_dbcan.py *.fna -c cluster --cgc_dis 2 --cgc_sig_genes all --db_dir /tres/DB/dbCAN2 prok --use_signalP all --out_dir output_dbcan2

############################# CAZYME COUNTS FOR EACH GENOME####################################################


#1.Use cutoff for bacteria in the hmmer file file, use E-value < 1e-18 and coverage > 0.35

# 1. Run: cat hmmer.out| awk '$5<1e-18&&$10>0.35' > hmmer.out.stringent (this allows you to get the same result as what is produced in our dbCAN2 webpage). Se which results got a signal peptide.
cat hmmer.out | awk '$5<1e-18&&$10>0.35' > hmmer.out.stringent
grep "Y" overview.txt | cut -f1 > extra.list

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
# 1. Filter data.
#1.1 Filter only hmms that were validate by another of the dbCAN2 tool with the cutoff_list.txt and see if the have the signal peptide with the extra.list. Put results into the extra.hmmer.out.stringent.cutoff file
awk '$6 > 1 { print }' overview.txt |cut -f1 > cutoff_list.txt 
cat hmmer.out.stringent | grep -f extra.list |  grep -f cutoff_list.txt > extra.hmmer.out.stringent.cutoff
#1.2
cut -f1,3 extra.hmmer.out.stringent.cutoff  | sort -n > proteins.stringent.txt
# 1.3
sed -i 's/\.hmm//' proteins.stringent.txt
# 1.4
cut -f1 proteins.stringent.txt > proteins.only.stringent.txt
# 1.5
uniq -c proteins.only.stringent.txt > proteins.stringent.count
# 1.6 
cat *.only.stringent.txt > all.only.stringent.txt
# 1.7
sort all.only.stringent.txt |  uniq -c  > all.uniq.stringent.count
# 1.8
sort all.only.stringent.txt |  uniq   > all.stringent.uniq
# 1.9
for s in $(cat all.stringent.uniq); do (grep -c -w $s proteins.only.stringent.txt >> proteins.stringent.full_count);done
# 1.10
paste all.stringent.uniq proteins.stringent.full_count > extra.cazy_counts.stringent.txt

rm hmmer.out.stringent cutoff_list.txt  proteins.stringent.txt proteins.only.stringent.txt proteins.stringent.count all.only.stringent.txt all.uniq.stringent.count all.stringent.uniq proteins.stringent.full_count extra.hmmer.out.stringent.cutoff extra.list