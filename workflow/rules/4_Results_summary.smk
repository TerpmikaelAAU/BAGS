rule summarize_genome:
    input:
        busco="data/busco/{genome}/short_summary.txt",
        antismash="data/Antismash/{genome}",
    output:
        summary="data/summary/{genome}_summary.tsv",
    threads: 1
    resources:
        mem_mb=1000,
        time="00:10:00",
    shell:
        """
        # 1. Extract the BUSCO score string (e.g., C:98.5%[S:98.5%,D:0.0%],F:0.5%,M:1.0%,n:4046)
        BUSCO_SCORE=$(grep -o "C:.*%" {input.busco} | head -n 1 || echo "N/A")
        
        # 2. Count the number of individual region .gbk files antiSMASH generated
        TOTAL_BGCs=$(ls {input.antismash}/*region*.gbk 2>/dev/null | wc -l || echo 0)
        
        # 3. Search those region files for specific cluster type signatures.
        # We use grep -il to count the matching files, not the matching lines.
        # The '|| echo 0' ensures the pipeline doesn't crash if a genome has 0 BGCs.
        NRPS_COUNT=$(grep -il "NRPS" {input.antismash}/*region*.gbk 2>/dev/null | wc -l || echo 0)
        PKS_COUNT=$(grep -il "PKS" {input.antismash}/*region*.gbk 2>/dev/null | wc -l || echo 0)
        TERPENE_COUNT=$(grep -il "terpene" {input.antismash}/*region*.gbk 2>/dev/null | wc -l || echo 0)
        
        # 4. Write everything to a tab-separated format
        echo -e "Genome\tBUSCO_Score\tTotal_BGCs\tNRPS_Regions\tPKS_Regions\tTerpene_Regions" > {output.summary}
        echo -e "{wildcards.genome}\t$BUSCO_SCORE\t$TOTAL_BGCs\t$NRPS_COUNT\t$PKS_COUNT\t$TERPENE_COUNT" >> {output.summary}
        """

rule compile_summaries:
    input:
        # Calls the dynamic function to get all the individual summary TSV files
        summaries=aggregate_results,
    output:
        final_report="data/final_BAGS_report.tsv",
    threads: 1
    resources:
        mem_mb=1000,
        time="00:10:00",
    shell:
        """
        # Copy the header from the very first summary file
        head -n 1 {input.summaries[0]} > {output.final_report}
        
        # Append the data rows (skipping the header line) from all summary files
        for f in {input.summaries}; do
            tail -n +2 $f >> {output.final_report}
        done
        """