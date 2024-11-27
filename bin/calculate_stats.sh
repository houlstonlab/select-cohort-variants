#!/bin/bash
pheno=$1
category=$2
file=$3

# TODO: write functions to calculate the rest of the statistics
# Function to query bcftools and get the maximum value
get_val() {
    local query=$1
    local file=$2
    local which=$3

    bcftools query -f "$query" $file | \
        awk '{print $1 + 0}' | \
        grep -v '^0$' | \
        sort -ug > tmp.txt
    
    if [ "$which" == "max" ]; then
        tail -1 tmp.txt
    elif [ "$which" == "min" ]; then
        head -1 tmp.txt
    else
        echo "Invalid option. Please choose either 'max' or 'min'."
        exit 1
    fi
    rm tmp.txt
}

# Genotype: Missing, GQ, DP, VAF
N_MSS=$( bcftools query -f '[ %GT]\n' ${file} | tr ' ' '\n' | sort -u | grep '\\.' | wc -l)
MIN_GQ=$( bcftools query -f '[ %GT=%GQ]\n' ${file}  | tr ' ' '\n' | grep -v '^$' | grep -v '0/0' | cut -d '=' -f 2 | sort -ug | head -1)
MIN_DP=$( bcftools query -f '[ %GT=%DP]\n' ${file}  | tr ' ' '\n' | grep -v '^$' | grep -v '0/0' | cut -d '=' -f 2 | sort -ug | head -1)
MIN_VF=$( bcftools query -f '[ %GT=%VAF]\n' ${file}  | tr ' ' '\n' | grep -v '^$' | grep -v '0/0' | cut -d '=' -f 2 | sort -ug | head -1)

# Variants: HWE, ExcHet, MAF
MIN_HWE=$(get_val "%HWE\n" ${file} "min")
MIN_HET=$(get_val "%ExcHet\n" ${file} "min")
MAX_MAF=$(get_val "%MAF\n" ${file} "max")

# MIN_HWE=$( bcftools query -f '%HWE\n' ${file} | sort -ug | head -1 )
# MIN_HET=$( bcftools query -f '%ExcHet\n' ${file} | sort -ug | head -1 )
# MAX_MAF=$( bcftools query -f '%MAF\n' ${file} | sort -ug | tail -1 )

# Annotations: gnomADe_AF, MAX_AF, CLIN_SIG, IMPACT, SpliceAI_pred_DS_AG, SpliceAI_pred_DS_AL, SpliceAI_pred_DS_DG, SpliceAI_pred_DS_DL
MIN_AF=$( bcftools +split-vep -s worst -c gnomADe_AF:Float -f '%gnomADe_AF\n' ${file} | sort -ug | tail -1 )
MIN_MAX=$( bcftools +split-vep -s worst -c MAX_AF:Float -f '%MAX_AF\n' ${file} | sort -ug | tail -1 )
MIN_CADD=$( bcftools +split-vep -s worst -c CADD_PHRED:Float -f '%CADD_PHRED\n' ${file} | sort -ug | head -1 )
MIN_DS=$( bcftools +split-vep -s worst -c SpliceAI_pred_DS_AG:Float,SpliceAI_pred_DS_AL:Float,SpliceAI_pred_DS_DG:Float,SpliceAI_pred_DS_DL:Float -f '%SpliceAI_pred_DS_AG %SpliceAI_pred_DS_AL %SpliceAI_pred_DS_DG %SpliceAI_pred_DS_DL\n' ${file} | awk '{ max = $1; for (i = 2; i <= NF; i++) {if ($i > max) { max = $i }}; print max }' | sort -ug | head -1 )
COUNT_CLN=$( bcftools +split-vep -s worst -c CLIN_SIG -f '%CLIN_SIG\n' ${file} | sort -u | grep -v pathogenic | wc -l )
COUNT_VEP=$( bcftools +split-vep -s worst -c IMPACT -f '%IMPACT\n' ${file} | sort -u | grep -v HIGH | wc -l )

# Print headers
echo -e "Pheno\tChrom\tCategory\tN_MSS\tMIN_GQ\tMIN_DP\tMIN_VF\tMIN_HWE\tMIN_HET\tMAX_MAF\tMIN_AF\tMIN_MAX\tMIN_DS\tCOUNT_CLN\tCOUNT_VEP"

# Print values
echo -e "$pheno\t$chrom\t$category\t$N_MSS\t$MIN_GQ\t$MIN_DP\t$MIN_VF\t$MIN_HWE\t$MIN_HET\t$MAX_MAF\t$MIN_AF\t$MIN_MAX\t$MIN_DS\t$COUNT_CLN\t$COUNT_VEP"