process FILTER {
    tag "${chrom}:${category}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/filtered/", mode: 'copy')

    input:
    tuple val(pheno), val(chrom),
          path(file), path(index),
          val(category)

    output:
    tuple val(pheno), val(chrom), val(category),
          path("${pheno}.${chrom}.${category}.vcf.gz"),
          path("${pheno}.${chrom}.${category}.vcf.gz.tbi"),
          path("${pheno}.${chrom}.${category}.tsv")
        
    script:
    if ( category == 'Pathogenic' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,gnomADe_AF:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'MAF > ${params.MAF} || HWE < ${params.HWE} || ExcHet < ${params.ExcHet}' | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'CLIN_SIG ~ "pathogenic" || CLIN_SIG ~ "likely_pathogenic"' | \
        bcftools view -e 'gnomADe_AF > ${params.PAF} || MAX_AF > ${params.PAF}' | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${pheno}.${chrom}.${category}.vcf.gz

        tabix ${pheno}.${chrom}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${pheno}.${chrom}.${category}.vcf.gz > ${pheno}.${chrom}.${category}.tsv
        """
    } else if ( category == 'Damaging' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,IMPACT,CADD_PHRED:Float,gnomADe_AF:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'MAF > ${params.MAF} || HWE < ${params.HWE} || ExcHet < ${params.ExcHet}' | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'IMPACT="HIGH"' | \
        bcftools view -i 'CADD_PHRED > ${params.CADD}' | \
        bcftools view -e 'gnomADe_AF > ${params.gnomADe_AF} || MAX_AF > ${params.MAX_AF}' | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${pheno}.${chrom}.${category}.vcf.gz

        tabix ${pheno}.${chrom}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${pheno}.${chrom}.${category}.vcf.gz > ${pheno}.${chrom}.${category}.tsv
        """
    } else if ( category == 'High' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,IMPACT,CADD_PHRED:Float,gnomADe_AF:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'MAF > ${params.MAF} || HWE < ${params.HWE} || ExcHet < ${params.ExcHet}' | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'IMPACT="HIGH"' | \
        bcftools view -e 'gnomADe_AF > ${params.gnomADe_AF} || MAX_AF > ${params.MAX_AF}' | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${pheno}.${chrom}.${category}.vcf.gz

        tabix ${pheno}.${chrom}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${pheno}.${chrom}.${category}.vcf.gz > ${pheno}.${chrom}.${category}.tsv
        """
    } else if ( category == 'Stop' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,IMPACT,Consequence,CADD_PHRED:Float,gnomADe_AF:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'MAF > ${params.MAF} || HWE < ${params.HWE} || ExcHet < ${params.ExcHet}' | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'Consequence~"stop_gained"' | \
        bcftools view -e 'gnomADe_AF > ${params.gnomADe_AF} || MAX_AF > ${params.MAX_AF}' | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${pheno}.${chrom}.${category}.vcf.gz

        tabix ${pheno}.${chrom}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${pheno}.${chrom}.${category}.vcf.gz > ${pheno}.${chrom}.${category}.tsv
        """
    } else if ( category == 'PTV' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,IMPACT,Consequence,CADD_PHRED:Float,gnomADe_AF:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'MAF > ${params.MAF} || HWE < ${params.HWE} || ExcHet < ${params.ExcHet}' | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'Consequence~"stop_gained" || Consequence~"frameshift_variant" || Consequence~"splice_acceptor_variant"' | \
        bcftools view -e 'gnomADe_AF > ${params.gnomADe_AF} || MAX_AF > ${params.MAX_AF}' | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${pheno}.${chrom}.${category}.vcf.gz

        tabix ${pheno}.${chrom}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${pheno}.${chrom}.${category}.vcf.gz > ${pheno}.${chrom}.${category}.tsv
        """
    } else if ( category == 'Scored' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,IMPACT,CADD_PHRED:Float,gnomADe_AF:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'MAF > ${params.MAF} || HWE < ${params.HWE} || ExcHet < ${params.ExcHet}' | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'IMPACT="HIGH" || IMPACT="MODERATE"' | \
        bcftools view -i 'CADD_PHRED > ${params.CADD}' | \
        bcftools view -e 'gnomADe_AF > ${params.gnomADe_AF} || MAX_AF > ${params.MAX_AF}' | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${pheno}.${chrom}.${category}.vcf.gz

        tabix ${pheno}.${chrom}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${pheno}.${chrom}.${category}.vcf.gz > ${pheno}.${chrom}.${category}.tsv
        """
    } else if ( category == 'NotScored' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,IMPACT,CADD_PHRED:Float,gnomADe_AF:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'MAF > ${params.MAF} || HWE < ${params.HWE} || ExcHet < ${params.ExcHet}' | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'IMPACT="HIGH" || IMPACT="MODERATE"' | \
        bcftools view -e 'gnomADe_AF > ${params.gnomADe_AF} || MAX_AF > ${params.MAX_AF}' | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${pheno}.${chrom}.${category}.vcf.gz

        tabix ${pheno}.${chrom}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${pheno}.${chrom}.${category}.vcf.gz > ${pheno}.${chrom}.${category}.tsv
        """
    } else if ( category == 'Splicing' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,SpliceAI_pred_DS_AG:Float,SpliceAI_pred_DS_AL:Float,SpliceAI_pred_DS_DG:Float,SpliceAI_pred_DS_DL:Float,gnomADe_AF:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'MAF > ${params.MAF} || HWE < ${params.HWE} || ExcHet < ${params.ExcHet}' | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'SpliceAI_pred_DS_AG > ${params.DS} || SpliceAI_pred_DS_AL > ${params.DS} || SpliceAI_pred_DS_DG > ${params.DS} || SpliceAI_pred_DS_DL > ${params.DS}' | \
        bcftools view -e 'gnomADe_AF > ${params.gnomADe_AF} || MAX_AF > ${params.MAX_AF}' | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${pheno}.${chrom}.${category}.vcf.gz

        tabix ${pheno}.${chrom}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${pheno}.${chrom}.${category}.vcf.gz > ${pheno}.${chrom}.${category}.tsv
        """
    }
}
