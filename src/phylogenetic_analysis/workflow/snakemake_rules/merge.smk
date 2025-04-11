rule merge_metadata:
    message:
        "merging metadata of lab and OS samples"
    input:
        os_metadata="data/{a_or_b}/os_metadata.tsv",
        lab_annotation="input_lab_samples/{a_or_b}/comb_ann.tsv",
        lab_nextclade="input_lab_samples/{a_or_b}/nextclade_results.tsv"
    output:
        metadata="data/{a_or_b}/metadata.tsv"
    shell:
        """
        python scripts/merge_ann.py \
            --os_metadata {input.os_metadata} \
            --lab_annotation {input.lab_annotation} \
            --lab_nextclade {input.lab_nextclade} \
            --merged_metadata {output.metadata}
        """

rule merge_sequences:
    message:
        "merging sequences of lab and OS samples"
    input:
        os_sequences="data/{a_or_b}/os_sequences.fasta",
        lab_sequences="input_lab_samples/{a_or_b}/sequences.fasta"
    output:
        sequences="data/{a_or_b}/sequences.fasta"
    shell:
        """
        cat {input.os_sequences} {input.lab_sequences} > {output.sequences}
        """
