# Copied from the monkeypox repo¹, with RSV-specific edits.
# ¹ https://github.com/nextstrain/monkeypox/blob/5c461dc7e90cd70c1f16b193f82fd1666d4c95e2/workflow/snakemake_rules/chores.smk

rule update_example_data_wildcards:
    """This updates the files under example_data/ based on latest available data from data.nextstrain.org.

    The subset of data is generated by an augur filter call which:
    - sets the subsampling size to 50
    - applies the grouping from the config
    """
    message:
        "Update example data"
    input:
        sequences = "data/{a_or_b}/sequences.fasta",
        metadata = "data/{a_or_b}/metadata.tsv",
    output:
        sequences = "example_data/{a_or_b}/sequences.fasta",
        metadata = "example_data/{a_or_b}/metadata.tsv",
    params:
        strain_id=config["strain_id_field"],
        group_by=config["filter"]["group_by"],
    shell:
        """
        augur filter \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --sequences {input.sequences} \
            --group-by {params.group_by} \
            --subsample-max-sequences 50 \
            --subsample-seed 0 \
            --output-metadata {output.metadata} \
            --output-sequences {output.sequences}
        """


rule update_example_data:
    input:
        "example_data/a/sequences.fasta",
        "example_data/b/sequences.fasta",
        "example_data/a/metadata.tsv",
        "example_data/b/metadata.tsv",
