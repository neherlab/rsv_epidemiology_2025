import pandas as pd
import argparse


def parse_arguments():
    """
    Parse command-line arguments.

    Returns:
        args: Parsed arguments including input and output file paths.
    """
    parser = argparse.ArgumentParser(description="Merge and format metadata files.")
    parser.add_argument(
        "-o",
        "--os_metadata",
        required=True,
        help="Path to the OS metadata file (TSV format).",
    )
    parser.add_argument(
        "-la",
        "--lab_annotation",
        required=True,
        help="Path to the lab metadata file (TSV format).",
    )
    parser.add_argument(
        "-ln",
        "--lab_nextclade",
        required=True,
        help="Path to the nextclade output for lab samples (TSV format).",
    )
    parser.add_argument(
        "-m",
        "--merged_metadata",
        required=True,
        help="Path to save the merged metadata file (TSV format).",
    )
    return parser.parse_args()


def merge_metadata(
    os_metadata_path, lab_annotation_path, lab_nextclade_path, output_path
):
    """
    Merge OS metadata with lab metadata files, format and align them as necessary.

    Args:
        os_metadata_path (str): Path to the OS metadata file.
        lab_annotation_path (str): Path to the lab metadata file.
        lab_nextclade_path (str): Path to the nextclade output for lab samples.
        output_path (str): Path to save the merged metadata file.
    """
    # Load the OS metadata
    os_metadata = pd.read_csv(os_metadata_path, sep="\t")
    os_metadata["source"] = "NCBI"
    ann_template_series = pd.Series(
        [None] * len(os_metadata.columns), index=os_metadata.columns
    )
    ann_template_series.region = "Europe"
    ann_template_series.country = "Switzerland"
    ann_template_series.division = "Basel"
    ann_template_series.host = "Homo sapiens"
    ann_template_series.institution = "USB"
    ann_template_series.source = "USB"
    lab_ann_df = pd.read_csv(lab_annotation_path, sep="\t")
    lab_ann_df.index = lab_ann_df.NGS_ID
    nextclade_df = pd.read_csv(lab_nextclade_path, sep="\t", index_col=1)
    lab_ann_df = lab_ann_df.loc[nextclade_df.index]
    res_ann_df = pd.concat([ann_template_series] * lab_ann_df.shape[0], axis=1).T
    res_ann_df.accession = lab_ann_df.index
    res_ann_df['date'] = lab_ann_df.apply(lambda row: f"{row['Year']}-{row['Months']:02d}-XX", axis=1).values

    for col in [
        "clade",
        "qc.overallScore",
        "qc.overallStatus",
        "alignmentScore",
        "alignmentStart",
        "alignmentEnd",
    ]:
        res_ann_df[col] = nextclade_df[col].values
    res_ann_df["genome_coverage"] = nextclade_df["coverage"].values
    res_ann_df["G_coverage"] = nextclade_df["coverage"].values
    res_ann_df["F_coverage"] = nextclade_df["coverage"].values
    res_ann_df["genbank_accession_rev"] = res_ann_df["accession"]
    res_ann_df["strain"] = res_ann_df["accession"]
    for col in ["Patient_Nr",	"Year", "Months"]:
        res_ann_df[col] = lab_ann_df[col].values
    merged_metadata = pd.concat([os_metadata, res_ann_df])
    # Save the merged metadata
    merged_metadata.to_csv(output_path, sep="\t", index=False)
    print(f"Merging completed and saved to {output_path}.")

if __name__ == "__main__":
    args = parse_arguments()
    merge_metadata(args.os_metadata, args.lab_annotation, args.lab_nextclade, args.merged_metadata)
