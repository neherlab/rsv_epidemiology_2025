import os
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO, SeqRecord, Seq

# important to use a non-interactive backend, otherwise will crash on cluster
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")

from create_allele_counts import load_allele_counts


alpha = np.array(["A", "C", "G", "T", "-", "N"], dtype="U1")
gap = 200  # space between segments


def plot_coverage_concatenated(sample, ac, figure_path, primer_boundaries=None):
    coverage_stat = {}
    plt.figure(figsize=(18, 8))
    offset = 0
    ticks = []
    plt.title("sample %s" % sample)
    colnum = 0
    print("Sample", sample)
    for ref, counts in sorted(
        ac, key=lambda x: x[1].shape[-1], reverse=False
    ):  # True):
        if ref == "":
            legKey = "New"
            colour = "r"
        else:
            legKey = "Old"
            colour = "b"
        # legKey = "New" if ref == '' else "Old"
        cov = coverage(counts)
        coverage_stat[ref] = {
            "mean_cov": np.mean(cov),
            "cov > 100": (cov > 100).mean(),
            "cov > 1000": (cov > 1000).mean(),
        }
        seg = ref.split("_")[-1]
        plt.plot(
            offset + np.arange(counts.shape[-1]),
            coverage(counts),
            c=colour,
            label=legKey,
        )  #'r')
        colnum = colnum + 1  # EBH
        ticks.append([offset + counts.shape[-1] / 2, seg])
        # offset+=counts.shape[-1]+gap
        if (cov < 100).mean() > 0.20:
            print(
                sample,
                ref,
                "has very low coverage: %2.1f" % cov.mean(),
                "fraction below 100: %1.2f" % (cov < 100).mean(),
            )
        elif (cov < 1000).mean() > 0.20:
            print(
                sample,
                ref,
                "has low coverage: %2.1f" % cov.mean(),
                "fraction below 100: %1.2f" % (cov < 100).mean(),
            )
        if primer_boundaries and ref in primer_boundaries:
            for p in primer_boundaries[ref]:
                y = 2 if int(p[1]) % 2 else 6
                plt.plot(
                    [
                        primer_boundaries[ref][p]["start"],
                        primer_boundaries[ref][p]["end"],
                    ],
                    [y, y],
                    lw=10,
                    c=(0.7, 0.7, 0.7),
                )

    # plt.xticks([t[0] for t in ticks], [t[1] for t in ticks])
    plt.yscale("log")
    plt.legend(loc="lower center")
    plt.savefig(figure_path)
    plt.close()
    return coverage_stat


def get_primer_mask(primer_file, ac, omit_ref=True):
    primer_masks = {}
    for ref, counts in ac:
        print(f"ref: {ref}")
        primer_masks[ref] = np.ones(counts.shape[-1], dtype=int).astype(bool)

    if primer_file:
        primers = pd.read_csv(primer_file, skipinitialspace=True)
        for pi, p in primers.iterrows():
            if p.segment in primer_masks:
                primer_masks[p.segment][p.start : p.end] = False
            elif omit_ref:
                primer_masks[ref][p.start : p.end] = False
            else:
                print(p.segment, "is not among the mapped segments")
    return primer_masks


def coverage(ac, window=None):
    if window is None:
        return ac.sum(axis=-2)
    else:
        if np.isscalar(window):
            return np.convolve(
                ac.sum(axis=-2), np.ones(int(window)) / float(window), "same"
            )
        else:
            return np.convolve(ac.sum(axis=-2), window, "same")


def consensus(ac, min_cov=1, mask_primers=False, primer_mask=None, freq_threshold=None):
    cov = coverage(ac)
    base_frequencies = np.max(ac, axis=-2) / np.sum(ac, axis=-2)
    consensus_seq = alpha[np.argmax(ac, axis=-2)]
    consensus_seq[cov < min_cov] = "N"
    if freq_threshold:
        consensus_seq[base_frequencies < freq_threshold] = "N"
    if mask_primers and primer_mask is not None:
        consensus_seq[
            np.logical_not(primer_mask)
        ] = "N"  # bc primer mask here is 'not a primer' mask
    return consensus_seq


if __name__ == "__main__":
    # Parse input args
    parser = argparse.ArgumentParser(
        description="plot coverage, diversity and output consensus sequence",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--sample", required=True, type=str, help="the sample to analyze"
    )
    parser.add_argument(
        "--out_dir", required=True, type=str, help="directory to output"
    )
    parser.add_argument(
        "--primers", type=str, help="file with primers to mask in diversity calculation"
    )
    parser.add_argument(
        "--min_cov", type=int, default=100, help="minimal coverage to call consensus"
    )
    parser.add_argument(
        "--all_counts",
        action="store_true",
        default=False,
        help="plot coverage/diversity for all count files found",
    )
    parser.add_argument(
        "--mask_primers",
        action="store_true",
        default=False,
        help="mask primer regions with N in consensus",
    )
    parser.add_argument(
        "--freq_threshold", type=float, default=None, help="minimal base frequency to be included in the consensus"
    )

    args = parser.parse_args()
    stats = {}
    ac, ins = load_allele_counts(args.sample, allCounts=args.all_counts)
    primer_masks = get_primer_mask(args.primers, ac)

    sample = args.sample.split("/")[-1]
    os.makedirs(args.out_dir + "/figures", exist_ok=True)
    stats = plot_coverage_concatenated(
        sample, ac, args.out_dir + "/figures/coverage.png"
    )

    seqs = []
    insertions_to_include = []
    for ref, counts in ac:
        if ref != "":
            continue
        consensus_seq = consensus(
            counts,
            min_cov=args.min_cov,
            mask_primers=args.mask_primers,
            primer_mask=primer_masks[ref],
            freq_threshold=args.freq_threshold
        )
        cov = coverage(counts)
        for pos in ins[ref]:
            if cov[pos] < args.min_cov:
                continue
            total_insertion = np.sum(
                [c.sum() for _, c in list(ins[ref][pos].items())]
            )
            total_freq = 1.0 * total_insertion / cov[pos]
            max_insertion = [pos, 0, 0]
            for insertion, c in list(ins[ref][pos].items()):
                ins_freq = 1.0 * c.sum() / cov[pos]
                if ins_freq > max_insertion[2]:
                    max_insertion = [pos, insertion, ins_freq]
                if ins_freq > 0.3:
                    print(
                        sample
                        + ": frequent insertion %s at position %d with frequency %f."
                        % (insertion, pos, ins_freq)
                    )
                elif ins_freq > 0.01:
                    print(
                        sample
                        + ": rare insertion %s at position %d with frequency %f."
                        % (insertion, pos, ins_freq)
                    )

            # id the most frequent insertion is more common than no insertion
            if 1 - total_freq < max_insertion[2]:
                insertions_to_include.append(max_insertion)

        seq = "".join(consensus_seq)
        if insertions_to_include:
            complete_seq = ""
            pos = 0
            for ins_pos, ins, freq in sorted(insertions_to_include):
                complete_seq += seq[pos:ins_pos] + ins
                pos = ins_pos
                print(
                    sample
                    + ": inserted %s at position %d with frequency %f."
                    % (ins, ins_pos, freq)
                )
            complete_seq += seq[pos:]
            seq = complete_seq

        if len(ac) == 1:
            seq_name = sample
        else:
            seq_name = sample + "_" + ref
        seqs.append(
            SeqRecord.SeqRecord(
                id=seq_name, name=seq_name, description="", seq=Seq.Seq(seq)
            )
        )
    SeqIO.write(seqs, args.out_dir + "/consensus.fasta", "fasta")

    df = pd.DataFrame(stats)
    df.to_csv(args.out_dir + "/statistics.csv")
