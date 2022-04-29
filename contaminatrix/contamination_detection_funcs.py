import pysam
import numpy as np
import pandas as pd
import sys
from collections import defaultdict, namedtuple


def merge_sites(canonical, alt):
    """Merges a canonical primer site with an alt site, producing an interval that encompasses both
    Parameters
    ----------
    canonical : dict
        The canonical primer site, provided as a dictionary of the bed file row
    alt : dict
        The alt primer site, provided as a dictionary of the bed file row
    Returns
    -------
    dict
        A dictionary of the merged site, where the dict represents a bed file row
    """
    # base the merged site on the canonical
    mergedSite = canonical

    # check the both the canonical and alt are the same direction
    if canonical["direction"] != alt["direction"]:
        print(
            "could not merge alt with different orientation to canonical",
            file=sys.stderr,
        )
        raise SystemExit(1)

    # merge the start/ends of the alt with the canonical to get the largest window possible
    if alt["start"] < canonical["start"]:
        mergedSite["start"] = alt["start"]
    if alt["end"] > canonical["end"]:
        mergedSite["end"] = alt["end"]

    return mergedSite


def getPrimerDirection(primerID):
    """Infer the primer direction based on it's ID containing LEFT/RIGHT
    Parameters
    ----------
    primerID : string
        The primer ID from the 4th field of the primer scheme
    """
    if "LEFT" in primerID:
        return "+"
    elif "RIGHT":
        return "-"
    else:
        print("LEFT/RIGHT must be specified in Primer ID", file=sys.stderr)
        raise SystemExit(1)


def read_bed_file(fn):
    """Parses a bed file and collapses alts into canonical primer sites
    Parameters
    ----------
    fn : str
        The bedfile to parse
    Returns
    -------
    list
        A list of dictionaries, where each dictionary contains a row of the parsed bedfile.
        The available dictionary keys are - Primer_ID, direction, start, end
    """

    # read the primer scheme into a pandas dataframe and run type, length and null checks
    primers = pd.read_csv(
        fn,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "Primer_ID", "PoolName"],
        dtype={
            "chrom": str,
            "start": int,
            "end": int,
            "Primer_ID": str,
            "PoolName": str,
        },
        usecols=(0, 1, 2, 3, 4),
        skiprows=0,
    )
    if len(primers.index) < 1:
        print("primer scheme file is empty", file=sys.stderr)
        raise SystemExit(1)
    if primers.isnull().sum().sum():
        print("malformed primer scheme file", file=sys.stderr)
        raise SystemExit(1)

    # compute the direction
    primers["direction"] = primers.apply(
        lambda row: getPrimerDirection(row.Primer_ID), axis=1
    )

    # separate alt primers into a new dataframe
    altFilter = primers["Primer_ID"].str.contains("_alt")
    alts = pd.DataFrame(
        columns=("chrom", "start", "end", "Primer_ID", "PoolName", "direction")
    )
    alts = pd.concat([alts, primers[altFilter]])
    primers = primers.drop(primers[altFilter].index.values)

    # convert the primers dataframe to dictionary, indexed by Primer_ID
    #  - verify_integrity is used to prevent duplicate Primer_IDs being processed
    bedFile = primers.set_index(
        "Primer_ID", drop=False, verify_integrity=True
    ).T.to_dict()

    # if there were no alts, return the bedfile as a list of dicts
    if len(alts.index) == 0:
        return list(bedFile.values())

    # merge alts
    for _, row in alts.iterrows():
        primerID = row["Primer_ID"].split("_alt")[0]

        # check the bedFile if another version of this primer exists
        if primerID not in bedFile:

            # add to the bed file and continue
            bedFile[primerID] = row.to_dict()
            continue

        # otherwise, we've got a primer ID we've already seen so merge the alt
        mergedSite = merge_sites(bedFile[primerID], row)

        # update the bedFile
        bedFile[primerID] = mergedSite

    # return the bedFile as a list
    return [value for value in bedFile.values()]


def generate_amplicons(bedfile):
    bed = read_bed_file(bedfile)
    primer_pairs = defaultdict(dict)
    for b in bed:
        scheme, pair, side = b["Primer_ID"].split("_")
        primer_pairs[pair][side] = b
    # this data structure is more useful for searching...
    amplicons = np.fromiter(
        (
            (
                k,
                v["LEFT"]["PoolName"],
                v["LEFT"]["end"],
                v["RIGHT"]["start"],  # just insert
                v["LEFT"]["start"],
                v["RIGHT"]["end"],  # contains primers
                v["LEFT"]["Primer_ID"],
                v["RIGHT"]["Primer_ID"],
            )
            for k, v in primer_pairs.items()
        ),
        dtype=[
            ("name", int),
            ("pool", int),
            ("insert_start", int),
            ("insert_end", int),
            ("start", int),
            ("end", int),
            ("left_primer", "U20"),
            ("right_primer", "U20"),
        ],
    )
    return amplicons


def generate_rg_overlaps(bedfile):
    amplicons = generate_amplicons(bedfile)
    pool_1 = [amplicon for amplicon in amplicons if amplicon["pool"] == 1]
    pool_2 = [amplicon for amplicon in amplicons if amplicon["pool"] == 2]
    overlap_coordinates = []
    for pool_1_amp in pool_1:
        for pool_2_amp in pool_2:
            overlap = range(
                max(pool_1_amp["insert_start"], pool_2_amp["insert_start"]) + 1,
                min(pool_1_amp["insert_end"], pool_2_amp["insert_end"]) - 1,
            )
            if len(overlap) != 0:
                overlap_coordinates.append((overlap.start, overlap.stop))
    return overlap_coordinates


def find_sus_positions(args):
    overlap_tuples = generate_rg_overlaps(args.bedfile)
    suspicious_positions = {}

    with pysam.AlignmentFile(
        args.bamfile, "rb", reference_filename=args.reference
    ) as bam:
        if len(bam.references) > 1:
            print("BAM file appears to be aligned to more than one reference")
            sys.exit(2)

        ref_name = bam.references[0]

        for overlap_start, overlap_end in overlap_tuples:
            for pileupcolumn in bam.pileup(
                ref_name, start=overlap_start, end=overlap_end
            ):
                if pileupcolumn.reference_pos + 1 not in range(
                    overlap_start, overlap_end
                ):
                    continue
                if not pileupcolumn.n >= args.min_coverage:
                    continue
                counts = {
                    1: {"A": 0, "C": 0, "G": 0, "T": 0},
                    2: {"A": 0, "C": 0, "G": 0, "T": 0},
                }
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:

                        RG = int(pileupread.alignment.get_tag("RG"))
                        if RG not in (1, 2):
                            print("Read skipped due to RG mismatch")
                            continue
                        base = pileupread.alignment.query_sequence[
                            pileupread.query_position
                        ]
                        counts[RG][base] += 1
                freqs = {
                    1: {
                        "total": sum(counts[1].values()),
                        "freqs": {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0},
                    },
                    2: {
                        "total": sum(counts[2].values()),
                        "freqs": {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0},
                    },
                }
                if all(count == 0 for count in counts[1].values()) or all(
                    (count == 0 for count in counts[2].values())
                ):
                    continue
                for read_group in (1, 2):
                    for base in ("A", "C", "G", "T"):
                        freqs[read_group]["freqs"][base] = (
                            counts[read_group][base] / freqs[read_group]["total"]
                        )
                if any(
                    abs(freq_1 - freq_2) >= args.min_difference
                    for freq_1, freq_2 in zip(
                        freqs[1]["freqs"].values(), freqs[2]["freqs"].values()
                    )
                ):
                    suspicious_positions[pileupcolumn.reference_pos + 1] = freqs
    return suspicious_positions


def print_positions(args, sus_positions):
    if not args.detailed_report:
        print("\n".join(str(x) for x in sus_positions.keys()), file=sys.stdout)
    else:
        print("Pos\tRG\trg_depth\tA\tC\tG\tT", file=sys.stdout)
        for position, details in sus_positions.items():
            print(
                "\t".join(
                    str(x)
                    for x in [
                        position,
                        1,
                        details[1]["total"],
                        details[1]["freqs"]["A"],
                        details[1]["freqs"]["C"],
                        details[1]["freqs"]["G"],
                        details[1]["freqs"]["T"],
                    ]
                ),
                file=sys.stdout,
            )
            print(
                "\t".join(
                    str(x)
                    for x in [
                        position,
                        2,
                        details[2]["total"],
                        details[2]["freqs"]["A"],
                        details[2]["freqs"]["C"],
                        details[2]["freqs"]["G"],
                        details[2]["freqs"]["T"],
                    ]
                ),
                file=sys.stdout,
            )

