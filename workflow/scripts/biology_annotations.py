"""
Shared biology annotations for repeat elements.

Pure helper functions used by reference-building scripts and downstream DGE
reports to stratify detected repeats by biology-driven properties that help
explain detection vs non-detection in polyA-selection assays.

Functions:

  classify_family_age(rep_name)
      Approximate age bracket (young / medium / old / unknown) from a
      RepeatMasker subfamily name. Useful at the gene_id granularity of the
      pipeline where gene_id == repName.

  classify_family_age_from_divergence(milli_div)
      Age bracket from RepeatMasker milliDiv (per-mille base divergence).
      More principled than the name-based version when milliDiv is available.

  strand_vs_gene(repeat_strand, gene_strands)
      Classify a repeat locus as sense / antisense / both / intergenic with
      respect to overlapping gene annotation on its own and opposite strand.

  longest_a_run(sequence)
      Length of the longest contiguous run of adenines in a sequence. Used
      as a proxy for internal polyA priming potential when passed a region
      downstream of a repeat in mRNA-sense orientation.

All functions are pure and side-effect free so they can be reused from
multiple scripts, imported in tests, and called row-by-row from streaming
readers over large rmsk/gtf inputs.
"""

import re


_L1_YOUNG = re.compile(r'^L1HS$|^L1PA[2-7](?:\D|$)')
_L1_MEDIUM = re.compile(r'^L1PA(?:[89]|1[0-6])(?:\D|$)')
_L1_OLD = re.compile(r'^L1PA1[7-9]|^L1PA[2-9]\d|^L1PB|^L1M[A-F]')

_ALU_YOUNG = re.compile(r'^AluY')
_ALU_MEDIUM = re.compile(r'^AluS')
_ALU_OLD = re.compile(r'^AluJ|^FLAM|^FRAM')

_SVA_YOUNG = re.compile(r'^SVA_[A-F]$')

_HERV_YOUNG = re.compile(r'^HERVK|^LTR5|^HERVH')


def classify_family_age(rep_name):
    """Return an age bracket for a repeat subfamily name.

    Returns one of "young", "medium", "old", or "unknown".

    Brackets are approximate and based on published subfamily ages for the
    major retrotransposon classes:

        L1HS, L1PA2-7        young    (roughly below 15 My)
        L1PA8-16             medium
        L1PA17+, L1PB, L1M*  old      (roughly above 50 My)
        AluY                 young
        AluS                 medium
        AluJ, FLAM, FRAM     old
        SVA_A-F              young    (below 25 My)
        HERVK, LTR5, HERVH   young    (loose call for ERV proviruses)

    Anything not matching returns "unknown" so callers can decide how to
    treat it.
    """
    if not rep_name:
        return "unknown"
    if _L1_YOUNG.match(rep_name):
        return "young"
    if _L1_MEDIUM.match(rep_name):
        return "medium"
    if _L1_OLD.match(rep_name):
        return "old"
    if _ALU_YOUNG.match(rep_name):
        return "young"
    if _ALU_MEDIUM.match(rep_name):
        return "medium"
    if _ALU_OLD.match(rep_name):
        return "old"
    if _SVA_YOUNG.match(rep_name):
        return "young"
    if _HERV_YOUNG.match(rep_name):
        return "young"
    return "unknown"


YOUNG_MILLIDIV_MAX = 50
MEDIUM_MILLIDIV_MAX = 200


def classify_family_age_from_divergence(milli_div):
    """Bucket a RepeatMasker milliDiv value into an age bracket.

    Returns one of "young", "medium", "old", or "unknown". milliDiv is
    expressed in per-mille substitutions (0 to ~1000). Thresholds map
    roughly to the subfamily brackets used in classify_family_age:

        milliDiv <= 50    young
        50 < milliDiv <= 200   medium
        milliDiv > 200    old

    Non-finite and negative inputs return "unknown".
    """
    if milli_div is None:
        return "unknown"
    try:
        d = float(milli_div)
    except (TypeError, ValueError):
        return "unknown"
    if d != d or d < 0:
        return "unknown"
    if d <= YOUNG_MILLIDIV_MAX:
        return "young"
    if d <= MEDIUM_MILLIDIV_MAX:
        return "medium"
    return "old"


def strand_vs_gene(repeat_strand, gene_strands):
    """Classify a repeat locus's strand relationship to overlapping genes.

    Parameters
    ----------
    repeat_strand : str
        Strand of the repeat, "+" or "-".
    gene_strands : iterable of str
        Strands of genes whose bodies overlap the repeat locus. Empty
        means no overlap.

    Returns
    -------
    str
        "intergenic"  no gene body overlap
        "sense"       at least one same-strand overlap, no opposite-strand
        "antisense"   only opposite-strand overlaps
        "both"        overlaps on both strands (bidirectional loci)
    """
    strands = set(gene_strands)
    if not strands:
        return "intergenic"
    opposite = "-" if repeat_strand == "+" else "+"
    has_same = repeat_strand in strands
    has_opposite = opposite in strands
    if has_same and has_opposite:
        return "both"
    if has_same:
        return "sense"
    return "antisense"


_A_RUN = re.compile(r'A+', re.IGNORECASE)


def longest_a_run(sequence):
    """Return the length of the longest contiguous run of As in sequence.

    Case insensitive. Returns 0 for empty or None input, or when no A is
    present. Any non-A character (including N) breaks the run.

    For detecting internal polyA priming risk, pass the genomic sequence
    downstream of the repeat (in mRNA-sense orientation: reverse-complement
    the reference if the repeat is on the minus strand).
    """
    if not sequence:
        return 0
    best = 0
    for match in _A_RUN.finditer(sequence):
        run = match.end() - match.start()
        if run > best:
            best = run
    return best
