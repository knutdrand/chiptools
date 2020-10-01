from collections import defaultdict, Counter, namedtuple
from itertools import dropwhile
from more_itertools import pairwise

GTFEntry = namedtuple("GTFEntry", ["chrom", "source", "feature", "left", "right", "score", "strand", "frame", "values"])

def parse_attributes(text):
    pairs = (tuple(pair.strip().split()) for pair in text.split(";") if pair.strip())
    return {key: value.strip('"') for key, value in pairs}

def parse_gtf_file(gtf_file):
    lines = (line.split("#")[0] for line in gtf_file)
    lines = (line.strip() for line in lines if line)
    parted = (line.split("\t") for line in lines)
    entries = (GTFEntry(chrom, source, feature, int(left)-1, int(right)-1, score, strand, frame, parse_attributes(values))
               for chrom, source, feature, left, right, score, strand, frame, values in parted)
    return entries

def get_junctions(exons):
    s = sorted(exons, key=lambda x: x.left)
    return ((first.chrom, first.right, second.left) for first, second in pairwise(s)
            if second.left>first.right+5)

def get_splice_sites(entries):
    genes = defaultdict(list)
    transcripts = defaultdict(list)
    exons = (entry for entry in entries if entry.feature=="exon")
    exons = (exon for exon in exons if "gene_id" in exon.values and "transcript_id" in exon.values)
    for exon in exons:
        transcripts[exon.values["transcript_id"]].append(exon)
        genes[exon.values["transcript_id"]].append(exon)
    splice_dict = defaultdict(list)
    for transcript, exons in transcripts.items():
        for junction in get_junctions(exons):
            splice_dict[junction].append(transcript)
    return splice_dict

def get_transcript_dict(entries):
    transcript_dict = defaultdict(list)
    transcripts = (e for e in entries if e.feature=="transcript")
    for transcript in transcripts:
        transcript_dict[transcript.values["gene_id"]].append(transcript.values["transcript_id"])
    return transcript_dict
