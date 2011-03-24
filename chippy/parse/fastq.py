#!/usr/bin/env python

from cogent import DNA
from cogent.parse.fasta import Info
from cogent.parse.fastq import MinimalFastqParser

from chippy.parse.light_seq import LightSeq

def FastqParser(data, numeric_qual=False, remove_ambig=True,
    trim_bad_bases=True, make_seq=None, strict=True):
    """yields name, seq from fastq file

    Arguments:
        - data is a data series (e.g. file, list)
        - numeric_qual: whether base quality score is returned or raw ASCII
        - remove_ambig: excludes any sequence with N
        - trim_bad_bases: trims terminal positions with B quality score
        - make_seq: a function that takes name, seq, qual and returns data,
          defaults to the LightSeq class
    """

    if make_seq is None:
        light_seq = LightSeq()

    ambig_count = 0
    parser = MinimalFastqParser(data, strict=strict)
    for name, seq, qual in parser:
        if remove_ambig:
            if 'N' in seq:
                seq = None
                qual = None
                seq_label = None
                qual_label = None
                ambig_count+=1
                continue

        if trim_bad_bases:
            # we find B and trim seq and qual
            bad = qual.find('B')
            if bad == 0:
                continue
            elif bad > 0:
                seq = seq[: bad]
                qual = qual[: bad]

        if len(seq) == 0:
            continue

        if numeric_qual:
            qual = [v-64 for v in map(ord, qual)]

        light_seq(seq, name, qual)
        yield light_seq

if __name__ == "__main__":
    parser = MinimalFastqParser('../data/s_6_sequence.txt')
    r = 0
    for record in parser:
        r += 1
    print r