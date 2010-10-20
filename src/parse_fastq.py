#!/usr/bin/env python

from cogent import DNA
from cogent.parse.fasta import Info

from light_seq import LightSeq

def MinimalFastqParser(data, strict=True):
    """yields name, seq, qual from fastq file

    Arguments:
        - strict: checks the quality and sequence labels are the same
    """
    if type(data) == str:
        data = open(data)

    # fastq format is very simple, defined by blocks of 4 lines
    line_num = -1
    record = []
    for line in data:
        line_num += 1
        if line_num == 4:
            if strict: # make sure the seq and qual labels match
                assert record[0][1:] == record[2][1:], \
                  'Invalid format: %s -- %s' % (record[0][1:], record[2][1:])
            yield record[0][1:], record[1], record[3]
            
            line_num = 0
            record = []
        
        record.append(line.strip())
    
    if record:
        if strict and record[0]: # make sure the seq and qual labels match
            assert record[0][1:] == record[2][1:], 'Invalid format'
        
        if record[0]: # could be just an empty line at eof
            yield record[0][1:], record[1], record[3]
        
    
    if type(data) == file:
        data.close()


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
        make_seq = LightSeq

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

        yield make_seq(seq, name, qual)

if __name__ == "__main__":
    parser = MinimalFastqParser('../data/s_6_sequence.txt')
    r = 0
    for record in parser:
        r += 1
    print r