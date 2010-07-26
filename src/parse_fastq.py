from cogent import DNA
from cogent.parse.fasta import Info

from light_seq import LightSeq

def MinimalFastqParser(data):
    """yields name, seq, qual from fastq file
    """
    if type(data) == str:
        data = open(data)

    seq = None
    qual = None
    seq_label = None
    qual_label = None
    num = 0
    for line in data:
        if line[0] == '@': # is seq block
            seq_label = line[1:].strip()
            qual_label = None
        elif line[0] == '+':# is qual block
            qual_label = line[1:].strip()
            assert qual_label == seq_label
        elif seq_label is not None and qual_label is None:
            seq = line.strip()
        elif seq_label is not None and qual_label is not None:
            qual = line.strip()
            yield seq_label, seq, qual
            seq = None
            qual = None
            seq_label = None
            qual_label = None



def _make_seq(seq, name, qual):
    """makes a DnaSeq object"""
    info = Info(qual=qual)
    seq = DNA.makeSequence(seq, name)
    seq.Info = info
    return name, seq

def FastqParser(data, remove_ambig=True, trim_bad_bases=True,
                make_seq=LightSeq):
    """yields name, seq from fastq file

    Arguments:
        - data is a data series (e.g. file, list)
        - remove_ambig: excludes any sequence with N
        - trim_bad_bases: trims terminal positions with B quality score
        - make_seq: a function that takes seq, name, qual and returns data
    """

    ambig_count = 0
    parser = MinimalFastqParser(data)
    for name, seq, qual in parser:
        if remove_ambig:
            if 'N' in seq:
                seq = None
                qual = None
                seq_label = None
                qual_label = None
                ambig_count+=1
                continue

        seq_object = make_seq(Seq=seq, Name=name, Quality=qual)

        if trim_bad_bases:
            # we find B and trim seq and qual
            bad = qual.find('B')
            if bad == 0:
                continue
            elif bad > 0:
                seq_object = seq_object[:bad]

        if len(seq_object.Seq) == 0:
            continue

        yield seq_object

