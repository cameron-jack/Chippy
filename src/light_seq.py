from numpy import array
from cogent.core import moltype

class LightSeq(object):

    def __init__(self, Seq = '', Name = None, Quality = ''):

        if isinstance(Seq, LightSeq):
            Seq = Seq.seq
            Name = Seq.Name
            Quality = Seq.Quality

        else:
            if Name is None and hasattr(Seq, 'Name'):
                Name = Seq.Name

            if not Seq.isupper():
                Seq = Seq.upper()

            if Quality is not '':
                assert len(Quality) == len(Seq)

        self.seq = Seq
        self.Name = Name
        self.Quality = Quality


    def __getitem__(self, index):
        """ability to index and slice this sequence object. Returns another
        LightSeq object"""

        seq = self.seq[index]
        name = self.Name
        if self.Quality is not '':
            quality = self.Quality[index]
        else:
            quality = ''

        return LightSeq(seq, name, quality)


    def __str__(self):
        """__str__ returns self.seq unmodified."""
        return 'Name:\t\t%s\nSequence:\t%s\nQuality:\t%s\n'%\
               (self.Name, self.seq, self.Quality)


    def get_numeric_quality(self, scheme='illumina'):
        """Quality scores are stored as ascii strings. This method returns a
        list of the quality scores as numbers depending on the scheme specified.
        Currently, only phred and illumina are supported. """

        if self.Quality is not '':
            if scheme is 'illumina':
                return [v-64 for v in map(ord, self.Quality)]
            if scheme is 'phred':
                return [v-33 for v in map(ord, self.Quality)]
        return None



#name = '@GAPC_0015:6:1:1259:10413#0/1'
#seq = 'AACACCAAACTTCTCCACCACGTGAGCTACAAAAG'
#quality = '````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF'

#sequence = LightSeq(Seq=seq, Name=name, Quality=quality)
#print sequence
#print sequence.get_numeric_quality()
#trimSeq = sequence[5:15]
#print trimSeq
#print trimSeq.get_numeric_quality()









