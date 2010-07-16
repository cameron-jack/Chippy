from numpy import array
from cogent.core import moltype

class LightSeq(object):

    def __init__(self, Seq = '', Name = None, Quality = ''):

        if isinstance(Seq, LightSeq):
            Seq = Seq.Seq
            Name = Seq.Name
            Quality = Seq.Quality

        else:
            if Name is None and hasattr(Seq, 'Name'):
                Name = Seq.Name

            if not Seq.isupper():
                Seq = Seq.upper()

            if Quality is not '':
                assert len(Quality) == len(Seq)

        self.Seq = Seq
        self.Name = Name
        self.Quality = Quality

    def __len__(self):
        return len(self.Seq)
    
    def __getitem__(self, index):
        """ability to index and slice this sequence object. Returns another
        LightSeq object"""

        seq = self.Seq[index]
        name = self.Name
        if self.Quality is not '':
            quality = self.Quality[index]
        else:
            quality = ''

        return LightSeq(seq, name, quality)


    def __str__(self):
        """__str__ returns self.Seq unmodified."""
        return 'Name:\t\t%s\nSequence:\t%s\nQuality:\t%s\n'%\
               (self.Name, self.Seq, self.Quality)


    def getNumericQuality(self, scheme='illumina'):
        """Quality scores are stored as ascii strings. This method returns a
        list of the quality scores as numbers depending on the scheme specified.
        Currently, only phred and illumina are supported. """

        if self.Quality is not '':
            if scheme is 'illumina':
                return [v-64 for v in map(ord, self.Quality)]
            if scheme is 'phred':
                return [v-33 for v in map(ord, self.Quality)]
        return None

    def toFasta(self):
        """return a fasta formatted string from the sequence object"""
        return '\n'.join(['>%s' % self.Name, self.Seq])

    def toFastq(self, quality=None):
        """return a fastq formatted string from the sequence object. Quality
        can be specified. This method will also update the Quality of the
        LightSeq object if it was empty"""

        if quality is None:
            if self.Quality is '':
                print 'Quality Scores need to specified for fastq format'
                raise ValueError
        else:
            assert len(quality) == len(self.Seq)
            self.Quality = quality

        return '\n'.join(['>%s' % self.Name, self.Seq, '>%s' % self.Name, self.Quality])







