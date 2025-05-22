import random
import DNA_interface

#Should Each Sequence have some base attributes that contribute to the overall fitness of the individual?
#If there is variation among the sequences in terms of how much they contribute to fitness they will evolve once reproduction of individuals happens
#TE's could have a negative affect on fitness that models how they interfere with regular genetic activity
#TE's spread in Genomes until they get hunted out?
#Do the genes that regulate TE activity hurt overall fitness?
#something like this
#__slots__ = ['variant', 'resistance to TE targeting', 'individual contribution to genome fitness',
#I have realized I'm attempting to procedurally generate DNA. This may be a mistake.

class Sequence(DNA_interface.DNAInterface):
    #Class of DNA Sequence#

    nucleotides = {"A","T","G","C"}
    complementary_nucleotides = {"G":"C","C":"G","A":"T","T":"A"}

    def __init__(self, sequence=''):

        if isinstance(sequence,str):
            self.sequence = [nucleotide.upper() for nucleotide in sequence]
        elif isinstance(sequence,list):
            self.sequence = [
                nucleotide.upper() if isinstance(nucleotide, str) else nucleotide for nucleotide in sequence
            ]
        else:
            raise TypeError("Sequence must be of type str or list")
        self._validate_sequence()

    @classmethod
    def generate_DNA(cls, length, gc_content=0.5):
        if not isinstance(length, int) or length < 0:
            raise ValueError("Length must be positive")
        nucleotides = ['A','T','G','C']
        g_prob = c_prob = gc_content/2
        a_prob = t_prob = (1- gc_content)/2
        weights = [a_prob, t_prob, g_prob, c_prob]

        sequence = [random.choices(nucleotides, weights=weights)[0] for _ in range(length)]
        return cls(sequence)

    def _validate_sequence(self):
        for nucleotide in self.sequence:
            if nucleotide not in self.nucleotides:
                raise ValueError(f"Invalid nucleotide: {nucleotide}. Valid options are {self.nucleotides}")

    def __str__(self): #Returns the sequence as a string
        return ''.join(self.sequence)

    def __len__(self): #return the length of a sequence
        return len(self.sequence)

    def add(self, other): #connect two sequences
        if isinstance(other, Sequence):
            return Sequence(self.sequence + other.sequence)
        elif isinstance(other, str):
            return Sequence(self.sequence + other)
        else:
            raise TypeError(f"Cannot add sequence with type {type(other)}")

    def split_at(self, position): #I have given up on the _, I don't entirely understand their purpose
        if not isinstance(position, int):
            raise TypeError("Position must be an integer")
        if position < 0 or position > len(self.sequence):
            raise IndexError(f"Position {position} must be within a sequence")

        first_part = self.sequence[:position]
        second_part = self.sequence[position:]

        first_sequence = Sequence(''.join(first_part))
        second_sequence = Sequence(''.join(second_part))
        return (first_sequence, second_sequence)

    def get_complement(self):
        complementary = ''
        for nucleotide in self.sequence:
            complementary += self.complementary_nucleotides[nucleotide]
        return Sequence(complementary)

    def count_nucleotides(self):
        counts = {"A":0, "C":0, "G":0, "T":0}
        for nucleotide in self.sequence:
            counts[nucleotide] += 1
        return counts

    def gc_content(self):
        if not self.sequence:
            return 0
        counts = self.count_nucleotides()
        return (counts["G"]+counts["C"])/len(self.sequence)*100

    def find_forward_motif(self, motif):
        positions = []
        motif = motif.upper()
        sequence_str = ''.join(self.sequence)

        for i in range(len(self.sequence)-len(motif)+1):
            if sequence_str[i:i+len(motif)] == motif:
                positions.append(i)
        return positions

    def find_reverse_motif(self, motif):
        reverse_motif = ''
        motif_length = len(motif)
        for nucleotide in reversed(motif):
            if nucleotide in self.complementary_nucleotides:
                reverse_motif += self.complementary_nucleotides[nucleotide]
        sequence_str = str(self)
        reverse_positions = []

        for i in range(len(self.sequence)-len(motif)+1):
            if sequence_str[i:i+len(motif)] == reverse_motif:
                reverse_positions.append(i+len(motif))
        return reverse_positions


    def insert (self, position, nucleotides): #updates existing sequence

        if isinstance(nucleotides, str):
            nuc_insert = [n.upper() for n in nucleotides]
        elif isinstance(nucleotides, Sequence):
            nuc_insert = nucleotides.sequence
        else:
            raise TypeError("Insert must be a string or a Sequence")

        self.sequence = self.sequence[:position] + nuc_insert + self.sequence[position:]
        return self


#Things to maybe add, a way of identifying which nucleotides are at certain places in the sequence
# - Merge recombination eventually?
#some insert function, although I want to keep my TEs as seperate objects so maybe not
#an insert, delete, and replace function might be useful for those "hunter" classes that could turn TEs off
