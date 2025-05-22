import DNA_Sequence
import DNA_interface
import random

class Genome(DNA_Sequence.Sequence):
    @classmethod #I'm not sure why this is necessary if you could help explain it
    def new_genome(cls, name, gene_count = 5, min_length=10, max_length=50, gc_content=0.5): #creates a new genome with n sequences that vary in length between two values

        if gene_count <= 0:
            raise ValueError("Sequence count must be greater than 0")
        if min_length <= 0 or max_length <= 0:
            raise ValueError("Sequence length must be greater than 0")
        if min_length > max_length:
            raise ValueError("Minimum must be smaller than maximum length")

        genome = cls(name)

        for i in range(gene_count):
            seq_length = random.randint(min_length, max_length) #I want to weight this length generator to better represent the distribution of genes but haven't done the research yet
            new_seq = DNA_Sequence.Sequence.random(seq_length, gc_content) #populates the sequences with bases
            genome.add_sequence(new_seq)
        return genome
#so far this just makes a sequence, I haven't figured out how I'm going to determine functional units or TEs and assign them fitness effects, I think it might just be through attributes but if you have any ideas let me know
    def __init__(self, name = '', sequences = []):
        self.name = name
        self.sequences = sequences or []
        self.sequence = []
        self._update_composite_sequence()
        self.genome_sequence = self.sequence

    def _update_composite_sequence(self):
        genome_sequence = []

        for seq in self.sequences:
            if isinstance(seq, DNA_Sequence.Sequence):
                genome_sequence.extend(seq.sequence)
            elif isinstance(seq, str):
                genome_sequence.extend([n.upper() for n in seq])
            elif isinstance(seq, list):
                genome_sequence.extend([n.upper() if isinstance(n, str) else n for n in seq])
            else:
                raise TypeError("Sequence type not supported")

        self.sequence = genome_sequence
        self.genome_sequence = genome_sequence

    def add_sequence(self, add_seq):
        if isinstance(add_seq, str) or isinstance(add_seq, list):
            add_seq = DNA_Sequence.Sequence(add_seq)
        elif not isinstance(add_seq, DNA_Sequence.Sequence):
            raise TypeError("Sequence type not supported")
        self.sequences.append(add_seq)
        self._update_composite_sequence()

    def get_sequence_at_index(self, index):
        if 0 <= index < len(self.sequences):
            return self.sequences[index]
        else:
            raise IndexError("Sequence index out of range")

    def get_sequence_by_position(self, position): #This returns which sequence is contains a certain position
        if position < 0 or position >= len(self.sequence):
            raise ValueError("Sequence position out of range")
        current_position = 0
        for seq in self.sequences: #looks at the sequences in order and sums their lengths until the position of interest lies within the next sequence
            seq_len = len(seq)
            if current_position <= position < current_position + seq_len:
                local_position = position - current_position #gives the position of our site of interest within its sequence
                return seq, local_position
            current_position += seq_len

    def print_sequence(self, line_length = 80): #prints the composite genome
        if not self.genome_sequence:
            print("No sequence")
            return
        if isinstance(self.genome_sequence[0], str):
            sequence_str = ''.join(self.genome_sequence)
        else:
            sequence_str = str(self.genome_sequence)

        for i in range(0, len(sequence_str), line_length):
            print(sequence_str[i:i+line_length])

    def __str__(self): #prints the genome name and its stats
        result = f"Genome: {self.name}\n"
        result += f"Total Length: {len(self.genome_sequence)} bp\n"
        result += f"Number of Sequences: {len(self.sequences)}\n"
        result += (f"GC Content: {self.gc_content():.2f}%")
        return result

    def insert_at_motif(genome, motif, insert): #Allows the class to locate a motif within the genome and insert a new sequence at that location
        genome_str = ''.join(genome.genome_sequence) #makes both the genome and the motif into strings so I can use the find motif method from the sequence class
        if isinstance (motif, DNA_Sequence.Sequence):
            motif_str = str(motif)
        else:
            motif_str = motif

        motif = motif.upper()
        position = genome_str.find(motif) #returns the positions where the motif is located

        if position == -1:
            print("Motif not found")
            return False, -1

        seq_obj, local_position = genome.get_sequence_by_position(position) #determines the sequence where the motif is located
        if seq_obj is None:
            print("Motif not found")
            return False, -1

        if isinstance(insert, str): #turns our insert into a sequence (used for testing but eventually the inserts will be TE sequences)
            insert = DNA_Sequence.Sequence(insert)

        original_sequences = genome.sequences.copy() #stores the original sequence list

        genome.sequences = []
        current_position = 0
        for i, seq in enumerate(original_sequences): #iterates through the original sequences reconstructing the genome until we reach the sequence with the split motif
            if current_position + len(seq) > position:
                break
            genome.add_sequence(seq)
            current_position += len(seq)

        first_part, second_part = seq_obj.split_at(local_position + len(motif)) #splits the sequence at the end of the motif

        genome.add_sequence(first_part)
        genome.add_sequence(insert)
        genome.add_sequence(second_part) #adds the split sequence and inserted sequence to the composite sequence

        for seq in original_sequences[i+1:]:
            genome.add_sequence(seq)
        genome._update_composite_sequence()

        return True, position
#TO DO
#I haven't tested how this handles multiple locations yet, I might write some code that assigns each motif a percentage chance of being a target and then uses a random number to determine which ones get activated and then loop the insert code
genome = Genome('Example')
seq2 = DNA_Sequence.Sequence("GGGCCGCAA")
genome.add_sequence(seq2)
insert = DNA_Sequence.Sequence("TTTTTTTTTTTTTT")
motif = "CCGC"

genome_str = ''.join(genome.genome_sequence)
print(f"Motif: {motif} in {genome_str}")

test_position = genome_str.find(motif)
print(f"Test Position: {test_position}")
print(f"Sequence type: {genome.genome_sequence}")
print(f"Sequence: {genome.genome_sequence[0:9]}")
for i, seq in enumerate(genome.sequences):
    print(f"Sequence: {i}: {seq}")

success, position = genome.insert_at_motif(motif, insert)

if success:
    print("\n" + str(genome))
    print("\nGenome Sequence:")
    genome.print_sequence()
else:
    print(f'\nfailed to find motif "{motif}" in genome')


