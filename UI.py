import DNA_interface
from DNA_Sequence import Sequence
from DNA_interface import DNAInterface


def give_me_stats_dna(dna_seq: DNA_interface.DNAInterface, motif):
    print(f"Sequence: {dna_seq}")
    print(f"complementary strand: {dna_seq.get_complement()}")
    print(f"GC Content equal to: {dna_seq.gc_content()}%")
    print(f"Locations of Motif:{motif} {dna_seq.find_forward_motif(motif)}")
    print(f"Locations of Reverse Motif:{motif} {dna_seq.find_reverse_motif(motif)} ")



dna_seq = Sequence.generate_DNA(1000, gc_content = 0.75)
#insert_1 = Sequence("CCCCC")
#insert_2 = Sequence("AAAAA")
#insert_3 = Sequence("GGGGG")

#dna_seq.insert(2, insert_1).insert(5, insert_2).insert(7, insert_3)

give_me_stats_dna(dna_seq, "ATG")






