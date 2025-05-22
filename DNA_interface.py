from abc import ABC, abstractmethod

class DNAInterface(ABC):

    @abstractmethod
    def generate_DNA(self, length, gc_content=0.5):
        pass
    @abstractmethod
    def __str__(self):
        pass
    @abstractmethod
    def __len__(self):
        pass
    @abstractmethod
    def add(self, other):
        pass
    @abstractmethod
    def split_at(self, position):
        pass
    @abstractmethod
    def get_complement(self):
        pass
    @abstractmethod
    def count_nucleotides(self):
        pass
    @abstractmethod
    def gc_content(self):
        pass
    @abstractmethod
    def find_forward_motif(self):
        pass
    @abstractmethod
    def find_reverse_motif(self):
        pass
    @abstractmethod
    def insert(self, position, nucleotides):
        pass
