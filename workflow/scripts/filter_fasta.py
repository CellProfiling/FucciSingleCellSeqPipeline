import sys, io
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
fasta=sys.argv[1]
output=sys.argv[2]
fasta_sequences = SeqIO.parse(open(fasta),'fasta')

with open(output,"w") as out:
    for seq_record in fasta_sequences:
        if seq_record.id.startswith("20") or seq_record.id.startswith("21") or seq_record.id.startswith("22"):
            out.write(seq_record.format("fasta"))
