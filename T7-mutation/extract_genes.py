from Bio import SeqIO

whole_sequence = list(SeqIO.parse("T7_genome.gb", "genbank"))[0]
slices_of_sequence = []

for f in whole_sequence.features:
  if f.type=='CDS':
    print(f)
    #print(f.location.start, f.location.end)
    slice_of_sequence = whole_sequence[f.location.start:f.location.end]
    #print(f.qualifiers)
    slice_of_sequence.id = f.qualifiers['name'][0]
    slice_of_sequence.name = f.qualifiers['name'][0]
    slice_of_sequence.description = "Piece of T7 genome overlapping " + f.qualifiers['name'][0]
    slices_of_sequence.append(slice_of_sequence)

#print(slices_of_sequence)

SeqIO.write(slices_of_sequence, "T7_genes.gb", "genbank")
