# ORFfinder
Open reading frame finder
Developed a computational tool using python to identify Open Reading Frames (ORFs) in DNA sequences by implementing an algorithm that locates start and stop codons. The tool processes FASTA-formatted files, analyzes both strands of DNA across six possible coding frames, and outputs formatted data on ORFs that are 100 nucleotides or longer. This project involved enhancing a sequence analysis module, handling overlapping genes, and ensuring correct frame alignment between start and stop codons.
To run the finder: SierraDonnORFfinder.py  < tass2.fa > testoutputFile -mG 300 -lG -s ATG -s GTG
