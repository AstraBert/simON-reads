## Composition

The names of the loci are hereby reported, alongside with the species they where taken from (to simulate at best an unknown species from Coleoptera, we mixed 2 genes from 2 different taxa that belong to the order) and with their length in bp. All the reference sequences are stored in reference.fasta and were taken, without further modification, from NCBI-Nucleotide. 

| GENE | SPECIES | LENGTH(bp) |
|----------------|----------------|----------------|
| 28S-rRNA | *Scepticus uniformis* | 660 |
| S7 | *Glaresis ecostata* | 875 |

## Generate data

To generate data in Linux, run:

```bash
simON_reads.py -i reference.fasta -snp 28S-rRNA:3:G>T,28S-rRNA:41:C>G,S7:0:A>G,S7:63:C>T -n 1000 > test.fastq
```

If you are in Windows:
```powershell
python3 .\simON_reads.py -i .\reference.fasta -snp 28S-rRNA:3:G>T,28S-rRNA:41:C>G,S7:0:A>G,S7:63:C>T -n 1000 > .\test.fastq
```

Always remember to redirect the stream to your desidered file, unless you want it to be printed on the standard output of your terminal.
