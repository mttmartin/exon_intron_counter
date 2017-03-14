# exon_intron_counter

Counts the number of reads which map to the exons and introns of genes and outputs this information to a CSV file.

## Installation
You can clone this repository with:
```
$ git clone https://github.com/mttmartin/exon_intron_counter.git
```

## Usage

On Linux/OS X:
```
$ ./exon_intron_counter.py --gtf GTF_FILE.gtf -i INPUT_BAM.bam -o OUTPUT_FILE.csv
```

On all systems the following should work:
```
$ python3 exon_intron_counter.py --gtf GTF_FILE.gtf -i INPUT_BAM.bam -o OUTPUT_FILE.csv
```

### Arguments
<pre>
required argumnets:
  -g GTF, --gtf GTF     GTF file

optional arguments:
  -i INPUT, --input INPUT
                        Input BAM file(default standard input)
  -o OUTPUT, --output OUTPUT
                        Output file(default standard output)
</pre>

## Output format
Output is in CSV format with the following columns: gene ID, exon count, intron count.

## License
Copyright (c) 2017 Matthew Martin

Licensed under GPLv3.
