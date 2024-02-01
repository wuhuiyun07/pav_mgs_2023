#!/usr/bin/env bash
diamond blastx 
        --threads {threads}
        -q {input.fa} 
        -d {input.db} 
        -o {output.tsv}
        --header
        --outfmt {params.fmt}
        --verbose
        --compress 1
        --log {log}