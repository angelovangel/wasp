# wasp

Align long or short reads to a reference entirely in the browser (minimap2 + samtools compiled to WebAssembly via [biowasm](https://biowasm.com)). No data leaves your computer. Suitable for up to bacterial genome-sized references. 

https://angelovangel.github.io/wasp/

## Performance and limits

Rough estimates, untested (10-core Mac, 16 GB RAM, Chrome):

| Reference | Depth | Time |
|---|---|---|
| Plasmid (5–200 kb) | 100–1000x | seconds |
| Bacterial genome (~5 Mb) | 50–100x | 1–5 min |

Memory is the limit. Expect crashes above ~2 Gb of kept reads (e.g. 5 Mb at 400x) or references above ~100 Mb. Keep bacterial genomes at 50–100x.
