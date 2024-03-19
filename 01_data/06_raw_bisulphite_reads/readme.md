Symlinks to raw sequence data files.

These are sorted by plate ID on the Nordborg NGS master list:
https://docs.google.com/spreadsheets/d/1XjO8zabj-1vlu-ex37MRnsnXB_c1U3_k-uXoeKaeGn0/edit#gid=26733257

sequence_run_info.csv gives information on the codes from the NGS sequencing facility:
- plate: plate ID from the NGS master list above
- request: NGS request number
- flowcell
- lane (nested within flow cell)
- date
- library (nested within request)
- str: concatenation of other columns for easier parsing
- index_set: Index set from the Unique Nextera Dual XT indices.
