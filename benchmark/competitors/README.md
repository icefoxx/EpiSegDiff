# Competitor benchmark log

Every attempted comparator run must append a row to `runs.tsv` even if
the process fails. Do not drop unfavorable or crashed methods.

Required fields:

```text
method
software_version
command
parameters
input_transformation
runtime_sec
exit_status
output_parsing_rule
notes
```

Count-assay targets: csaw, MACS3+DiffBind, THOR.

Methylation targets: DSS, dmrseq, metilene (MethCP if installable).
