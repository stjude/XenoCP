#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: Tests a bam for existence, sorting/indexing, and stats available via flagstat
baseCommand: qc_bam.sh

hints:
  SoftwareRequirement:
    packages:
      qclib:
        specs: ["qclib.sh"]
  DockerRequirement:
    dockerPull: "ghcr.io/stjude/xenocp:latest"

inputs:
  bam:
    type: File
    inputBinding:
      position: 1
  status:
    type: string
    default: indexed
    inputBinding:
      position: 2
  flagstat:
    type: File
    inputBinding:
      position: 3

outputs: []
