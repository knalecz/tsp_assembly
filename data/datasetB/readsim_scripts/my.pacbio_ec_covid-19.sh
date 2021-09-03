#!/bin/bash

#set -xv 


cd /home/knalecz/Pulpit/QuantumComputing/quantum_gsp/data/datasetC

for l in 5000; do 
  for c in 10; do 

    /opt/readsim-1.6/src/readsim.py sim fa \
    --ref covid19.fna \
    --pre covid19.pacbio_ec.reads.cov10.rl5000 \
    --rev_strd off \
    --tech pacbio_ec --read_mu $l --cov_mu $c --read_dist uniform

  done;
done;

