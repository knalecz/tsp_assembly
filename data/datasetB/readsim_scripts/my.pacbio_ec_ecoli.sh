#!/bin/bash

#set -xv 


cd /home/knalecz/Pulpit/QuantumComputing/quantum_gsp/data/datasetC

for l in 5000; do 
  for c in 5; do 

    /opt/readsim-1.6/src/readsim.py sim fa \
    --ref NC_000913.50k.fna \
    --pre NC_000913.50k.pacbio_ec.reads.L50000.cov5.rl5000 \
    --rev_strd off \
    --tech pacbio_ec --read_mu $l --cov_mu $c --read_dist uniform

  done;
done;

for l in 5000; do 
  for c in 10; do 

    /opt/readsim-1.6/src/readsim.py sim fa \
    --ref NC_000913.50k.fna \
    --pre NC_000913.50k.pacbio_ec.reads.L50000.cov10.rl5000 \
    --rev_strd off \
    --tech pacbio_ec --read_mu $l --cov_mu $c --read_dist uniform

  done;
done;

for l in 5000; do 
  for c in 15; do 

    /opt/readsim-1.6/src/readsim.py sim fa \
    --ref NC_000913.50k.fna \
    --pre NC_000913.50k.pacbio_ec.reads.L50000.cov15.rl5000 \
    --rev_strd off \
    --tech pacbio_ec --read_mu $l --cov_mu $c --read_dist uniform

  done;
done;

