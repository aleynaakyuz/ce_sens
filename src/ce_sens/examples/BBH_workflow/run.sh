#!/usr/bin/env

python BBH_work.py --detectors CE30A CE30B \
 --workflow-name test --config-files /home/aakyuz/runs7/PE/workflow_BBH/Sep14/runs/light/CE3030/5/workflow.ini \
 --inference-config /home/aakyuz/runs7/PE/workflow_BBH/Sep14/runs/light/CE3030/5/margtime.ini \
 --output-dir /home/aakyuz/runs7/PE/workflow_BBH/Sep14/runs/light/CE3030/5/workflow_outputs --submit-now \
 --psd-paths /home/aakyuz/runs7/ce_strain/cosmic_explorer_30km_strain.txt /home/aakyuz/runs7/ce_strain/cosmic_explorer_30km_strain.txt \
 --injections /home/aakyuz/runs7/PE/workflow_BBH/Sep14/injections/light/CE3030/5/*
