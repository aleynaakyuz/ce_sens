#!/usr/bin/env

python BNS_work.py --injections /home/aakyuz/BNS_workflow_inj/* --workflow-name test \
 --config-files /home/aakyuz/runs7/PE/workflow_BNS/static/CE4020/workflow.ini \
 --inference-config /home/aakyuz/runs7/PE/workflow_BNS/static/CE4020/relative.ini --output-dir /home/aakyuz/runs7/PE/workflow_BNS/static/CE4020/workflow_outputs --submit-now \
 --psd-paths /home/aakyuz/runs7/ce_strain/cosmic_explorer_strain.txt /home/aakyuz/runs7/ce_strain/cosmic_explorer_20km_strain.txt \
 --detectors CE40 CE20
