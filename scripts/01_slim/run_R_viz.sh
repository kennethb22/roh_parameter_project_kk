#!/bin/bash
#
#   +-----------------------+
#   |  USE:                 |
#   |    - MEDIUM queue     |
#   |    - 1 CPU + 8 Gb     |
#   +-----------------------+
#
#  Replace the USER name in this script with your username and
#  call your project whatever you want
#
#  This script must be made executable like this
#    chmod +x my_script
#
#  Submit this script to the queue with a command like this
#    run_script my_script.sh
#
#  My preferred setup before running:
#    -- script to be run in /home/projectdir/scripts
#    -- project directory (of same name as script) in /home/
#    -- /input/ and /output/ subdirs within project dir


## --------------------------------
## Load R + run script 
module load R/4.1.0

R CMD BATCH /home/aubkbk001/roh_param_project/01_slim/scripts/slim_output_roh_viz_ASC.R

mail -s 'R viz finished' kirkseykb1@appstate.edu <<< 'R viz finished'
