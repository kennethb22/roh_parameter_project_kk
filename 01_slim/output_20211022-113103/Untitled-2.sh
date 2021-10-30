#!/bin/bash
#
# initialize_script_variables.sh
# 
# This script initializes variables that are used in multiple scripts in the 
# ROH calling workflow. It should be included in all scripts for each step of
# the workflow. Put the following line in the scripts:
#
#     source /home/aubkbk001/roh_param_project/init_script_vars.sh
#
# replacing aubkbk001 with your username.

##  Set username
USER=aubkbk001

## Set project name
PROJECT=roh_param_project

## Set user home directory
USER_HOME_DIR=/home/${USER}/



