#!/bin/bash

# -----------------------------------------------------------------------------
# Copy output files to user's home directory. If an output directory exists
# in the user's home directory, rename the existing output directory, then
# copy the new output directory.
# -----------------------------------------------------------------------------

if [ -d "$HOME_STEP_DIR/output" ]; then
    TIMESTAMP=$(date "+%Y%m%d-%H%M%S")
    mv ${HOME_STEP_DIR}/output ${HOME_STEP_DIR}/output_${TIMESTAMP}
fi

cp -r ${OUTPUT_DIR} ${HOME_STEP_DIR}/output
