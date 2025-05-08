#!/bin/bash

conda activate copernicusmarine
Rscript code/00_getCMEMS.R
conda deactivate