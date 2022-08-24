#!/bin/bash

prefix=$1
python Scripts/make_config.py --meta meta.tab --exp_type generic --template Templates/template_TETranscripts.json --prefix $prefix
