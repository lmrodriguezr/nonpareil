#!/bin/bash

DIR=$(dirname -- $0)/Nonpareil
cd $DIR && echo "devtools::document()" | R --vanilla

