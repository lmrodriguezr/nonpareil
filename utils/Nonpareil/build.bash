#!/bin/bash

echo "
library(inlinedocs)
package.skeleton.dx('$(dirname -- $0)/');
" | R --vanilla

