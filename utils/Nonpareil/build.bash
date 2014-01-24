#!/bin/bash

echo "
library(inlinedocs)
package.skeleton.dx('$(dirname -- $0)/');
install.packages('$(dirname -- $0)/', repos=NULL);
" | R --vanilla

