#!/bin/bash

DIR=$(dirname -- $0)/Nonpareil
echo "
library(inlinedocs)
package.skeleton.dx('${DIR}/');
" | R --vanilla

cat "${DIR}/man/nonpareil-package.Rd" | tr -d '\r' \
  | grep -v '^}$' | grep -v '^\\author{' \
  | grep -v '^Maintainer' \
  | perl -pe 's/^\\keyword/}\n\\author{Luis M. Rodriguez-R <lmrodriguezr\@gmail.com> [aut, cre]}\n\n\\keyword/' \
  | perl -lwe '$/=\0; $_=<>; s/^\\details{\n+([^}].*\n+)*}\n+//mg; print' \
  > o && mv o "${DIR}/man/nonpareil-package.Rd"
for i in "curve" "set" "curve.batch" ; do
  cat "${DIR}/man/Nonpareil.${i}.Rd" \
    | perl -pe 's/(.*)###<<< dontrun/\\dontrun{ $1 }/' \
    > o && mv o "${DIR}/man/Nonpareil.${i}.Rd"
done

