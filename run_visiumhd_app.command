#!/bin/bash
# Double-click to launch the Visium HD Analysis app (data from inst/scripts/visiumhd_osta_workflow.R).
cd "$(dirname "$0")"
if ! command -v Rscript >/dev/null 2>&1; then
  for c in /usr/local/bin/Rscript /opt/homebrew/bin/Rscript /Library/Frameworks/R.framework/Resources/bin/Rscript; do
    [ -x "$c" ] && export PATH="$(dirname "$c"):$PATH" && break
  done
fi
Rscript run_app.R analysis "${1:-../visiumhd_isee_ready.RDS}" --platform=visiumhd
read -n 1 -s -r -p "Press any key to close this window..."
