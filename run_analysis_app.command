#!/bin/bash
# Double-click this in Finder to launch the CosMx Analysis app.
# Opens a Terminal window, runs the app, and opens your browser once ready.
cd "$(dirname "$0")"

# Finder-launched Terminal windows don't always inherit the same PATH as an
# interactive shell (depends on how R was installed / which rc file sets
# PATH), so `Rscript` can silently fail to be found here even though it
# works fine from a normal terminal. Try the common install locations
# before giving up.
if ! command -v Rscript >/dev/null 2>&1; then
  for candidate in /usr/local/bin/Rscript /opt/homebrew/bin/Rscript \
                   /Library/Frameworks/R.framework/Resources/bin/Rscript; do
    if [ -x "$candidate" ]; then
      export PATH="$(dirname "$candidate"):$PATH"
      break
    fi
  done
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript not found on PATH, and not in any common install location."
  echo "Is R installed? (https://www.r-project.org)"
  read -n 1 -s -r -p "Press any key to close this window..."
  exit 1
fi

Rscript run_app.R analysis "$@"
status=$?

# Keep the window open regardless of outcome -- a Terminal window that
# closes itself the instant the script exits (success or failure) is why
# errors can go unseen.
if [ $status -ne 0 ]; then
  echo ""
  echo "App exited with an error (status $status). See the output above."
fi
read -n 1 -s -r -p "Press any key to close this window..."
