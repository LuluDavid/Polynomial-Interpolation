#!/bin/sh

# Commit all pending modifications
git add -A
git commit -m "Auto-commit all pending changes"
# Delete all local cached remote resources
git rm -r --cached .
# Re-push all local resources applying .gitignore
git add .
git commit -m ".gitignore is now working"
git push