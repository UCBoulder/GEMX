#Store local code updates in run dir.
REPO_ROOT=$(git rev-parse --show-toplevel)
cp $REPO_ROOT/bin/codeChanges.txt .