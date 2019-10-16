#!/bin/bash

set -e
errors=0

# Run unit tests
python snvly/snvly_test.py || {
    echo "'python python/snvly/snvly_test.py' failed"
    let errors+=1
}

# Check program style
pylint -E snvly/*.py || {
    echo 'pylint -E snvly/*.py failed'
    let errors+=1
}

[ "$errors" -gt 0 ] && {
    echo "There were $errors errors found"
    exit 1
}

echo "Ok : Python specific tests"
