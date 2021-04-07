#!/bin/bash
if [ -z "$(ls -A automlsa2-examples)" ]; then
    git submodule update --init
fi

if [ -d "autotest" ]; then
    rm -rf autotest
fi

automlsa2 -t 4 --allow_missing 1 --dir automlsa2-examples/db --query automlsa2-examples/query/queries.fas -- autotest
