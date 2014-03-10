#!/bin/bash

# Given a regex that object names have to match,
# gives you a list of object names that match it from ws-listobj

if [ $# -lt 2 ]; then
    echo "Usage: $0 [Workspace] [Regex]" 
    exit 1;
fi

WORKSPACE="$1";
REGEX="$2";

ws-listobj -w "$WORKSPACE" | grep -P '^\s*\d+\s+[^\s]*'"$REGEX"'[^\s]*\s' | grep -o -P '\s+[^\s]*'"$REGEX"'[^\s]*\s'