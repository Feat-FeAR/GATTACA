#!/usr/bin/env bash

# NOTE: Do **not** echo anything here other than the cli script.
# Everything that is echo-ed here will (or could) be captured and evaluated
# as code.

if [ "$1" == "getcli" ]; then
    # This is a request for the GATTACA script to be evaluated outside
    # the docker. We return it to be evaluated.
    echo "$(cat ./src/cli.sh)"
    exit 0
fi

# Note: Echo-ing is fine from now on.

# If we get to here, this is a normal call.
/usr/bin/env bash ./src/entry.sh "$@"
