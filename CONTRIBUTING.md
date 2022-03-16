# Some pointers for devs

Thank you for wanting to contribute to this repository. This file contains some guidelines for new submissions to the project, as well as well as a guide on how to work on the project locally.

## Environment polluters in the docker
R makes .Rproj and more session files that are automatically loaded upon startup. This is bad as they contain variables that migth have an effect on the analysis. This means two things:
    - Always run Rscript in the docker using the `--vanilla` option, so that these files are not loaded or saved.
    - Always `.gitignore` and `.dockerignore` these files (`*/.Rprofile` and `*/*.Rproj`) so that they do not propagate to the remote repository.
    - Remember to add the `-rm`

## Working and testing locally
It is recommended to rebuild the docker each time you need to test the package, even when working locally. See the next section.
A small testing suite can be run by executing `GATTACA test`. However, one should prepare integration tests on real data locally. For more information and help on this points, please contact us (e.g. send an email at luca.visentin (at) unito.it).

## Rebuilding the Docker
Docker saves intermediate containers in the cache, so that changes in the dockerfile will not (often) cause the container to rebuild completely.
Therefore, changing the files is `/src/` will not trigger a full, lengthy rebuild.
This allows development to be tested in the docker environment, and not locally.

Docker builds made during development should always be tagged with `bleeding`.

The usual workflow to work locally while rebuilding the docker instance is:
    - Make changes to the source files in `/src`;
    - Rebuild the docker with `docker build -t cmalabscience/gattaca:bleeding`;
    - Run `GATTACA -v bleeding ...` when testing.

## Rebuilding the annotations
The GATTACA internal annotations can be rebuilt with the scripts in the `/scripts/`
folder, and more specifically the `rebuild_annotations.sh` script. The process is
really *really* rough. You need both Python and R installed. There is no specific
list of required packages to rebuild the data (a good place to start is with the
`logger`, `progress` and `BiocManager` for R and `tqdm`, `requests`, `bs4` and
`jellyfish` for Python).

Just run `rebuild_annotations.sh` while your working directory is `/scripts/` and
hope for the best.

If annotation collisions are ever detected, it will crash. There is no agreed
method as of right now to deal with collisions.

This is probably the worst part of GATTACA, so contributions are welcome in making it better.

## Opening pull requests
Once you make a change, open a pull request. Be sure to fill out the default checklist.
