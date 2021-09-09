# Some pointers for devs

## Environment polluters
R makes .Rproj and more session files that are automatically loaded upon startup. This is bad as they contain `renv` activation lines that may be sourced before dependencies are loaded. This means two things:
    - Always run R and Rscript using the `--vanilla` option, so that these files are not loaded or saved.
    - Always `.gitignore` and `.dockerignore` these files (`*/.Rprofile` and `*/*.Rproj`) so that they do not propagate.
Together we can stop the spread of these files and the contamination of more environments!

## Rebuilding the Docker
Docker saves intermediate containers in the cache, so that changes in steps further on in the dockerfile will cause the container to rebuild completely.

What does this mean for us? It means that the very long installation process doesn't need to be restarted if the installation files, namely `install_dep`, `install.R` and `renv.lock` don't change.

So, feel free to change any other file without having to worry about lengthy recompilation. This also means that tests can (and should) be run in the container and not in interactive mode (like Rstudio).

## Releases
To release a new build, follow these steps:
1. On github, stage a new release. Take note of the tag (ex: 1.0.0).
2. Download or pull the release commit (`git fetch` and `git merge <commit hash>`).
3. Trigger a docker rebuild with `docker build gattaca:<version>`. Replace `<version>` with the release version (ex: `docker build gattaca:1.0.0`). It could be useful to run the rebuild without the cache (using the `--no-cache` option).
4. Tag the container for release: `docker tag gattaca:<version> mrhedmad/gattaca:<version>` (ex: `docker tag gattaca:1.0.0 mrhedmad/gattaca:1.0.0`).
5. Push the container to Docker hub: `docker push`
6. ???
7. Profit!

If pushing returns an error, it might be because you have not logged in: run `docker login`.

We could also setup an automatic build, but that requires cash :3
