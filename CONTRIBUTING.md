# Some pointers for devs

## Environment polluters in the docker
R makes .Rproj and more session files that are automatically loaded upon startup. This is bad as they contain `renv` activation lines that may be sourced before dependencies are loaded in the docker. This means two things:
    - Always run Rscript using the `--vanilla` option, so that these files are not loaded or saved.
    - Always `.gitignore` and `.dockerignore` these files (`*/.Rprofile` and `*/*.Rproj`) so that they do not propagate.

## Working and testing locally
When working on the environment locally, restore a local environment by running `renv::restore()` from inside `./GATTACA`. `renv` will create a .Rprofile file in the root folder which will be ran every time Rstudio loads a project in the folder. Simply load a project inside the root folder and Rstudio will handle it. 

## Rebuilding the Docker
Docker saves intermediate containers in the cache, so that changes in steps further on in the dockerfile will cause the container to rebuild completely.

What does this mean for us? It means that the very long installation process doesn't need to be restarted if the installation files, `dock_install.R` and `renv.lock` don't change.

So, feel free to change any other file without having to worry about lengthy recompilation. This also means that tests can (and should) be run in the container and not in interactive mode (like Rstudio).

Docker builds made during development should always be tagged with `bleeding`.

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
