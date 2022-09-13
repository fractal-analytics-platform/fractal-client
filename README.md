# Fractal
Fractal is a framework to process high content screening data at scale and prepares it for interactive visualization.

Fractal provides distributed workflows that convert TBs of image data into OME-Zarr files. The platform then processes the 3D image data by applying tasks like illumination correction and maximum intensity projection. The pyramidal OME-Zarr files enable interactive visualization in the napari viewer.
We are building towards integrating object segmentation (nuclei, cells, organoids) and feature measurements into Fractal.

Fractal is currently an **early alpha build**. We currently support only Yokogawa CV7000 image data as an input. Also, we're in the process of refactoring the workflow management into a client-server architecture. Thus, proceed with it at your own risk, there will still be significant breaking changes. While we don't have any releases or stable versions and thus don't provide user support, we welcome questions and discussions. Open an issue to get in touch.


![Fractal_multiwell_plate](https://user-images.githubusercontent.com/18033446/177169496-09784413-6ba9-4041-80e2-c70a70e0a5d9.gif)

Shortened movie of browsing an OME-Zarr file generated with Fractal in napari, using the [napari-ome-zarr plugin](https://github.com/ome/napari-ome-zarr). Actual loading times vary and can be a bit slower than in this GIF.

## Contributors
Fractal was conceived in the Liberali Lab at the Friedrich Miescher Institute for Biomedical Research and in the Pelkmans Lab at the University of Zurich (both in Switzerland). The project lead is with [@gusqgm](https://github.com/gusqgm) & [@jluethi](https://github.com/jluethi). The project was originally led by [@dvischi](https://github.com/dvischi).
The core development is done under contract by [@mfranzon](https://github.com/mfranzon), [@tcompa](https://github.com/tcompa) & [jacopo-exact](https://github.com/jacopo-exact) from eXact lab S.r.l. <exact-lab.it>.

-----------------------------

## How to install

Fractal is currently split into three repositories, corresponding to the client, server and tasks. There exist several ways of installing these three components, depending on whether you want to use local/github/pypi versions of the packages. Independently on this choice, it is convenient to run all the following commands from a fresh virtual environment (a standard python venv, or a conda one).

WARNING: As the three packages evolve, the following instructions may need to be updated!

### Development install (all from local folders)

This approach requires some small changes to the `pyproject.toml` file of two client/server repositories (`fractal` and `fractal-server`).
For simplicity, let us assume that the three repositories are in the same folder, and that we already issued `git checkout main` for the three of them (where `main` can be replaced by any relevant branch).

Instructions:

1. Install `poetry` version `1.2.0b2` (e.g. `pip install poetry==1.2.0b2`).

Note that we will switch to `1.2.0` (or patches) as soon as it becomes more stable (e.g. after this fix is released: https://github.com/python-poetry/poetry-core/pull/466), and hopefully we will then stick with it.

2. From the `fractal-server` folder, run
```
poetry add --editable ../fractal-tasks-core
```
and
```
poetry add --editable ../fractal
```
(see e.g. [here](https://github.com/python-poetry/poetry/discussions/1135) for more info about the `--editable` option)

If you now inspect the `pyproject.toml` file, you will notice lines like
```
fractal-client = {path="../fractal", develop=true}
...
fractal-tasks-core = {path = "../fractal-tasks-core", develop=true}
```
3. Still from the `fractal-server` folder, run
```
poetry install
```
4. From the `fractal` folder, run
```
poetry add --editable ../fractal-server/ --group dev
```
5. From the `fractal` folder, run
```
poetry install
```


### Development install (all from github branches)

Similar to the previous one, but `path` points to a github repo and you can also specify the branch
Full instructions will be added later.

### Fully pypi

Assuming that `pip` is available on your system, you can do:
```
pip install fractal-server fractal-client fractal-tasks-core
```
and you will get the most recent PyPI release of these packages.


### TO CHECK

What about `graphviz`? Is it necessary to install it on a side (with `sudo apt-get install graphviz` or equivalent commands on other systems)?


## OLD: Setting up a local SLURM cluster

*** WARNING: This part has not been tested with the new parsl version of fractal! ***

To test slurm integration it is possible to run Fractal on a slurm node or to test it locally, creating a small slurm cluster using Docker/Podman

### Docker setup

To create your testbed with a slurm cluster you have to install Docker :
[Guide](https://docs.docker.com/get-docker/)

Then in your working environment install docker-compose: [Guide](https://docs.docker.com/compose/install/)

Now just deploy the docker-slurm.yaml file which contains four containers:

- Master node with a jupyter hub; in this way will be easy launch scripts / create ones, and monitor slurm cluster.
- Three worker nodes.

To deploy it just use docker-compose :
```
docker-compose -f docker-slurm.yaml up -d
```

### Podman setup
As docker case, it requires Podman and podman-compose installed. Here the two references:

Podman : [Guide](https://podman.io/getting-started/installation)

podman-compose : [Guide](https://github.com/containers/podman-compose)

To deploy the .yaml file just use podman-compose :

```
podman-compose -f docker-slurm.yaml up -d
```

## Misc

To easily load all channels of the 2D example in Napari you can just run this piece of code in a cell of a Jupyter notebook or run it as an executable.

```
import napari
import zarr

x = zarr.open('path/to/201703012-NUPfollowup.zarr')
viewer = napari.Viewer()
[viewer.add_image(zarr.open(f'path/to/201703012-NUPfollowup.zarr/C22/{n}', 'r'), name=f'{n}', blending="additive") for n in x.C22]
```

## Development

Development takes place on Github. You are welcome to submit an issue and open
pull requests.

### Releasing

Before release checklist:

- [ ] The `main` branch is checked out
- [ ] You reviewed dependencies and dev dependencies and the lock file is up to
      date with `pyproject.toml`.
- [ ] The current `HEAD` of the main branch passes all the tests
- [ ] Use
```
poetry run bumpver update --dry --[patch|minor] --tag-commit --commit
```
to test updating the version bump
- [ ] If the previous step looks good, use
```
poetry run bumpver update --[patch|minor] --tag-commit --commit
```
to actually bump the version and commit the changes locally.
- [ ] Test the build with
```
poetry build
```
- [ ] If the previous step was successful, push the version bump and tags:
```
git push && git push --tags
```
- [ ] Finally, publish the updated package to pypi with
```
poetry publish --dry-run
```
removing the `--dry-run` when you made sure that everything looks good.
