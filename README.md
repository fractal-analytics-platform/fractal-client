# Fractal

[![PyPI version](https://badge.fury.io/py/fractal-client.svg)](https://badge.fury.io/py/fractal-client)

Fractal is a framework to process high content screening data at scale and prepares it for interactive visualization.

Fractal provides distributed workflows that convert TBs of image data into OME-Zarr files. The platform then processes the 3D image data by applying tasks like illumination correction and maximum intensity projection. The pyramidal OME-Zarr files enable interactive visualization in the napari viewer.
We are building towards integrating object segmentation (nuclei, cells, organoids) and feature measurements into Fractal.

Fractal is currently an **early alpha build**. We currently support only Yokogawa CV7000 image data as an input. Also, we're in the process of refactoring the workflow management into a client-server architecture. Thus, proceed with it at your own risk, there will still be significant breaking changes. While we don't have any releases or stable versions and thus don't provide user support, we welcome questions and discussions. Open an issue to get in touch.


![Fractal_multiwell_plate](https://user-images.githubusercontent.com/18033446/177169496-09784413-6ba9-4041-80e2-c70a70e0a5d9.gif)

Shortened movie of browsing an OME-Zarr file generated with Fractal in napari, using the [napari-ome-zarr plugin](https://github.com/ome/napari-ome-zarr). Actual loading times vary and can be a bit slower than in this GIF.

### Contributors
Fractal was conceived in the Liberali Lab at the Friedrich Miescher Institute for Biomedical Research and in the Pelkmans Lab at the University of Zurich (both in Switzerland). The project lead is with [@gusqgm](https://github.com/gusqgm) & [@jluethi](https://github.com/jluethi). The project was originally led by [@dvischi](https://github.com/dvischi).
The core development is done under contract by [@mfranzon](https://github.com/mfranzon), [@tcompa](https://github.com/tcompa) & [jacopo-exact](https://github.com/jacopo-exact) from eXact lab S.r.l. <exact-lab.it>.

*Installation instructions below have not been fully updated yet*

-----------------------------

### Requirements and configuration

- Install poetry, if you don't already have it. Here the official [guide](https://python-poetry.org/docs/)

- Move into the the folder
```
cd mwe_fractal/
```
- In your working (virtual) environment install all the dependencies with

```
poetry install
```
(takes few minutes)
- Install `graphviz`. This can be done via `conda install -c anaconda graphviz`, for instance, or via `sudo apt-get install graphviz` (or equivalent commands on other systems).

- Define some global configuration parameters in `fractal/fractal_config.py`. The essential ones are `partition` (the name of the SLURM partition on your cluster) and `worker_init` (which typically includes the activation of a virtual environment). An example is:
```
# Parameters of parsl.executors.HighThroughputExecutor
max_workers = 32
# Parameters of parsl.providers.SlurmProvider
# Note that worker_init is a command which is included at the beginning of
# each SLURM submission scripts
nodes_per_block = 1
cores_per_node = 16
mem_per_node_GB = 64
partition = "main"
worker_init = "source /opt/easybuild/software/Anaconda3/2019.07/"
worker_init += "etc/profile.d/conda.sh\n"
worker_init += "conda activate fractal"
```

## `fractal_cmd.py` playbook

On the previous terminal window, create the first project.

Move into the fractal folder.

First of all have a look on all the possible
commands:

```
python fractal_cmd.py
```

Typing ```--help``` after each command like:

```
python fractal_cmd.py project --help
```
you will see the arguments you had to pass.

### Let's create new project:

```
python fractal_cmd.py project new mwe-test path_output dstest
```
   - In which mwe-test is the name on the project
   - Then the path in which I want to put all the config files of my project (in this case I create a subfolder in the current folder, called ```path_output```)
   - Last one is the name of the first dataset of the project (assuming that each project has least one dataset)

- Now you should see two files:
   - One called fractal.json in the fractal/ directory
   - The second one, into the project folder you have choose called ```<project_name>.json```, in our case ```mwe-test.json```

```fractal.json``` stores all the general information like:
 - projects names, their path and the datasets name associated, then for each project the user that have created it
 - tasks, in this way tasks could be visible by everyone
 - users and groups

```mwe-test.json``` stores the specific information regarding current project:
 - datasets, their name, resources associated and the type
 - workflow (at this step workflow are project-specific)

With

```
python fractal_cmd.py project list
```
now you can see your projects.

### Datasets
Add resources to the dataset we have already created.

```
python fractal_cmd.py dataset add-resources mwe-test dstest absolute/path/of/resource
```

arguments are:
 - mwe-test, name of project
 - dstest, name of dataset
 - path to the folder in which are the files. For example the path to the folder in which there are tif images.

 Then you can update the type of the dataset, for example set it as ```tif``` or ```png```

```
python fractal_cmd.py dataset update-type mwe-test dstest tif
```

Now you can see the dataset object using list command passing as argument the project name:

```
python fractal_cmd.py dataset list mwe-test
```

### Tasks

***At the moment the task name and the unix executable which contains the logic of the tasks should have same name***

The tasks executable are in ```tasks/``` folder

Add a task; tasks are into the tasks folder. To add one or more just copy the filename without extension, example:

```
python fractal_cmd.py task add yokogawa_to_zarr zarr zarr well
```
The last three arguments are the input/output type and the parallelization level. Each task should act at the level of a `plate` or `well`, while `none` means that the task is not parallelized over any of its arguments.

For the moment the "depends_on" argument is not used.

```
python fractal_cmd.py task list
```
Now you should see the new task.

### Workflow
Let's create a workflow :

```
python fractal_cmd.py workflow new mwe-test wftest create_zarr_structure
```

 - project name
 - workflow name
 - tasks to add, in this case just one


*** Note: at the moment, all workflows must start with the `create_zarr_structure` task. ***


### Workflow execution

First, we need to specify a set of workflow-dependent parameters, which by now should be stored in a dedicated JSON file. An example of this file (for a workflow that performs yokogawa-to-zarr conversion and maximum intensity projection) reads
```
{
"workflow_name": "uzh_1_well_2x2_sites",
"dims": [2, 2],
"coarsening_factor_xy": 3,
"coarsening_factor_z": 1,
"num_levels": 5
}
```

After saving these parameters in `wf_params.json`, we execute the workflow through

```
python fractal_cmd.py workflow apply mwe-test wftest dstest dstest resource_in path_output wf_params.json
```

***PAY ATTENTION: `resource_in` should be a folder that is already inserted as resource in the input dataset, otherwise fractal_cmd will raise an error. Instead, the path_out folder will be created, if it does not exist.***

This command will use [parsl](https://parsl.readthedocs.io/en/stable/index.html) to run the tasks on a SLURM cluster. When the `fractal_cmd.py` script is complete, you will see a `runinfo` folder, in the current directory, with subfolder for each run (`000`, `001`, `002`..). This subfolder includes all logs from parsl, whose structure will be described in detail at a later stage.


### Setting up a local SLURM cluster

*** WARNING: This part has not been tested with the new parsl version of fractal! ***

To test slurm integration it is possible to run Fractal on a slurm node or to test it locally, creating a small slurm cluster using Docker/Podman

#### Docker setup

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

## Commands cli

```
fractal project new [project_name] [path] [dataset_name] - creates new project in path and create the first dataset object

fractal project list - lists registered projects

fractal dataset list [project_name] - lists datasets in project

fractal dataset new [project_name] [dataset_name] [resources] [dataset_type] - create a new dataset, add resources and the dataset type

fractal dataset add-resources [project_name] [dataset_name] [resources] - add resources to an existing dataset

fractal dataset update-type [project_name] [dataset_name] [dataset_type] - update type of dataset

fractal task list - list all the existing task

fractal task add [task_name] [input_type] [output_type] [parallelization_level] - create a new task

fractal workflow new [project_name] [workflow_name] [tasks] - create a new workflow to the project

fractal workflow list [project_name] - list all the workflow to that project

fractal workflow add-task [project_name] [workflow_name] [tasks_names] - add new tasks to the workflow

fractal workflow apply [project_name] [workflow_name] [dataset_name] [resource_in] [resource_out] [wf_parameters_json_file] - run workflow using the parameters into the .json

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
