# Fractal
## Minimal working example

### Requirements

- Install poetry, if you don't already have it. Here the official [guide](https://python-poetry.org/docs/)

- Move into the the folder
```
cd fractal/
```
- In your working (virtual) environment install all the dependencies with

```
poetry install
```
(takes few minutes)

- Open a new terminal window, activate the working environment then start the luigi deamon

```
luigid
```

## Playbook

On the previous terminal window, create the first project.

First of all have a look on all the possible
commands:

```
python fractal.py
```

Typing ```--help``` after each command like:

```
python fractal.py project --help
```
you will see the arguments you had to pass.

### Let's create new project:

```
python fractal.py project new mwe-test $PWD/../test-proj dstest
```
   - In which mwe-test is the name on the project
   - Than the path in which I want to put all the config files of my project (in this case I create a folder in the parent of the current folder, called ```test-proj```)
   - Last one is the name of the first dataset of the project (assuming that each project had to have at least one dataset)

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
python fractal.py project list
```
now you can see your projects.

### Datasets
Add resources to the dataset we have already created.

```
python fractal.py dataset add-resources mwe-test dstest absolute/path/of/resource
```

arguments are:
 - mwe-test, name of project
 - dstest, name of dataset
 - path to the folder in which are the files. For example the path to the folder in which there are tif images.

 Then you can update the type of the dataset, for example set it as ```tif```

```
python fractal.py dataset update-type mwe-test dstest tif
```

Now you can see the dataset object using list command passing as argument the project name:

```
python fractal.py dataset list mwe-test
```

### Tasks

***At the moment the task name and the unix executable which contains the logic of the tasks should have same name***

The tasks executable are in ```tasks/``` folder

Add a task; tasks are into the tasks folder. To add one or more just copy the filename without extension, example:

```
python fractal.py task add compression_tif tif tif
```

There is also the "Depends_on" argument which is optional, but for the moment we do not want to have dependencies on other tasks.

```
python fractal.py task list
```
Now you should see the new task.

### Workflow
Let's create a workflow :

```
python fractal.py workflow new mwe-test wftest compression_tif
```

 - project name
 - workflow name
 - tasks to add, in this case just one


Finally test the workflow.
In this folder there is a file called ```test_apply.json``` which is the template of the file needed to apply command.

### test_apply.json
This file containes the meta-information to correctly run the workflow
After that all the corrects names and paths are inserted, it is necessary add other few arguments:
-  ```delete``` argument refers to the possibility to delete or not the resource input. If it is set to ```True```, in the compression case, it deletes the images after the compression.

- ```scheduler```, which could be ```local``` or ```slurm```. If local, it will use just Luigi scheduler as backend.

- ```other_params```, which is a dictionary in which stores all the extra parameters that could be useful for a task.


```
python fractal.py apply test_apply.json
```

***PAY ATTENTION:
resource_in should be a folder that is already
inserted as resource in the input dataset, otherwise it raises an error.
Instead, the output folder, will be created, if it not exists.***

When the workflow is finished, you should see a ```log/``` folder in which there are a ```.txt``` file in which there is the standard error, if something goes wrong should be written here.

## Slurm integration
To use slurm as backend just specify ```slurm``` as scheduler in the ```test_apply.json``` file. This will create a slurm job for each Luigi task you have in the workflow. Luigi workers wait until the slurm job has finished then they return the status.

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

## Commands cli

```
fractal project new [project_name] [path] [dataset_name] - creates new project in path and create the first dataset object

fractal project list - lists registered projects


fractal dataset list [project_name] - lists datasets in project

fractal dataset new [project_name] [dataset_name] [resources] [dataset_type] - create a new dataset, add resources and the dataset type

fractal dataset add-resources [project_name] [dataset_name] [resources] - add resources to an existing dataset

fractal dataset update-type [project_name] [dataset_name] [dataset_type] - update type of dataset


fractal task list - list all the existing task

fractal task add [task_name] [input_type] [output_type] [depends_on] - create a new task


fractal workflow new [project_name] [workflow_name] [tasks] - create a new workflow to the project
fractal workflow list [project_name] - list all the workflow to that project
fractal workflow add-task [project_name] [workflow_name] [tasks_names] - add new tasks to the workflow


fractal apply [filename] - run workflow using the parameters into the .json
```
