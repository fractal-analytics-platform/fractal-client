import copy
import datetime
import os
import random
import re
import uuid
from glob import glob
from subprocess import PIPE  # nosec
from subprocess import Popen  # nosec

import jinja2
import luigi
import zarr
from devtools import debug


class ConversionTaskWrap(luigi.Task):

    """
    This class is a wrapper for the conversions task. It initialiazes the
    zarr folder then call the conversion task to populate the new zarr file
    with the arrays.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.complete_flag = False

    #: name of the task, it must be the same of one of the tasks into
    #: tasks folder
    task_name = luigi.parameter.Parameter(significant=True)

    #: name of the workflow
    wf_name = luigi.parameter.Parameter(significant=True)

    #: input path
    in_path = luigi.parameter.Parameter(significant=True)

    #: output path
    out_path = luigi.parameter.Parameter(significant=True)

    #: delete the input files
    delete_in = luigi.parameter.Parameter()

    #: scheduler, local or slurm
    sclr = luigi.parameter.Parameter()

    #: extension of the input files
    ext = luigi.parameter.Parameter()

    #: slurm configs parameters
    slurm_param = luigi.parameter.DictParameter()

    #: extra parameters rquired by tasks
    other_param = luigi.parameter.DictParameter()

    #: path in which are stored the tasks executable
    tasks_path = os.getcwd() + "/tasks/"

    #: complete or not the luigi task
    done = False

    def metadata(self, filename):
        """
        function to extract all the metadata stored
        in image's filename. Return a list with all params.

        :param filename: name of the image
        :type filename: str
        """
        f = filename.rsplit(".", 1)[0]
        plate = f.split("_")[0]
        well = re.findall(r"_(.*)_T", f)[0].split("_")[-1]
        site = re.findall(r"F(.*)L", f)[0]
        chl = re.findall(r"[0-9]C(.*)", f)[0].split(".")[0].split("_")[0]
        time_s = re.findall(r"T(.*)F", f)[0]
        z_ind = re.findall(r"Z(.*)C", f)[0]
        return [plate, well, time_s, chl, z_ind, site]

    def unique(self, list_):
        """
        Giving a list with double values,
        return a list with unique values.
        """

        unique_list = []

        for x in list_:
            if x not in unique_list:
                unique_list.append(x)
        ch = []
        for x in unique_list:
            ch.append(x)
        return ch

    def complete(self):
        """
        Method from base class, if return False
        luigi task is running, if True it ends.
        """
        if self.done:
            return True
        return False

    def output(self):

        """
        Method from base class to write logs of the
        luigi task
        """

        f_log = (
            f"./log/{self.task_name}_"
            + str(uuid.uuid4())[:8]
            + str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S"))
            + ".txt"
        )

        return luigi.LocalTarget(f_log)

    def do_proc(self, cmd, c):

        """
        This function is used when local scheduler is selected.
        It takes as input the cmd and a channel.
        It executes the task executable,
        passing the channel as parameter.

        :param cmd: command bash
        :type cmd: str
        :param c: channel
        :type c: str
        """

        process = Popen(cmd + [f"{c}"], stdout=PIPE, stderr=PIPE)

        stdout, stderr = process.communicate()
        if not stderr:
            debug("--No errors--\n", stdout.decode())
        else:
            debug("--Error--\n", stderr.decode())

        with self.output().open("w") as outfile:
            outfile.write(f"{stderr}\n")
        return process

    def run(self):
        """
        Method from base class. Here create the hierarchy of
        the zarr folder than launch the conversion task.
        """
        rows_cols = self.other_param["dims"]
        rows = rows_cols.split(",")[0]
        cols = rows_cols.split(",")[1]
        chl = []
        well = []
        plate = []
        for i in glob(self.in_path + "*." + self.ext):
            try:
                plate.append(self.metadata(os.path.basename(i))[0])
            except IndexError:
                pass

        plate_unique = self.unique(plate)

        if self.sclr == "local":
            run = ["python"]

            for plate in plate_unique:
                group_plate = zarr.group(self.out_path + f"{plate}.zarr")
                well = [
                    self.metadata(os.path.basename(fn))[1]
                    for fn in glob(self.in_path + f"{plate}_*." + self.ext)
                ]
                well_unique = self.unique(well)

                well_rows_columns = [
                    ind for ind in sorted([(n[0], n[1:]) for n in well_unique])
                ]

                group_plate.attrs["plate"] = {
                    "acquisitions": [
                        {"id": id_, "name": name}
                        for id_, name in enumerate(plate_unique)
                    ],
                    "columns": [
                        {"name": well_row_column[1]}
                        for well_row_column in well_rows_columns
                    ],
                    # takes unique rows from (row,col) tuples
                    "rows": [
                        {"name": u_row}
                        for u_row in set(
                            [
                                well_row_column[0]
                                for well_row_column in well_rows_columns
                            ]
                        )
                    ],
                    "wells": [
                        {
                            "path": well_row_column[0]
                            + "/"
                            + well_row_column[1],
                        }
                        for well_row_column in well_rows_columns
                    ],
                }

                for row, column in well_rows_columns:

                    group_well = group_plate.create_group(f"{row}/{column}/")

                    chl = [
                        self.metadata(os.path.basename(fn))[3]
                        for fn in glob(
                            self.in_path
                            + f"{plate}*_{row+column}*."
                            + self.ext
                        )
                    ]

                    chl_unique = self.unique(chl)
                    group_well.attrs["well"] = {
                        "images": [
                            {"path": f"{int(chl)-1}"} for chl in chl_unique
                        ],
                        "version": "0.3",
                    }

                    for ch in ["0"]:
                        group_field = group_well.create_group(
                            f"{int(ch)}/"
                        )  # noqa: F841

                        group_field.attrs["multiscales"] = [
                            {
                                "version": "0.3",
                                "axes": [
                                    {"name": "c", "type": "channel"},
                                    {
                                        "name": "z",
                                        "type": "space",
                                        "unit": "micrometer",
                                    },
                                    {"name": "y", "type": "space"},
                                    {"name": "x", "type": "space"},
                                ],
                                "datasets": [{"path": "0"}],
                                "metadata": {
                                    "kwargs": {"axes_names": ["y", "x"]}
                                },
                            }
                        ]

                        cmd = copy.copy(run)
                        cmd.extend(
                            [
                                self.tasks_path + self.task_name + ".py",
                                self.in_path,
                                self.out_path,
                                f"{plate}.zarr/" + f"{well}_{ch}/",
                                self.delete_in,
                                rows,
                                cols,
                                self.ext,
                            ]
                        )
                        process = Popen(
                            cmd + [f"{ch}"], stdout=PIPE, stderr=PIPE
                        )
                        stdout, stderr = process.communicate()
                        if not stderr:
                            print("--No errors--\n")
                        else:
                            print("--Error--\n", stderr.decode())

                        with self.output().open("w") as outfile:
                            outfile.write(f"{stderr}\n")

        elif self.sclr == "slurm":

            srun = ""
            # loop over plate, each plate could have n wells
            # debug(plate_unique)
            for plate in plate_unique:
                group_plate = zarr.group(self.out_path + f"{plate}.zarr")
                well = [
                    self.metadata(os.path.basename(fn))[1]
                    for fn in glob(self.in_path + f"{plate}_*." + self.ext)
                ]
                well_unique = self.unique(well)

                well_rows_columns = [
                    ind for ind in sorted([(n[0], n[1:]) for n in well_unique])
                ]

                group_plate.attrs["plate"] = {
                    "acquisitions": [
                        {"id": id_, "name": name}
                        for id_, name in enumerate(plate_unique)
                    ],
                    # takes unique cols from (row,col) tuples
                    "columns": sorted(
                        [
                            {"name": u_col}
                            for u_col in set(
                                [
                                    well_row_column[1]
                                    for well_row_column in well_rows_columns
                                ]
                            )
                        ],
                        key=lambda key: key["name"],
                    ),
                    # takes unique rows from (row,col) tuples
                    "rows": sorted(
                        [
                            {"name": u_row}
                            for u_row in set(
                                [
                                    well_row_column[0]
                                    for well_row_column in well_rows_columns
                                ]
                            )
                        ],
                        key=lambda key: key["name"],
                    ),
                    "wells": [
                        {
                            "path": well_row_column[0]
                            + "/"
                            + well_row_column[1],
                        }
                        for well_row_column in well_rows_columns
                    ],
                }
                # debug(well_rows_columns)
                for row, column in well_rows_columns:

                    group_well = group_plate.create_group(f"{row}/{column}/")

                    chl = [
                        self.metadata(os.path.basename(fn))[3]
                        for fn in glob(
                            self.in_path
                            + f"{plate}*_{row+column}*."
                            + self.ext
                        )
                    ]

                    chl_unique = self.unique(chl)
                    # debug(chl_unique)
                    group_well.attrs["well"] = {
                        "images": [
                            {
                                "path": "0"
                            }  # multiscale level, until pyramids just 0
                        ],
                        "version": "0.3",
                    }

                    cores = str(
                        self.slurm_param["cores"]
                    )  # "4" #slurm_param["cores"]
                    mem = str(
                        self.slurm_param["mem"]
                    )  # "1024" #slurm_param["mem"]
                    nodes = str(self.slurm_param["nodes"])

                    loader = jinja2.FileSystemLoader(searchpath="./")
                    env = jinja2.Environment(loader=loader)  # nosec
                    t = env.get_template("job.default.j2")
                    job = (
                        self.wf_name
                        + "_"
                        + self.task_name
                        + str(random.randrange(0, 101, 5))  # nosec
                    )

                    levels = ["0", "1", "2", "3", "4"]

                    group_field = group_well.create_group("0/")  # noqa: F841

                    group_field.attrs["multiscales"] = [
                        {
                            "version": "0.3",
                            "axes": [
                                {"name": "c", "type": "channel"},
                                {
                                    "name": "z",
                                    "type": "space",
                                    "unit": "micrometer",
                                },
                                {"name": "y", "type": "space"},
                                {"name": "x", "type": "space"},
                            ],
                            "datasets": [
                                {"path": f"{level}"} for level in levels
                            ],
                        }
                    ]

                srun += " ".join(
                    [
                        "python",
                        self.tasks_path + self.task_name + ".py ",
                        "-i " + self.in_path,
                        "-o " + self.out_path,
                        f"-z {plate}.zarr/" + "$RO/" + "$CO/0/",
                        "-d " + self.delete_in,
                        "-r " + rows,
                        "-c " + cols,
                        "-e " + self.ext,
                        "-C " + "${input[@]}",
                    ]
                )

                with open(f"./jobs/{job}", "w") as f:
                    f.write(
                        t.render(
                            job_name="test"
                            + str(random.randrange(0, 101, 5)),  # nosec
                            nodes=nodes,
                            cores=cores,
                            mem=mem + "MB",
                            array=len(well_unique) - 1,
                            channels=str(tuple(chl_unique)).replace(",", ""),
                            wells=str(
                                [
                                    unique_r
                                    for unique_r in set(
                                        [r for r in well_rows_columns]
                                    )
                                ]
                            )
                            .replace("'", "")
                            .replace("[", "")
                            .replace("]", "")
                            .replace(",", "")
                            .replace("(", '"')
                            .replace(")", '"'),
                            command=srun,
                        )
                    )

                cmd = ["sbatch", f"./jobs/{job}"]
                # debug(cmd)
                process = Popen(cmd, stderr=PIPE)
                stdout, stderr = process.communicate()
                if not stderr:
                    print("--No errors--\n")
                else:
                    print("--Error--\n", stderr.decode())

        self.done = True
