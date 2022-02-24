import datetime
import json
import os
import re
import sys
import uuid
from functools import partial
from glob import glob
from multiprocessing import Pool
from subprocess import PIPE  # nosec
from subprocess import Popen  # nosec

import luigi
import zarr
from devtools import debug

# TODO add task to clean old logs.
# class CleanLogs(luigi.Task):
#     pass

# TODO handle in-memory target

# _data = {}

# class MemoryTarget(luigi.Target):
#     _data = {}
#     def __init__(self, path):
#         self.path = path
#     def exists(self):
#         return self.path in _data
#     def put(self, value):
#         _data[self.path] = value
#     def get(self):
#         return _data[self.path]


class CompressionTaskWrap(luigi.Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.complete_flag = False

    # accepts_messages = True

    task_name = luigi.parameter.Parameter(significant=True)

    in_path = luigi.parameter.Parameter(significant=True)

    out_path = luigi.parameter.Parameter(significant=True)

    delete_in = luigi.parameter.Parameter()

    sclr = luigi.parameter.Parameter()

    other_param = luigi.parameter.DictParameter()

    tasks_path = os.getcwd() + "/tasks/"

    done = False

    def complete(self):
        if self.done:
            return True
        return False

    def output(self):

        f_log = (
            f"./log/{self.task_name}_"
            + str(uuid.uuid4())[:8]
            + str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S"))
            + ".txt"
        )

        return luigi.LocalTarget(f_log)

    def do_comp(self, cmd, interval):

        process = Popen(  # nosec
            cmd + [f"{interval[0]}"] + [f"{interval[1]}"],  # nosec
            stderr=PIPE,  # nosec
        )  # nosec

        stdout, stderr = process.communicate()
        debug(process.communicate())
        with self.output().open("w") as outfile:
            outfile.write(f"{stderr}\n")
        return process

    def run(self):

        self.done = True

        if self.sclr == "local":
            cmd = ["python"]
        else:
            cmd = ["srun", "python"]

        batch_size = self.other_param["batch_size"]
        l_file = glob(self.in_path + "*.tif")

        batch = len(l_file) // batch_size

        tmp_s = 0
        intervals = []
        for tmp_e in range(batch, len(l_file) + 1, batch):
            intervals.append((tmp_s, tmp_e))
            tmp_s = tmp_e

        if self.sclr == "local":
            cmd = ["python"]
        else:
            cmd = ["srun", "python"]

        cmd.extend(
            [
                self.tasks_path + self.task_name + ".py",
                self.in_path,
                self.out_path,
                self.delete_in,
            ]
        )

        debug(intervals)
        p = Pool()
        do_comp_part = partial(self.do_comp, cmd)
        p.map_async(do_comp_part, intervals)
        p.close()
        p.join()

        debug(cmd)

        self.done = True


class ConversionTaskWrap(luigi.Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.complete_flag = False

    retry_count = 5
    # accepts_messages = True

    task_name = luigi.parameter.Parameter(significant=True)

    in_path = luigi.parameter.Parameter(significant=True)

    out_path = luigi.parameter.Parameter(significant=True)

    delete_in = luigi.parameter.Parameter()

    sclr = luigi.parameter.Parameter()

    other_param = luigi.parameter.DictParameter()

    tasks_path = os.getcwd() + "/tasks/"

    done = False

    def metadata(self, filename):
        ##
        f = filename.rsplit(".", 1)[0]
        f_s = f.split("_")
        plate = f_s[0]
        well = f_s[1]
        site = re.findall(r"F(.*)L", f_s[2])[0]
        chl = re.findall(r"C(.*)", f_s[2])[0]
        return [well, site, chl, plate]

    def unique(self, list1):

        unique_list = []

        for x in list1:
            if x not in unique_list:
                unique_list.append(x)
        ch = []
        for x in unique_list:
            ch.append(x)
        return ch

    def complete(self):
        if self.done:
            return True
        return False

    def output(self):

        f_log = (
            f"./log/{self.task_name}_"
            + str(uuid.uuid4())[:8]
            + str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S"))
            + ".txt"
        )

        return luigi.LocalTarget(f_log)

    def do_proc(self, cmd, c):

        process = Popen(cmd + [f"{c}"], stdout=PIPE, stderr=PIPE)

        stdout, stderr = process.communicate()
        debug(process.communicate())
        with self.output().open("w") as outfile:
            outfile.write(f"{stderr}\n")
        return process

    def run(self):

        chl = []
        well = []
        plate = []
        for i in glob(self.in_path + "*.tif"):
            chl.append(self.metadata(os.path.basename(i))[2])
            well.append(self.metadata(os.path.basename(i))[0])
            plate.append(self.metadata(os.path.basename(i))[3])

        chl_unique = self.unique(chl)
        well_unique = self.unique(well)
        plate_unique = self.unique(plate)
        if self.sclr == "local":
            cmd = ["python"]
        else:
            cmd = ["srun", "python"]

        # at the moment we consider just one well.
        # then it will be necessary iterate over plate/s
        # and well/s (2 for loops)

        # for loop plates
        group_plate = zarr.group(self.out_path + f"{plate_unique[0]}.zarr")
        group_plate.attrs["plate"] = {
            "acquisition": [
                {"id": id_, "name": name}
                for id_, name in enumerate(plate_unique)
            ],
            "columns": [],
            "rows": [],
            "wells": [
                {"path": well_unique[0], "rowIndex": None, "columnIndex": None}
            ],
        }

        # for loop wells and each well have n channels
        group_well = group_plate.create_group(f"{well_unique[0]}")
        group_well.attrs["well"] = {
            "images": [
                {"acquisition": 0, "path": path}  # id_ of plate
                for path in chl_unique
            ]
        }

        cmd.extend(
            [
                self.tasks_path + self.task_name + ".py",
                self.in_path,
                self.out_path
                + f"{plate_unique[0]}.zarr/"
                + f"{well_unique[0]}",
                self.delete_in,
            ]
        )

        p = Pool()
        do_proc_part = partial(self.do_proc, cmd)
        p.map_async(do_proc_part, chl_unique)
        p.close()
        p.join()

        debug(cmd)

        self.done = True


DICT_TASK = {
    "compression_tif": CompressionTaskWrap,
    "tif_to_zarr": ConversionTaskWrap,
}


class WorkflowTask(luigi.Task):

    flags = luigi.parameter.DictParameter()

    _complete = False

    def run(self):

        task_names = self.flags["tasks"].keys()
        in_paths = [r_in for r_in in self.flags["arguments"]["resource_in"]]
        out_paths = [
            r_out for r_out in self.flags["arguments"]["resource_out"]
        ]
        delete_ins = [d for d in self.flags["arguments"]["delete"]]
        sclr = self.flags["arguments"]["scheduler"]
        other_params = self.flags["arguments"]["other_params"]

        for i, task_name in enumerate(task_names):
            yield DICT_TASK[task_name](
                task_name=task_name,
                in_path=in_paths[i],
                out_path=out_paths[i],
                delete_in=delete_ins[i],
                sclr=sclr,
                other_param=other_params,
            )

        self._complete = True

    def complete(self):
        return self._complete


if __name__ == "__main__":

    luigi.build(
        [WorkflowTask(flags=json.loads(sys.argv[1]))],
        workers=1,
        parallel_scheduling=True,
    )
