import json
import sys

import luigi
from wrappers.compressionTaskWrap import CompressionTaskWrap
from wrappers.conversionTaskWrap import ConversionTaskWrap

# TODO add task to clean old logs.
# class CleanLogs(luigi.Task):
#     pass


DICT_TASK = {
    "compression_tif": CompressionTaskWrap,
    "yokogawa_to_zarr": ConversionTaskWrap,
}


class WorkflowTask(luigi.Task):

    """
    WorkflowTask class takes a dictionary as input in which
    are stored all the inputs for tasks.
    It is an extension of the luigi.Task class.
    """

    flags = luigi.parameter.DictParameter()

    _complete = False

    def run(self):

        """
        Method from base class. Here return a generator of tasks.
        """

        task_names = self.flags["tasks"].keys()
        in_path = self.flags["arguments"]["resource_in"]
        wf_name = self.flags["arguments"]["workflow_name"]

        out_path = self.flags["arguments"]["resource_out"]
        delete_in = self.flags["arguments"]["delete"]
        sclr = self.flags["arguments"]["scheduler"]
        ext = self.flags["arguments"]["ext"]
        other_params = self.flags["arguments"]["other_params"]
        slurm_params = self.flags["arguments"]["slurm_params"]

        for task_name in task_names:
            yield DICT_TASK[task_name](
                task_name=task_name,
                wf_name=wf_name,
                in_path=in_path,
                out_path=out_path,
                delete_in=delete_in,
                sclr=sclr,
                ext=ext,
                other_param=other_params,
                slurm_param=slurm_params,
            )

        self._complete = True

    def complete(self):

        """
        Method from base class to check if run method has finished.
        """

        return self._complete


if __name__ == "__main__":

    luigi.build(
        [WorkflowTask(flags=json.loads(sys.argv[1]))],
        workers=1,
        parallel_scheduling=True,
    )
