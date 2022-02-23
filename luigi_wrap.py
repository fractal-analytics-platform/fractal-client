import os
import sys
from subprocess import PIPE  # nosec
from subprocess import Popen
from devtools import debug
import luigi
import datetime

#TODO add task to clean old logs.
# class CleanLogs(luigi.Task):
#     pass

#TODO handle in memory target

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



class TaskWrap(luigi.Task):

    # accepts_messages = True 

    task_name = luigi.parameter.Parameter(significant=True)

    in_path =  luigi.parameter.Parameter(significant=True)

    out_path = luigi.parameter.Parameter(significant=True)

    delete_in = luigi.parameter.Parameter()

    tasks_path = os.getcwd() + "/tasks/"

    done = False

    def output(self):
        
        f_log = f'./log/{self.task_name}_' + str(datetime.datetime.now().strftime(
                                                '%Y_%m_%d_%H_%M_%S')) + '.txt'

        return luigi.LocalTarget(f_log)


    def requires(self):
        return None
   

    def complete(self):
        if self.done:
            return True
        return False


    def run(self):

        process = Popen(
            ["python",  self.tasks_path + self.task_name + ".py", self.in_path, 
             self.out_path, self.delete_in],
            stdout=PIPE,
            stderr=PIPE,
        )
        
        self.done = True
        
        stdout, stderr = process.communicate()

        with self.output().open("w") as outfile:
            outfile.write(f"{stderr}\n")
        


class WorkflowTask(luigi.WrapperTask):

    flags = luigi.parameter.DictParameter()

    def requires(self):

    #TODO add here task dependencies 

        task_names = self.flags['tasks'].keys()
        in_paths = [r_in for r_in in self.flags["arguments"]["resource_in"]]
        out_paths = [r_out for r_out in self.flags["arguments"]["resource_out"]]
        delete_ins = [d for d in self.flags["arguments"]["delete"]]
        for i, task_name in enumerate(task_names):
             yield TaskWrap(task_name=task_name,in_path=in_paths[i],
                                out_path=out_paths[i],delete_in=delete_ins[i])




if __name__ == '__main__':
    luigi.run(['WorkflowTask', '--workers',
     '1','--flags', sys.argv[1], '--local-scheduler'])
