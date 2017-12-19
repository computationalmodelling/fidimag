from __future__ import print_function
import os
import time
import numpy as np
from multiprocessing import Process, Queue, Lock

lock = Lock()


class TaskState(object):
    """
    Each task has three states: "Done!", "Started!", "Unstarted!"
    """

    def __init__(self, taskname):
        self.taskname = taskname
        self.state = {}
        self.load()

    def load(self):
        if not os.path.exists(self.taskname):
            return
        f = open(self.taskname, 'r')
        data = f.read()
        f.close()

        for line in data.splitlines():
            k, v = line.split(':')
            self.state[k.strip()] = v.strip()

    def save_state(self):
        f = open(self.taskname, 'w')
        for (k, v) in self.state.items():
            f.write(u'%s : %s\n' % (k, v))
        f.close()

    def update_state(self, k, v, save=True):
        """
        v = 0: Done
        v = 1: Started but Unfinished
        v = 2: Unstarted
        """
        key = self.dict2str(k)

        if save:
            self.load()

        if v==0:
            self.state[key] = 'Done!'
        elif v==1:
            self.state[key] = 'Started!'
        else:
            self.state[key] = 'Unstarted!'

        if save:
            self.save_state()

    def dict2str(self, d):
        res = []
        for k in d:
            res.append(k)
            res.append(str(d[k]))
        return '_'.join(res)


    def is_unstarted(self, k):
        key = self.dict2str(k)
        if key not in self.state:
            self.update_state(k, 2, False)
            return True

        if key in self.state:
            if not 'Done' in self.state[key] and not 'Started' in self.state[key]:
                self.update_state(k, 2, False)
                return True
        return False

    def is_done(self, k):
        key = self.dict2str(k)
        if key not in self.state:
            self.update_state(k, 2, True)
            return False
        if key in self.state:
            if 'Done' in self.state[key]:
                return True
        return False

class BatchTasks(object):

    def __init__(self, fun, processes=4, taskname='task', waiting_time=1):
        self.fun = fun
        self.tasks = [{}]
        self.parameters = []
        self.current_directory = os.getcwd()

        self.ts = TaskState(taskname + '.txt')
        self.waiting_time = waiting_time
        self.dims = []

        self.processes = processes

        self.process_res = []

    def add_parameters(self, name, values):
        new_tasks = []
        self.parameters.append(name)
        self.dims.append(len(values))

        for task in self.tasks:
            for v in values:
                t = dict(task)
                t[name] = v
                new_tasks.append(t)

        self.tasks = list(new_tasks)

    def generate_directory(self, task):
        base = self.current_directory
        for name in self.parameters:
            base = os.path.join(base, name + '_' + str(task[name]))

        return base

    def run_single(self):

        while not self.task_q.empty():
            task = self.task_q.get()

            dirname = self.generate_directory(task)
            if not os.path.exists(dirname):
                os.makedirs(dirname)

            lock.acquire()
            self.ts.update_state(task, 1)
            lock.release()

            os.chdir(dirname)
            self.fun(**task)
            os.chdir(self.current_directory)

            lock.acquire()
            self.ts.update_state(task, 0)
            lock.release()

            time.sleep(self.waiting_time)

    def start(self, run_started=False):

        self.task_q = Queue()
        for task in self.tasks:
            if run_started:
                if not self.ts.is_done(task):
                    self.task_q.put(task)
            else:
                if self.ts.is_unstarted(task):
                    self.task_q.put(task)

        self.ts.save_state()

        self.threads = []

        for _ in range(self.processes):
            t = Process(target=self.run_single)
            t.start()
            self.threads.append(t)
            time.sleep(self.waiting_time)

        for t in self.threads:
            t.join()

    def post_process(self, fun):
        for task in self.tasks:
            dirname = self.generate_directory(task)
            os.chdir(dirname)
            try:
                self.process_res.append(fun(**task))
            except RuntimeError:
                print('error:', task)
            os.chdir(self.current_directory)

    def get_res(self, key=None, value=None):
        res = []
        par = []

        if len(self.parameters) == 1:
            for i, task in enumerate(self.tasks):
                par.append(task.values()[0])
                res.append(self.process_res[i])
        elif len(self.parameters) == 2:
            pass
            # I can not remember what's the meaning of the following lines.
            """
            for i,task in enumerate(self.tasks):
                if key in task and task[key]==value:
                    res.append(self.process_res[i])
                    tmp_task=dict(task)
                    del tmp_task[key]
                    par.append(tmp_task.values()[0])
            """
        else:
            raise NotImplementedError(
                'Only support one- and two- parameter case!')

        if key is None and value is None and len(self.parameters) == 2:
            v0 = self.parameters[0]
            v1 = self.parameters[1]
            for i, task in enumerate(self.tasks):
                par.append((task[v0], task[v1]))
                res.append(self.process_res[i])

        return np.array(par), np.array(res)


def task(p1, p2):
    print('current directory:', os.getcwd())
    res = 'p1=' + str(p1) + '  p2=' + str(p2)
    #res= 1/0
    with open('res.txt', 'w') as f:
        f.write(res)
    time.sleep(1)


if __name__ == "__main__":
    tasks = BatchTasks(task, 2)
    tasks.add_parameters('p1', ['a', 'b', 'c'])
    tasks.add_parameters('p2', range(1, 5))
    tasks.start()
