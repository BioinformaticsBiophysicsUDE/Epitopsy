'''
Created on Jul 9, 2011

@author: chris
'''

import multiprocessing

class ConsumeRotation(multiprocessing.Process):
    '''
    This class is derived from multiprocessing.Process.
    I use this class for the 'processRotation'-function
    of the 'ParallelCombinationRotator'-class.
    '''
    def __init__(self, task_queue, result_queue):
        multiprocessing.Process.__init__()
        self.task_queue = task_queue
        self.result_queue = result_queue
    def run(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # None is a "Poison pill" and means shutdown
                self.task_queue.task_done()
                break
            answer = next_task()
            self.task_queue.task_done()
            self.result_queue.put(answer)
        return