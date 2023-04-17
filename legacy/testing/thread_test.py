import threading

import logging
logger = logging.getLogger(__name__)

class MyThread(threading.Thread):
    def __init__(self, threadId, name, test_list):
        threading.Thread.__init__(self)
        self.threadId = threadId
        self.name = name
        self.test_list = test_list
    
    def run(self):
        logger.debug(f"Starting {self.name}")
        thread_function(self.test_list)
        logger.debug(f"Exiting {self.name}")
    
def thread_function(test_list):
    logger.debug(f"test_list = {test_list}")

def main():
    # Create new threads
    thread1 = MyThread(1, "Thread-1", ["1","2","3"])
    thread2 = MyThread(2, "Thread-2", ["4", "5", "6"])

    # Start new threads
    thread1.start()
    thread2.start()

    logger.debug("Exiting main thread.")

if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()