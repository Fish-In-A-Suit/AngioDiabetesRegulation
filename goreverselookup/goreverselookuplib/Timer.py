import time
from datetime import timedelta
import logging
logger = logging.getLogger(__name__)

class Timer:
    def __init__(self):
        self.start_time = time.time()
    
    def set_start_time(self):
        """
        Sets a new reference start time.
        """
        self.start_time = time.time()
    
    def get_elapsed_seconds(self) -> int:
        """
        Returns the amount of seconds unformatted (contains decimal places)
        """
        return time.time() - self.start_time
    
    def get_elapsed_time(self) -> str:
        """
        Gets elapsed time in hh mm ss format.
        """
        sec = int(self.get_elapsed_seconds())
        td = timedelta(seconds=sec)
        return str(td)
    
    def print_elapsed_time(self, useLogger: bool = True, prefix: str = "Elapsed: "):
        """
        Prints the elapsed time in hh mm ss format. 
        
        Args:
          - useLogger: if True, then logger.info is used. If false, then print is used.
          - prefix: the string you want to use as a prefix
        """
        if useLogger:
            logger.info(f"{prefix}{self.get_elapsed_time()}")
        else:
            print(f"{prefix}{self.get_elapsed_time()}")