# Showcases Timer.py functionality
# TODO: time measurement
# DONE: time comparison

from goreverselookuplib.Timer import Timer

timestamp_one = "2023-06-29 15:30:00"
timestamp_two = "2023-06-29 16:00:00"

result = Timer.compare_time(timestamp_one, timestamp_two)
print(result)  # Will print True if timestamp_two > timestamp_one