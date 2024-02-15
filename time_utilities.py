import time
import datetime

def posix_sec_now():
    return time.time()

def datetime_now():
    return posix_sec_to_datetime(posix_sec_now())

def datetime_now_label():
    return str(datetime_now()).split('+')[0].replace(' ', '_').replace('-', '_').replace(':', '_').replace('.', '_').replace('+', '_')

def posix_sec_to_datetime(posix_sec):
    return datetime.datetime.fromtimestamp(float(posix_sec)).replace(tzinfo=datetime.timezone.utc)

def posix_msec_to_datetime(posix_msec):
    return posix_sec_to_datetime(float(posix_msec) / 1000)
