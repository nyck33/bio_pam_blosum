import os
from Bio import pairwise2
from Bio import Align


import redis
from rq import Worker, Queue, Connection

listen = ['high', 'default', 'low']

redis_url = os.getenv('REDISTOGO_URL', 'redis://localhost:6379')

conn2 = redis.from_url(redis_url)

if __name__ == '__main__':
    with Connection(conn2):
        worker = Worker(map(Queue, listen))
        worker.work()