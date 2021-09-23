import redis


class Jobqueue:
    def __init__(self):
        self.redis_conn = redis.from_url("redis://localhost:6379")
        self.queue = "queue"

    def __len__(self):
        return self.redis_conn.llen(self.queue)

    def add(self, subID):
        return self.redis_conn.rpush(self.queue, subID)

    def top(self):
        elem = self.redis_conn.lindex(self.queue, 0)
        if elem:
            elem = elem.decode()
        return elem

    def pop(self):
        self.redis_conn.lpop("queue")

    def in_queue(self, subID):
        for elem in self.redis_conn.lrange(self.queue, 0, -1):
            if elem == subID:
                return True
        return False

    def all(self):
        return [i.decode("utf-8") for i in self.redis_conn.lrange(self.queue, 0, -1)]
