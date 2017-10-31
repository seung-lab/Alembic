from sockjs.tornado import SockJSConnection, SockJSRouter
from tornado import web, ioloop, httpserver
from collections import OrderedDict

import threading
import sys
import json
import numpy as np
import os

class Controller:
    """Simple class to set & get points from neuroglancer
    """
    state = {}
    clients = set()

    def __init__(self, tornado_thread):
        self.router = tornado_thread.router
        # if self.get() is None:
        #     self.set(np.array([]))

    def get(self):
        if 'points' in self.state['layers']:
            if 'points' in self.state['layers']['points']:
                return self.state['layers']['points']['points']
        return None

    def set(self, points):
        pts = [points.tolist()]
        self.state['layers']['points'] = {'type':'point', 'points':pts}
        # import pdb; pdb.set_trace()
        self.router.broadcast(self.clients.copy(), json.dumps(self.state))

    # websockets connections
    class Connection(SockJSConnection):

        def on_open(self, info):
            """
            info is an object which contains caller IP address, query string
            parameters and cookies associated with this request"""
            # Whcren new client comes in, will add it to the clients list
            Controller.clients.add(self)

        def on_message(self, json_state):
            """
            This will call initialize_current_state or on_current_state_change depening on if it is
            the first message recieved.
            """
            print(json_state)
            Controller.state = json.JSONDecoder(object_pairs_hook=OrderedDict).decode(json_state)

        def on_close(self):
            # If client disconnects, remove him from the clients list
            self.clients.remove(self)


# Tornado & Tk need to run on separate threads
class TornadoThread(threading.Thread):
    def __init__(self, port):
        # super(TornadoThread, self).__init__()
        # self._stop = threading.Event()
        threading.Thread.__init__(self)
        self.daemon = True

        self.router = SockJSRouter(Controller.Connection)
        socketApp = web.Application(self.router.urls)
        self.http_server = httpserver.HTTPServer(socketApp, ssl_options={
            "certfile": os.path.join(sys.path[0], 'certificate.crt'),
            "keyfile": os.path.join(sys.path[0], 'privateKey.key')
        })
        self.http_server.bind(port) #port
        self.http_server.start(1)

    def run(self):
        print("IOLoop starting")
        ioloop.IOLoop.instance().start()

def create(port=9999):
    thread = TornadoThread(port)
    thread.start()
    return Controller(thread.router)

if __name__ == '__main__':
    create()