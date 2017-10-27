from sockjs.tornado import SockJSConnection, SockJSRouter
from tornado import web, ioloop, httpserver
from collections import OrderedDict

import threading
import sys
import json
import numpy as np

clients = set()
n_messages = 0
current_state = None

class PointsController:
    """Simple class to set & get points from neuroglancer
    """
    def __init__(self):
        self.val = 0 

    def get(self):
        global current_state
        if 'points' in current_state['layers']['points']:
            return current_state['layers']['points']['points']
        else:
            return None

    def set(self, points):
        global current_state
        pts = [points.tolist()]
        current_state['layers']['points'] = {'type':'point', 'points':pts}
        broadcast()

# websockets connections
class Connection(SockJSConnection):
    def on_open(self, info):
        """
        info is an object which contains caller IP address, query string
        parameters and cookies associated with this request"""
        # When new client comes in, will add it to the clients list
        clients.add(self)

    def on_message(self, json_state):
        """
        This will call initialize_current_state or on_current_state_change depening on if it is
        the first message recieved.
        """
        global current_state
        print(json_state)
        current_state = json.JSONDecoder(object_pairs_hook=OrderedDict).decode(json_state)
        global n_messages
        if not n_messages:
            print('current_state initialized')
        n_messages += 1

    def on_close(self):
        # If client disconnects, remove him from the clients list
        clients.remove(self)

    def initialize_current_state(self, current_state):
        """
        This is called once the connection is stablished
        """
        return current_state

    def on_current_state_change(self, current_state):
        """
        This is called every time there is a new current_state available
        (except the very first time).
        """
        return current_state


# Tornado & Tk need to run on separate threads
class TornadoThread(threading.Thread):
    def __init__(self):
        # super(TornadoThread, self).__init__()
        # self._stop = threading.Event()
        threading.Thread.__init__(self)
        self.daemon = True

        socketApp = web.Application(router.urls)
        http_server = httpserver.HTTPServer(socketApp, ssl_options={
            "certfile": "./certificate.crt",
            "keyfile": "./privateKey.key",
        })
        http_server.bind(9998) #port
        http_server.start(1)

    def run(self):
        print("IOLoop starting")
        ioloop.IOLoop.instance().start()


router = SockJSRouter(Connection)
# In order for the webbrowser to connect to this server
# add to the url 'current_stateURL':'http://localhost:9999'
def broadcast():
    """
    Use this method to broadcast a new current_state to all connected clients.
    Without the need to wait for an `on_current_state_change`.
    """
    global current_state
    router.broadcast(clients.copy(), json.dumps(current_state))

thread = TornadoThread()
thread.start()
