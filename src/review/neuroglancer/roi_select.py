import sys
from os.path import expanduser, join
import threading
from tornado import web, ioloop, httpserver
from sockjs.tornado import SockJSConnection, SockJSRouter
import json
from collections import OrderedDict
import numpy as np

clients = set()
n_messages = 0

def read_points(fn):
    return np.genfromtxt(fn, delimiter=",", dtype=int).tolist()

def write_points(fn, points):
    print 'write_points: ' + str(points)
    np.savetxt(fn, points, delimiter=',', fmt='%d')

def get_points(state):
    print 'get_points'
    if 'points' in state['layers']['points']:
        points = state['layers']['points']['points']
        return [[int(round(pt)) for pt in point] for point in points]
    else:
        return None

filename = "temp.csv"
current_points = []
if len(sys.argv) > 1:
    filename = sys.argv[1]
    current_points = read_points(filename)

class Connection(SockJSConnection):
    def on_open(self, info):
        """
        info is an object which contains caller IP address, query string
        parameters and cookies associated with this request"""
        # When new client comes in, will add it to the clients list
        clients.add(self)

    def on_message(self, json_state):
        """
        This will call initialize_state or on_state_change depening on if it is
        the first message recieved.
        """
        current_state = json.JSONDecoder(object_pairs_hook=OrderedDict).decode(json_state)
        global n_messages

        if not n_messages: #first message ever
            self.initialize_state(current_state)
        else:
            self.on_state_change(current_state)
        n_messages += 1


    def on_close(self):
        # If client disconnects, remove him from the clients list
        clients.remove(self)

    def initialize_state(self, state):
        """
        This is called once the connection is stablished
        """
        state['layers']['points'] = {'type':'point', 'points':current_points}
        broadcast(state)

    def on_state_change(self, state):
        """
        This is called every time there is a new state available
        (except the very first time).
        """
        # store position
        global current_points
        points = get_points(state)
        if not points == current_points:
            current_points = points
            write_points(filename, points)

# In order for the webbrowser to connect to this server
# add to the url 'stateURL':'http://localhost:9999'
router = SockJSRouter(Connection)
def broadcast(state):
    """
    Use this method to broadcast a new state to all connected clients.
    Without the need to wait for an `on_state_change`.
    """
    router.broadcast(clients.copy(), json.dumps(state))

socketApp = web.Application(router.urls)
http_server = httpserver.HTTPServer(socketApp, ssl_options={
    "certfile": "./certificate.crt",
    "keyfile": "./privateKey.key",
})
http_server.bind(9999) #port
http_server.start(1)
try:
    ioloop.IOLoop.instance().start()
except KeyboardInterrupt:
    ioloop.IOLoop.instance().stop()