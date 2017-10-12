import sys
from os.path import expanduser, join, isfile
import threading
from tornado import web, ioloop, httpserver
from sockjs.tornado import SockJSConnection, SockJSRouter
import json
from collections import OrderedDict
import numpy as np

clients = set()
n_messages = 0

def read_points(fn):
    if isfile(fn):
        print 'read_points: ' + fn
        return np.genfromtxt(fn, delimiter=",", dtype=int).tolist()
    else:
        return []

def write_points(fn, points):
    if points:
        print 'write_points: ' + fn
        np.savetxt(fn, points, delimiter=',', fmt='%d')

def get_points(state):
    if 'points' in state['layers']['points']:
        points = state['layers']['points']['points']
        return [[int(round(pt)) for pt in point] for point in points]
    else:
        return None

def get_z(state):
    return int(state['navigation']['pose']['position']['voxelCoordinates'][2])

def get_filename(z):
    return str(z) + ".csv"

def set_points(state, points):
    state['layers']['points'] = {'type':'point', 'points':points}
    broadcast(state)

current_z = 0

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
        state = json.JSONDecoder(object_pairs_hook=OrderedDict).decode(json_state)
        global n_messages

        if not n_messages: #first message ever
            self.initialize_state(state)
        else:
            self.on_state_change(state)
        n_messages += 1


    def on_close(self):
        # If client disconnects, remove him from the clients list
        clients.remove(self)

    def initialize_state(self, state):
        """
        This is called once the connection is stablished
        """
        current_z = get_z(state)
        filename = get_filename(current_z)
        points = read_points(filename)
        set_points(state, points)

    def on_state_change(self, state):
        """
        This is called every time there is a new state available
        (except the very first time).
        """
        # store position
        global current_z
        points = get_points(state)
        z = get_z(state)
        if (current_z != 0) & (z != current_z):
            print str(current_z) + ' -> ' + str(z)
            filename = get_filename(current_z)
            write_points(filename, points)
            filename = get_filename(z)
            points = read_points(filename)
            set_points(state, points)
        current_z = z

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
